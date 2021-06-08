#---------------------------------------------------------------------------------------------------
# Code for running a monte carlo markov using a model structure and parameters specified in an excel
# document. (See README and test_parameters.xlsx for addional explaination and example)
#
# Brandon Chan | January 2021
#---------------------------------------------------------------------------------------------------
# Import packages/libraries
#---------------------------------------------------------------------------------------------------
import pandas as pd # Dataframe structure and manipulation
import numpy as np # Scientific computing functionality (array and matrix operations and some stats things)
import math # For additional math functions
from time import perf_counter, strftime # For timing the runtime of the code
from arch.bootstrap import IIDBootstrap # Bootstrap analysis for ICER CI calculation
# Importing required functions outlined in model_helpers.py
from model_helpers import check_excel_file, get_gamma, get_beta, calculate_residual, set_transition, check_row_sums, normalize_transitions, check_model_population 

class ModelSpec:
    '''
    A container to load in an excel sheet specification of a model and store parameters. 
    Performs a check on the validitiy of the specified model via the check_excel_file() function.
    Is used as an input for various functions of the analysis framework so that relevant model
    information is used. Also allows for in-code changes to model parameters if desired.

    Initialization args: filepath = full filepath to an excel sheet specification of a model arm
                         model_name = string that denotes the name of the model

    TODO:// add model specification checks when in-code modifications are made to parameters 
    '''
    def __init__(self, filepath, model_name):
        self.model_name = model_name
        self.excel_file = filepath

        input_file = pd.ExcelFile(filepath) # read in excel file
        check_excel_file(input_file)

        self.structure = pd.read_excel(input_file, 'transitions') 
        self.cost_specs = pd.read_excel(input_file, 'costs')
        self.util_specs = pd.read_excel(input_file, 'utilities')
        self.simulation_parameters = pd.read_excel(input_file, 'specification', header=None, index_col=0)

        # Specification of variables regarding states in model
        unique_states = self.structure['start_state'].unique().tolist()
        self.state_mapping = {i : unique_states.index(i) for i in unique_states}
        self.num_states = len(unique_states)

        # Initialize model parameters based on spreadsheet defined values
        self.max_iterations = int(self.simulation_parameters.loc['max_iterations'].values[0])
        self.cycle_length = self.simulation_parameters.loc['cycle_length'].values[0]
        self.time_horizon_days = self.simulation_parameters.loc['time_horizon'].values[0] * 365
        self.num_cycles = int(self.time_horizon_days / self.cycle_length)
        self.name_start_state = self.simulation_parameters.loc['name_start_state'].values[0]
        self.discount_rate = self.simulation_parameters.loc['discount_rate'].values[0]

    # TODO:// Consider adding mode functions like saving, etc.


class ModelData:
    '''
    A container to store the population, cost, and utility outputs for a single arm of a model.
    Is used as an input for convenience functions to calculate metrics of interest.

    On initialization, requires 3 numpy arrays (population, cost, and utility) outputs to be passed in
    
    Args: pdata = numpy array of shape [iteration_number x states x timesteps] representing the state to state movement in the model
          cdata = numpy array of shape [iteration_number x states x timesteps] representing the cost per state per cycle 
          udata = numpy array of shape [iteration_number x states x timesteps] representing the utility per state per cycle
    '''
    def __init__(self, pdata, cdata, udata):
        # Checks that all data have same dims (in two dimensions...)
        for i in range(0, 2):
            if pdata.shape[i] != cdata.shape[i] or pdata.shape[i] != udata.shape[i] or cdata.shape[i] != udata.shape[i]:
                raise ValueError('Error: Dimensions of all 3 input arrays must be the same')

        self.pop_data = pdata
        self.cost_data = cdata
        self.util_data = udata

        self.num_states = pdata.shape[0]
        self.num_cycles = pdata.shape[1]
        self.num_iterations = pdata.shape[2]

        # Auto-consense data to per-state and per-iteration
        self.cycle_cost_data = np.sum(self.cost_data, axis=1)
        self.cycle_util_data = np.sum(self.util_data, axis=1)

        self.iteration_cost_data = np.sum(self.cycle_cost_data, axis=1)
        self.iteration_util_data = np.sum(self.cycle_util_data, axis=1)        


def run_model(model_specification: ModelSpec):
    '''
    Function that takes a loaded model specification and runs the model based on the stored parameters

    Input: model_specification = a ModelSpec object with a loaded model specification
    Output: results_log = a 3D numpy array of the dimensions [iteration_number x states x timesteps] where timesteps = number of cycles
                          this represnets the "movement" throughout the model states at each cycle for each iteration.
    '''
    if isinstance(model_specification, ModelSpec) == False:
        raise TypeError('Error: Expected input model_specification to be of type markov_modeling.ModelSpec')
    
    #---------------------------------------------------------------------------------------------------
    # Parameter initialization
    #---------------------------------------------------------------------------------------------------
    transitions_df = model_specification.structure
    num_states = model_specification.num_states
    max_iterations = model_specification.max_iterations
    cycle_length = model_specification.cycle_length
    num_cycles = model_specification.num_cycles
    name_start_state = model_specification.name_start_state
    state_mapping = model_specification.state_mapping

    #---------------------------------------------------------------------------------------------------
    # Generate matrix representation of model and identify/log state transitions that require resampling
    # or additional calcuations when simulating
    # 
    # Dimensions of matrix = [num_states x num_states]
    # Rows map to starting state, columns map to target state
    #---------------------------------------------------------------------------------------------------
    # specify empty transition matrix
    transition_matrix = np.zeros((num_states, num_states)) 

    # Keep track of specific indicies in the transition matrix that need to be updated "in-simulation"
    # ie. time-dependent transitions and transitions that get resampled/recalculated every iteration
    # Intended to be stored as a list of dictionaries. Dictionaries contain key-value pairs that denote
    # type of transtion, index (i,j) of matrix corresponding to transtion, and appropriate parameters ie. a, b, shape, scale, etc.
    resample_indicies = []
    resample_r2p_indicies = []
    time_dependent_indicies = []
    probabilistic_time_dependent_indicies = []
    residual_indicies = []

    # Iterate through specified transtions and initialize constant values of matrix. Transitions that vary 
    # (ie. time dependent or resampled) are assigned in model iteration step. 
    # Log indicies of transitions that need to be updated to quickly index the correct position in the transition matrix
    # Does not do any initialization of Dirichlet transitions. That is handled separately
    # below
    for t in transitions_df.loc[transitions_df.type != 'dirichlet'].itertuples():
        start_state_index = state_mapping[t[1]] # mapped row number of start state
        end_state_index = state_mapping[t[2]] # mapped column number of end state
        t_type = t[3] # type of transition 
        params = t[4:] # parameters
        
        if t_type == 'constant':
            transition_matrix[start_state_index, end_state_index] = set_transition('constant', transition=params[0])
        elif t_type in ['beta', 'gamma']:
            resample_indicies += [{'start_state':t[1], 'end_state':t[2], 
                                   'i':start_state_index, 'j':end_state_index, 
                                   'type':t_type, 
                                   'a':params[0], 'b':params[1]}]
        elif t_type in ['beta_r2p']:
            resample_r2p_indicies += [{'start_state':t[1], 'end_state':t[2], 
                                       'i':start_state_index, 'j':end_state_index, 
                                       'type':t_type, 
                                       'a':params[0], 'b':params[1], 't':params[2]}]
        elif t_type in ['time_dependent_weibull', 'time_dependent_gompertz']:
            time_dependent_indicies += [{'start_state':t[1], 'end_state':t[2], 
                                         'i':start_state_index, 'j':end_state_index, 
                                         'type':t_type, 
                                         'const':params[0], 'ancillary':params[1]}]
        elif t_type in ['probabilistic_time_dependent_weibull', 'probabilistic_time_dependent_gompertz']:
            # Is treated the same as a time-dependent weibull or gompertz, hence why the type in the dictionary is adjusted to chop off the "probabilistic" part
            # sampled const and ancillary are default the listed "mean" const and ancillary
            probabilistic_time_dependent_indicies += [{'start_state':t[1],'end_state':t[2],
                                                       'i':start_state_index, 'j':end_state_index,
                                                       'type':t_type[14:], 
                                                       'const':params[0], 'ancillary':params[1],
                                                       'se_const':params[2], 'se_ancillary':params[3],
                                                       'sampled_const':params[0], 'sampled_ancillary':params[1]}]
        elif t_type == 'residual':
            residual_indicies += [{'start_state':t[1], 'end_state':t[2], 
                                   'i':start_state_index, 'j':end_state_index, 
                                   'type':t_type, 'params':params}]
        else:
            raise ValueError('Invalid transition type provided:',str(t_type),'Please double check excel file specification')

    # Set up and log transitions that utilize the Dirichlet for transition probabilities.
    # These by nature are similar to beta transitions as they are resampled per-iteration
    # Logs indicies of the transition matrix and the parameters of the distribution per Dirichlet "group"
    dirichlet_indicies = {}
    dirichlet_transitions = transitions_df.loc[transitions_df.type == 'dirichlet']
    if len(dirichlet_transitions) > 0:
        # For every start state that requires a dirichlet type
        for start_state in dirichlet_transitions.start_state.unique():
            outbound_transitions = transitions_df.loc[transitions_df.start_state == start_state]
            dirichlet_indicies[start_state] = []
            # Iterate through each of the outbound states
            for t in outbound_transitions.itertuples():
                start_state_index = state_mapping[t[1]] # mapped row number of start state
                end_state_index = state_mapping[t[2]] # mapped column number of end state
                t_type = t[3] # type of transition 
                params = t[4:] # parameters
                dirichlet_indicies[start_state] += [{'start_state':t[1], 'end_state':t[2], 
                                                    'i':start_state_index, 'j':end_state_index, 
                                                    'type':t[3], 
                                                    'a':params[0]}]

    #---------------------------------------------------------------------------------------------------
    # Run the simulation (for a single arm)
    # TODO: could wrap this into a numba function for increase in speed
    #---------------------------------------------------------------------------------------------------
    # Create result logging array/dataframe 
    # Shape: [iteration_number x states x timesteps] where timesteps = number of cycles
    results_log = np.zeros((max_iterations, num_states, num_cycles))

    print('beginning iterations...')
    iteration_times = np.zeros((max_iterations, 1))
    for iteration in range(0, max_iterations):
        start_iteration_time = perf_counter() # begin timer for logging run time of model iteration

        # Initialize the starting population/proportion for each state at beginning of each model iteration
        # I.e. starting with a zero column-vector of dimension [number_of_states x 1]
        # The initial "proportions" are assigned in accordance to the corresponding row of the vector
        # This corresponds with the row/name assignments of the transition matrix.
        population = np.zeros((num_states, 1)) 
        population[state_mapping[name_start_state]] = 1 
        
        # Initialize with population at time 0 (ie. first cycle everyone is in tx_1 or whereever presumablly)
        results_log[iteration, :, 0] = population.reshape(num_states)
        
        # Initialize as 1 at every iteration becasue population at 0 is always the same at the beginning of
        # any given iteration?
        cycle = 1

        # Resample transition probailities if needed (ie. from distributions) - Update matrix as appropriate
        for t in resample_indicies:
            transition_matrix[t['i'], t['j']] = set_transition(t['type'], a=t['a'], b=t['b'])
        
        # Resample transition probabilities that need the rate to prob conversion. 
        for t in resample_r2p_indicies:
            transition_matrix[t['i'], t['j']] = set_transition(t['type'], a=t['a'], b=t['b'], t=t['t'], cycle_length=cycle_length)  
        
        # Resample transition probabilities from dirichlet if needed - Updates matrix as appropriate
        for d in dirichlet_indicies:
            # 1) Compile parameterization of dirichlet (could move up?)
            parameters = []
            for entry in dirichlet_indicies[d]:
                parameters += [entry['a']]
            # 2) Sample from dirichlet
            dirichlet_output = np.random.dirichlet(parameters) 
            # 3) Update transition probability at appropriate indicied in the transiton matrix
            for t, tp in zip(dirichlet_indicies[d], dirichlet_output):
                transition_matrix[t['i'], t['j']] = tp

        # Resample time-dependent transition probabilities if needed (ie. from distributions)
        # First point of modification if we wish to handle per-interation probabilistic estimates
        for t in probabilistic_time_dependent_indicies:
            t['sampled_const'] = np.random.normal(t['const'], t['se_const'])
            t['sampled_ancillary'] = np.random.normal(t['ancillary'], t['se_ancillary']) 

        # TODO: Contemplate the functionality of this. It currently has pretty specific rules of use
        # Hacky way of dealing with multiple probabilistic Time-dependent states that share the same 
        # sampled parameters... Uses the first instance of the sampled parameters 
        # Based on the assumption that states share the same start-state and end-state name once numbers
        # are stripped.
        # (Note. this could be folded up to the above loop. Keeping separate for now...)
        first_instance_sampled = {}
        for t in probabilistic_time_dependent_indicies:
            start_state_stripped = ''.join([i for i in t['start_state'] if not i.isdigit()])
            end_state_stripped = ''.join([i for i in t['end_state'] if not i.isdigit()])
            pair_id = start_state_stripped + end_state_stripped
            if pair_id in first_instance_sampled:
                t['sampled_const'] = first_instance_sampled[pair_id]['sampled_const']
                t['sampled_ancillary'] = first_instance_sampled[pair_id]['sampled_ancillary']
            else:
                first_instance_sampled[pair_id] = {'sampled_const':t['sampled_const'], 
                                                   'sampled_ancillary':t['sampled_ancillary']}

        # For every timestep until max is reached
        while cycle < num_cycles:
            # Adjust time-dependent transition probabilities based on timestep if needed
            for t in time_dependent_indicies:
                transition_matrix[t['i'], t['j']] = set_transition(t['type'], const=t['const'], ancillary=t['ancillary'], cycle=cycle, cycle_length=cycle_length)
            
            # Adjust proabilistic time-dependent transition probabilities based on timestep if needed (Triggers the same set_transition "pathway" 
            # as a non-probabilistic time-dependent, but passes the sampled parameters instead)
            for t in probabilistic_time_dependent_indicies:
                transition_matrix[t['i'], t['j']] = set_transition(t['type'], const=t['sampled_const'], ancillary=t['sampled_ancillary'], cycle=cycle, cycle_length=cycle_length)

            # Calculate residual transition probailities if needed - Update matrix as appropriate
            for t in residual_indicies:
                transition_matrix[t['i'], t['j']] = 0 # AFTER RESAMPLING RESIDUAL IS RECALCULATED FROM ZERO (but only zero specific transition i to j)
                transition_matrix[t['i'], t['j']] = calculate_residual(transition_matrix, t['i'])
       
            # Normalize to a valid matrix and perform sum check
            transition_matrix = normalize_transitions(transition_matrix) 
            check_row_sums(transition_matrix)
            
            # Calculate the movement for the "population" in each state to another
            # population is already a column vector, dont need to transpose.
            population = transition_matrix.T @ population
            
            # Check: does population sum to 1? (assuming a round to 5 significant digits hold)
            check_model_population(population)
            
            # Update results_log after ever timestep
            results_log[iteration, :, cycle] = population.reshape(num_states) 
            
            cycle += 1 # move to next cycle

        iteration_times[iteration] = [perf_counter() - start_iteration_time] # log time required to run interation

    print('model done...')
    print('total time:',round(iteration_times.sum(), 2),'seconds || mean time per iteration:', round(iteration_times.mean(), 2),'seconds') 

    return results_log


def calculate_costs(population_array, model_specification: ModelSpec, mode='uncorrected'):
    '''
    Function to apply costs to an existing population array. Requires a number of details regarding the
    original model to properly apply.

    Inputs: population_array = [num_iterations x num_states x num_cycles] shaped numpy array representing
                               model output with respect to population movement
            model_specification = instance of the ModelSpec object with a loaded specification
            mode = the method in which to calculate the utility output. default is uncorrected. 
                   "trapezoid" can be specified to use the trapezoid method for half-cycle corrections
    Output: results_cost = [num_iterations x num_states x num_cycles] representing the cost at each state, 
                           in each cycle, relative to the popoulation in a given state.
            *Note: if mode='trapezoid' the output result_cost will have a shape of [num_iterations x num_states x num_cycles-1]
                   this is because of the way it is calculated.
    '''
    if isinstance(model_specification, ModelSpec) == False:
        raise TypeError('Error: Expected input model_specification to be of type markov_modeling.ModelSpec')
    if len(model_specification.state_mapping) != population_array.shape[1]:
        raise ValueError('Number of states represented in population array (', population_array.shape[1], ') does not match the number of states in specification (', len(model_specification.state_mapping), ')')

    discount_rate = model_specification.discount_rate
    cycle_length = model_specification.cycle_length
    state_mapping = model_specification.state_mapping
    costs = model_specification.cost_specs

    iteration_costs = np.zeros((population_array.shape[0], population_array.shape[1]))
    copy_costs = {}
    for t in costs.itertuples():
        state_index = state_mapping[t[1]]
        c_type = t[2]
        if c_type == 'beta':
            for i in range(0,iteration_costs.shape[0]):
                iteration_costs[i, state_index] = get_beta(t[3], t[4])
        elif c_type == 'gamma':
            for i in range(0,iteration_costs.shape[0]):
                iteration_costs[i, state_index] = get_gamma(t[3], t[4])
        elif c_type == 'static':
            iteration_costs[:, state_index] = t[3]
        elif c_type == 'copy':
            copy_costs[state_index] = [state_mapping[t[3]]] # cost for index (noted as copy) points to index of target
            iteration_costs[:, state_index] = np.nan # initialize as nan
        else:
            raise ValueError('Error: Bad cost type specification', c_type)

    # Fill copied costs
    for c in copy_costs:
        iteration_costs[:, c] = iteration_costs[:, copy_costs[c]].reshape(iteration_costs.shape[0])
    if np.isnan(iteration_costs).any(): #TODO: Could make this output something more detailed. I.e. index of error cost
        raise ValueError('Error: NaN cost specified. Likely a copy type error. Please check costs sheet')

    # Apply sampled costs to population array
    if mode == 'trapezoid':
        # Calculate state membership via the trapezoid method
        population_array_trap = (population_array + np.roll(population_array, -1)) / 2
        population_array_trap = population_array_trap[:, :, :-1]

        results_costs = np.zeros(population_array_trap.shape)
        for state in state_mapping:
            idx = state_mapping[state]  
            results_costs[:,idx,:] = population_array_trap[:,idx,:] * iteration_costs[:,idx][:,np.newaxis]

        # Apply discount rate with an extra half-cycle tacked on (is the half cycle needed?)
        # cost * (1 / ((1+discount_rate)**year))
        for i in range(0,results_costs.shape[2]):
            year = math.floor(((i*cycle_length)+(cycle_length/2))/365)
            results_costs[:,:,i] = results_costs[:,:,i] * (1 / ((1+discount_rate)**year))

        return results_costs
    else:
        results_costs = np.zeros(population_array.shape)

        for state in state_mapping:
            idx = state_mapping[state]  
            results_costs[:,idx,:] = population_array[:,idx,:] * iteration_costs[:,idx][:,np.newaxis]
                
        # Apply discount rate
        # cost * (1 / ((1+discount_rate)**year))
        for i in range(0,results_costs.shape[2]):
            year = math.floor((i*cycle_length)/365)
            results_costs[:,:,i] = results_costs[:,:,i] * (1 / ((1+discount_rate)**year))

        return results_costs


def calculate_utilities(population_array, model_specification: ModelSpec, mode='uncorrected'):
    '''
    Function to apply utilities to an existing population array. Requires a number of details regarding the
    original model to properly apply.

    Inputs: population_array = [num_iterations x num_states x num_cycles] shaped numpy array representing
                               model output with respect to population movement
            model_specification = instance of the ModelSpec object with a loaded specification
            mode = the method in which to calculate the utility output. default is uncorrected.
                   "trapezoid" can be specified to use the trapezoid method for half-cycle corrections
    Output: results_utilities = [num_iterations x num_states x num_cycles] representing the utility at each state, 
                                in each cycle, relative to the popoulation in a given state.
            *Note: if mode='trapezoid' the output result_cost will have a shape of [num_iterations x num_states x num_cycles-1]
                   this is because of the way it is calculated.
    '''
    if isinstance(model_specification, ModelSpec) == False:
        raise TypeError('Error: Expected input model_specification to be of type markov_modeling.ModelSpec')
    if len(model_specification.state_mapping) != population_array.shape[1]:
        raise ValueError('Number of states represented in population array (', population_array.shape[1], ') does not match the number of states in specification (', len(model_specification.state_mapping), ')')

    discount_rate = model_specification.discount_rate
    cycle_length = model_specification.cycle_length
    state_mapping = model_specification.state_mapping
    utilities = model_specification.util_specs

    iteration_utils = np.zeros((population_array.shape[0], population_array.shape[1]))
    copy_utils = {}
    for t in utilities.itertuples():
        state_index = state_mapping[t[1]]
        u_type = t[2]
        if u_type == 'beta':
            for i in range(0,iteration_utils.shape[0]):
                iteration_utils[i, state_index] = get_beta(t[3], t[4])
        elif u_type == 'beta_disutility':
            for i in range(0,iteration_utils.shape[0]):
                iteration_utils[i, state_index] = 1-get_beta(t[3], t[4])
        elif u_type == 'gamma':
            for i in range(0,iteration_utils.shape[0]):
                iteration_utils[i, state_index] = get_gamma(t[3], t[4])
        elif u_type == 'gamma_disutility':
            for i in range(0,iteration_utils.shape[0]):
                iteration_utils[i, state_index] = 1-get_gamma(t[3], t[4])
        elif u_type == 'static':
            iteration_utils[:, state_index] = t[3]
        elif u_type == 'copy':
            copy_utils[state_index] = [state_mapping[t[3]]] # cost for index (noted as copy) points to index of target
            iteration_utils[:, state_index] = np.nan # initialize as nan
        else:
            raise ValueError('Error: Bad utility type specification', u_type)
    
    # Fill copied utilities
    for u in copy_utils:
        iteration_utils[:, u] = iteration_utils[:, copy_utils[u]].reshape(iteration_utils.shape[0])
    if np.isnan(iteration_utils).any(): #TODO: Could make this output something more detailed. I.e. index of error utility
        raise ValueError('Error: NaN utility specified. Likely a copy type error. Please check utility sheet')

    # Apply sampled utilities to population array.
    if mode == 'trapezoid':
        # Calculate state membership via the trapezoid method
        population_array_trap = (population_array + np.roll(population_array, -1)) / 2
        population_array_trap = population_array_trap[:, :, :-1]

        results_utilities = np.zeros(population_array_trap.shape)
        for state in state_mapping:
            idx = state_mapping[state]
            results_utilities[:,idx,:] = population_array_trap[:,idx,:] * iteration_utils[:,idx][:,np.newaxis]

        # Apply discount rate with an extra half-cycle tacked on (is the half cycle needed?)
        # cost * (1 / ((1+discount_rate)**year))
        for i in range(0,results_utilities.shape[2]):
            year = math.floor(((i*cycle_length)+(cycle_length/2))/365)
            results_utilities[:,:,i] = results_utilities[:,:,i] * (1 / ((1+discount_rate)**year))

        # Adjust utilities to per-year
        results_utilities = results_utilities * (cycle_length/365)
        
        return results_utilities
    else:
        results_utilities = np.zeros(population_array.shape)

        for state in state_mapping:
            idx = state_mapping[state]
            results_utilities[:,idx,:] = population_array[:,idx,:] * iteration_utils[:,idx][:,np.newaxis]

        # Apply discount rate
        # cost * (1 / ((1+discount_rate)**year))
        for i in range(0,results_utilities.shape[2]):
            year = math.floor((i*cycle_length)/365)
            results_utilities[:,:,i] = results_utilities[:,:,i] * (1 / ((1+discount_rate)**year))

        # Adjust utilities to per-year
        results_utilities = results_utilities * (cycle_length/365)
        
        return results_utilities


def calculate_icer(base: ModelData, treat: ModelData, calculate_ci=False):
    '''
    Convenience function to calculate an ICER and optional CI using bias-corrected and accelerated
    method. 

    Inputs: base = ModelData object containing the data outputs (population, costs, utils) from the base case arm
            treat = ModelData object containing the data outputs from the treatment (comparator) arm

    Outputs: ICER = the calculated ICER
                CI = the calculated CI of the ICER using the bias-corrected and accelerated bootstrap implmentation
                    in the ARCH package. (This is only returned when calculate_ci is set to True)
    '''
    if isinstance(base, ModelData) == False or isinstance(treat, ModelData) == False:
        raise TypeError('Error: Expected inputs base and treat to be of type markov_modeling.ModelData')

    # Calculate delta utility and delta cost between the iterations of the treatment and base case arms
    delta_mean_utility = treat.iteration_util_data.mean() - base.iteration_util_data.mean() #Average util of treat - average util of basecase
    delta_mean_cost = treat.iteration_cost_data.mean() - base.iteration_cost_data.mean() #Average cost of treat - average cost of basecase

    ICER = delta_mean_cost/delta_mean_utility

    if calculate_ci == True:
        results_data = pd.DataFrame({'cost_treat': treat.iteration_cost_data,
                                    'cost_base': base.iteration_cost_data,
                                    'utility_treat': treat.iteration_util_data,
                                    'utility_base': base.iteration_util_data})

        def func_icer(x):
            '''
            Function to return the ICER of the average of the boot strap sample.
            Input: x = rows of results_dataframe that were indexed from the bootstrap sample
                    has 4 columns for cost_treat, cost_base, utility_treat, utility_base
            Output: ICER calculated from the iterations present in x
            
                    mean(cost_treat) - mean(cost_base)
            ICER = ----------------------------------------
                mean(utility_treat) - mean(utility_base)
            '''
            return (x['cost_treat'].mean() - x['cost_base'].mean()) / (x['utility_treat'].mean() - x['utility_base'].mean())

        bs = IIDBootstrap(results_data) #use a "dummy" of array indicies to sample from. Needed to correctly calculate ICER of the average
        ci = bs.conf_int(func_icer, 1000, method='bca') #bias-corrected and accelerated method

        return ICER, ci

    return ICER
