#---------------------------------------------------------------------------------------------------
# Code for running a monte carlo markov using a model structure and parameters specified in an excel
# document. (See README and test_parameters.xlsx for addional explaination and example)
#
# Brandon Chan | July 2020
#---------------------------------------------------------------------------------------------------
# Import packages/libraries
#---------------------------------------------------------------------------------------------------
import pandas as pd # Dataframe structure and manipulation
import numpy as np # Scientific computing functionality (array and matrix operations and some stats things)
import math # For additional math functions
import os
from time import perf_counter, strftime
from openpyxl import load_workbook
from model_helpers import * 

def run_model(filepath, save=False, model_name='model'):
    '''
    Function that reads a specified excel sheet and runs the model specified by it
    Input: filepath = path of excel doc that defines model
    Optional inputs: save = boolean (True/False) flag to automatically save the model outputs in the
                            model_outputs/ directory. Default is Flase
                     model_name = Used when save is equal to True. Specifies a string to be used to
                                    as a in-filename descriptor for the saved model outputs. Default is "model" 
    Output: returns 3 npy arrays of the population, cost, and utility over each model iteration, 
            state, and cycle.
    '''
    #---------------------------------------------------------------------------------------------------
    # File IO and parameter initialization
    #---------------------------------------------------------------------------------------------------
    FILE_PATH = os.path.abspath(os.path.dirname(__file__))
    input_file = pd.ExcelFile(FILE_PATH + "\\model_specifications\\" + filepath) # read in excel file
    check_excel_file(input_file) # run checks on file formatting and specification of model

    # Read in each speadsheet in the file
    transitions_df = pd.read_excel(input_file, 'transitions')
    costs_df = pd.read_excel(input_file, 'costs')
    utilities_df = pd.read_excel(input_file, 'utilities')
    specification_df = pd.read_excel(input_file, 'specification', header=None, index_col=0)

    # Specification of variables regarding states in model
    unique_states = transitions_df['start_state'].unique().tolist() #list(set(start_state_names+end_state_names))
    num_states = len(unique_states)

    # Initialize model parameters based on spreadsheet defined values
    max_iterations = int(specification_df.loc['max_iterations'].values[0])
    cycle_length = specification_df.loc['cycle_length'].values[0]
    time_horizon_days = specification_df.loc['time_horizon'].values[0] * 365 # If not 360 / 365 adjust accordingly
    num_cycles = int(time_horizon_days / cycle_length)
    name_start_state = specification_df.loc['name_start_state'].values[0]
    discount_rate = specification_df.loc['discount_rate'].values[0]

    print('file and parameters loaded...')

    #---------------------------------------------------------------------------------------------------
    # Generate matrix representation of model
    #
    # Dimensions of matrix = [num_states x num_states]
    # Rows map to starting state, columns map to target state
    #---------------------------------------------------------------------------------------------------

    # Use dict to map from names to numeric index in array 
    state_mapping = {i : unique_states.index(i) for i in unique_states}
    #specify empty transition matrix
    transition_matrix = np.zeros((num_states,num_states)) 

    # Keep track of specific indicies in the transition matrix that need to be updated "in-simulation"
    # ie. time-dependant transitions and transitions that get resampled every iteration
    # Intended to be stored as a list of dictionaries
    resample_indicies = []
    time_dependant_indicies = []

    # Iterate through specified transtions and initialize values of matrix 
    # Log indicies of transitions that need to be updated to quickly index the correct position in the transition matrix
    for t in transitions_df.itertuples():
        start_state_index = state_mapping[t[1]] # mapped row number of start state
        end_state_index = state_mapping[t[2]] # mapped column number of end state
        t_type = t[3] # type of transition 
        params = t[4:] # parameters
        
        if t_type in ['beta', 'gamma']:
            resample_indicies += [{'start_state':t[1],'end_state':t[2],'i':start_state_index, 'j':end_state_index, 'type':t_type, 'params':params}]
        elif t_type == 'time-dependant':
            time_dependant_indicies += [{'start_state':t[1],'end_state':t[2],'i':start_state_index, 'j':end_state_index, 'type':t_type, 'params':params}]
        elif t_type == 'residual':
            residual_indicies += [{'start_state':t[1],'end_state':t[2],'i':start_state_index, 'j':end_state_index, 'type':t_type, 'params':params}]

        # Initializes the transition matrix with values
        # Optional? - no, we probally want to assign any static transitions here as they would otherwise not be set
        transition_matrix[start_state_index, end_state_index] = set_transition(t_type,params)

    # rescale row-wise to ensure row sums (ie. all transitions out of a state sum to 1)
    transition_matrix = normalize_transitions(transition_matrix) 

    #---------------------------------------------------------------------------------------------------
    # Run the simulation (for a single arm)
    #---------------------------------------------------------------------------------------------------
    # Create result logging array/dataframe 
    # Shape: [iteration_number x states x timesteps] where timesteps = number of cycles
    results_log = np.zeros((max_iterations, num_states, num_cycles))

    print('beginning iterations...')
    iteration_times = np.zeros((max_iterations,1))
    for iteration in range(0,max_iterations):
        start_iteration_time = perf_counter()
        # Initialize the starting population/proportion for each state at beginning of each model iteration
        # I.e. starting with a zero column-vector of dimension [number_of_states x 1]
        # The initial "proportions" are assigned in accordance to the corresponding row of the vector
        # This corresponds with the row/name assignments of the transition matrix.
        population = np.zeros((num_states,1)) 
        population[state_mapping[name_start_state]] = 1 
        '''
        ^ FLAGGED FOR CLARIFICATION 
        - Likely needs some work or initial specification from the excel file built-in
        - TODO: add specification for model start state (ie. which state(s) do we assign the population to at t0?)
        '''
        
        # Initialize with population at time 0 (ie. first cycle everyone is in tx_1 or whereever presumablly)
        results_log[iteration, :, 0] = population.reshape(num_states)
        
        # Initialize as 1 at every iteration becasue population at 0 is always the same at the beginning of
        # any given iteration?
        cycle = 1
        time = 0

        # Resample transition probailities if needed (ie. from distributions) - Update matrix as appropriate
        for t in resample_indicies:
            transition_matrix[t['i'],t['j']] = set_transition(t['type'],t['params'])  
        transition_matrix = normalize_transitions(transition_matrix) # normalize sampling
        check_row_sums(transition_matrix) # Check if row sums are == 1
        
        #----------------------------------------------------------------------------------
        # For every timestep until max is reached
        while cycle < num_cycles:
            # Adjust time-dependant transition probabilities based on timestep if needed
            for t in time_dependant_indicies:
                transition_matrix[t['i'],t['j']] = set_transition(t['type'],t['params']+[time, cycle_length])
            transition_matrix = normalize_transitions(transition_matrix) # Update matrix (certian i,j based on variable time)? - TODO this "targeted" functionality
            check_row_sums(transition_matrix) # Check if row sums are == 1
            
            # Initialize a temperary zero vector to log the updated population proportions
            new_population = np.zeros(population.shape)
            
            #------------------------------------------------------------------------------
            # Calculate the movement for the "population" in each state to another
            # ie. sum them up per state
            # population as a col vector basically 
            for i in range(0,population.shape[0]):
                movements = transition_matrix[i,:] # isolate possible transitions for the specific state
                
                #--------------------------------------------------------------------------
                # apply/redistribute the population based on the current value of the state and the defined 
                # transition probabilities
                for j in range(0,movements.shape[0]):
                    new_population[j] += population[i]*movements[j]
            
            population = new_population # Assign updated population-proportion numbers
            
            # Check: does population sum to 1? (assuming a round to 5 significant digits hold)
            check_model_population(population)
            
            # Update results_log after ever timestep
            results_log[iteration, :, cycle] = population.reshape(num_states) 
            
            time += cycle_length # increment time based on cycle length
            cycle += 1 # next cycle

        iteration_times[iteration] = [perf_counter() - start_iteration_time]
    print('model done...')
    print('total time:',round(iteration_times.sum(),2),'seconds || mean time per iteration:', round(iteration_times.mean(),2),'seconds')      
    #---------------------------------------------------------------------------------------------------
    # Costing and utility component
    #---------------------------------------------------------------------------------------------------
    # TODO: implement condensable and substate calcualtions
    condensable_state_mappings = {}
    multiple_toxicity_states = {}

    # Theroetically "faster" to move this into the proportion loops?
    results_log_costs = np.zeros(results_log.shape)
    results_log_utilities = np.zeros(results_log.shape)

    print('calculating costs and utilities...')
    for state in state_mapping:
        idx = state_mapping[state]
        
        # TODO: add loose string matching or precompute associated indicied for "repeat states" (ie. rest or remission)
        # When mapping to the utility/cost tables, it will match based on the "condensed" state name
        if state in condensable_state_mappings:
            state = condensable_state_mappings['state']
            
        # TODO: add "division" of treatment states for cost or utility? How to effectivley do this...
        if state in multiple_toxicity_states:
            multiple_toxicity_states['state']
            results_log[:,idx,:] # is proportion of pts in that patient at every time iterations
            
            results_log_costs[:,idx,:] = results_log[:,idx,:] * costs_df.loc[costs_df.state==state].cost.values[0] #cost_array[row,iteration]
            results_log_utilities[:,idx,:] = results_log[:,idx,:] * utilities_df.loc[utilities_df.state==state].utility.values[0] #utility_array[row,iteration]
        else:
            # multiply every (proportion) entry with the cost/utility associated with that state at that iteration
            # assign result to utility/cost array
            results_log_costs[:,idx,:] = results_log[:,idx,:] * costs_df.loc[costs_df.state==state].cost.values[0] #cost_array[row,iteration]
            results_log_utilities[:,idx,:] = results_log[:,idx,:] * utilities_df.loc[utilities_df.state==state].utility.values[0] #utility_array[row,iteration]
    print('done costs and utilities...')

    # Apply discount rate
    # cost * (1 / ((1+discount_rate)**year))
    for i in range(0,results_log.shape[2]):
        year = math.floor((i*cycle_length)/365) #+ 1 (if year 0 = 1) # if not 360 days. change.
        results_log_costs[:,:,i] = results_log_costs[:,:,i] * (1 / ((1+discount_rate)**year))
        results_log_utilities[:,:,i] = results_log_utilities[:,:,i] * (1 / ((1+discount_rate)**year))

    #---------------------------------------------------------------------------------------------------
    # Save/return result components in a raw form (ie. as 3D arrays)
    # 1) population
    # 2) costs
    # 3) utility
    #---------------------------------------------------------------------------------------------------
    if save:
        book = load_workbook(FILE_PATH + "\\model_specifications\\" + filepath)
        if 'state_mappings' not in book.sheetnames:
        # TODO:// more powerful matching. ie. update instead of ignore if different
            with pd.ExcelWriter(FILE_PATH + "\\model_specifications\\" + filepath, engine='openpyxl') as writer:
                writer.book = book
                state_mapping_df = pd.DataFrame.from_dict(state_mapping, orient='index')
                state_mapping_df.to_excel(writer, sheet_name='state_mappings')

        name = model_name + '_' + strftime("%d-%m-%Y")
        np.save('./model_outputs/'+name+'_population.npy', results_log)
        np.save('./model_outputs/'+name+'_costs.npy', results_log_costs)
        np.save('./model_outputs/'+name+'_utilities.npy', results_log_utilities)
        print('results saved as: ./model_outputs/'+name+'... .npy')
    print('')

    return results_log, results_log_costs, results_log_utilities