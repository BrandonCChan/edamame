#---------------------------------------------------------------------------------------------------
# Helper functions for cost-effectivness analysis
# Brandon Chan - July 2020
#---------------------------------------------------------------------------------------------------
# Import packages
#---------------------------------------------------------------------------------------------------
import numpy as np 
import pandas as pd
import math 
import matplotlib.pyplot as plt
import seaborn as sns
from markov_modeling import run_model
from model_helpers import get_beta, get_gamma
from arch.bootstrap import IIDBootstrap # Bootstrap analysis

def cea_analysis(base_filepath, treatment_filepath, save=False, analysis_name='', ce_threshold=20000):
    '''
    Function to run cost effectivness analysis. 
    Inputs: base_filepath = filepath to specification of the base case / standard of care arm
            treatment_filepath = filepath to the specification of the treatment arm
    Optional Inputs: save = boolean (true/false) flag for saving model outputs. Default is False
                     analysis_name = string to specify a descriptor to add to the filenames of saved outputs
                     ce_threshold = integer to denote the willingness to pay threshold.
    '''
    #---------------------------------------------------------------------------------------------------
    # Read in excel workbook and assign each sheet to a dataframe
    #---------------------------------------------------------------------------------------------------
    specification_df = pd.read_excel(base_filepath, 'specification', header=None, index_col=0)
    cycle_length = specification_df.loc['cycle_length'].values[0]

    _, cost_base, util_base = run_model('../model_specifications/test_parameters_base.xlsx', save=save, model_name='base_'+analysis_name)
    _, cost_treat, util_treat = run_model('../model_specifications/../model_specifications/test_parameters_treat.xlsx', save=save, model_name='treatment_'+analysis_name)

    # calculate sum at each time point (assign to bottom row)
    # result array shape: [num_iterations, num_timesteps]
    # each row is the sum of the costs/utilitys at time point t of an iteration
    cost_sum_per_cycle_base = np.sum(cost_base, axis=1)
    utility_sum_per_cycle_base = np.sum(util_base, axis=1)

    cost_sum_per_cycle_treat = np.sum(cost_treat, axis=1)
    utility_sum_per_cycle_treat = np.sum(util_treat, axis=1)

    # Calculate the total cost/utility of an iteration
    cost_sum_per_iteration_base = np.sum(cost_sum_per_cycle_base, axis=1)
    utility_sum_per_iteration_base = np.sum(utility_sum_per_cycle_base, axis=1)
    cost_sum_per_iteration_treat = np.sum(cost_sum_per_cycle_treat, axis=1)
    utility_sum_per_iteration_treat = np.sum(utility_sum_per_cycle_treat, axis=1)

    # Difference in costs/utilities when comparing each arm per iteration
    c = cost_sum_per_iteration_treat - cost_sum_per_iteration_base
    u = utility_sum_per_iteration_treat - utility_sum_per_iteration_base

    # Dataframe for ICER CI / Bootstrap Calculations
    results_dataframe = pd.DataFrame({'cost_treat':cost_sum_per_iteration_treat,
                                      'cost_base':cost_sum_per_iteration_base,
                                      'utility_treat':utility_sum_per_iteration_treat,
                                      'utility_base':utility_sum_per_iteration_base})

    delta_mean_utility = utility_sum_per_iteration_treat.mean() - utility_sum_per_iteration_base.mean() #Average util of treat - average util of basecase
    delta_mean_cost = cost_sum_per_iteration_treat.mean() - cost_sum_per_iteration_base.mean() #Average cost of treat - average cost of basecase

    def function_icer(x):
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

    bs = IIDBootstrap(results_dataframe) #use a "dummy" of array indicies to sample from. Needed to correctly calculate ICER of the average
    ci = bs.conf_int(function_icer, 1000, method='bca') #bias-corrected and accelerated method

    print('ICER with 95% CI:')
    print(round(delta_mean_cost/delta_mean_utility,2),'(',round(ci[0][0],2),'to',round(ci[1][0],2),')')

    # Make the CE plane plot
    plt.axhline(0, color='black')
    plt.axvline(0, color='black')
    plt.ylabel(r'$\Delta$C', rotation=0, fontsize=14, labelpad=14)
    plt.xlabel(r'$\Delta$E', fontsize=14)
    sns.scatterplot(x=u, y=c)

    # Plots some threshold line...
    x = np.linspace(min(u), max(u), 100)
    y = ce_threshold*x
    plt.plot(x, y, '-r')

    plt.ylim(min(c)-50, max(c)+50) # Add padding to axis limits

    plt.show()


def apply_new_costs(population_array, discount_rate, cycle_length, state_mapping, costs):
    '''
    Function to apply new costs to an existing population array. Requires a number of details regarding the
    original model to properly apply.
    Inputs: population_array = [num_iterations x num_states x num_cycles] shaped numpy array representing
                               model output with respect to population movement
            discount_rate = the discount rate you want to apply to the costs (represented as a number 0 < n < 1)
            cycle_length = cycle length in days of the original model run
            state_mapping = dictionary mapping the name of a state to an index number
            costs = dataframe of cost data to be applied to the population array
    Output: results_cost = [num_iterations x num_states x num_cycles] representing the cost at each state, 
                           in each cycle, relative to the popoulation in a given state.
    '''
    #  Checks: population_array.shape[1] == number of states in mapping
    #          all state mappings exist in costs
    results_costs = np.zeros(population_array.shape)

    iteration_costs = np.zeros((population_array.shape[0], population_array.shape[1]))
    for t in costs.itertuples():
        state_index = state_mapping[t[1]]
        c_type = t[2]
        if c_type == 'beta':
            for i in range(0,iteration_costs.shape[0]):
                iteration_costs[:,state_index] = get_beta(t[3], t[4])
        elif c_type == 'gamma':
            for i in range(0,iteration_costs.shape[0]):
                iteration_costs[:,state_index] = get_gamma(t[3], t[4])
        elif c_type == 'static':
            iteration_costs[:,state_index] = t[3]
        else:
            raise ValueError('Error: Bad cost type specification', c_type)

    for state in state_mapping:
        idx = state_mapping[state]  
        results_costs[:,idx,:] = population_array[:,idx,:] * iteration_costs[:,idx][:,np.newaxis]
            
    # Apply discount rate
    # cost * (1 / ((1+discount_rate)**year))
    for i in range(0,results_costs.shape[2]):
        year = math.floor((i*cycle_length)/365)
        results_costs[:,:,i] = results_costs[:,:,i] * (1 / ((1+discount_rate)**year))

    return results_costs

def apply_new_utilities(population_array, discount_rate, cycle_length, state_mapping, utilities):
    '''
    Function to apply new utilities to an existing population array. Requires a number of details regarding the
    original model to properly apply.
    Inputs: population_array = [num_iterations x num_states x num_cycles] shaped numpy array representing
                               model output with respect to population movement
            discount_rate = the discount rate you want to apply to the costs (represented as a number 0 < n < 1)
            cycle_length = cycle length in days of the original model run
            state_mapping = dictionary mapping the name of a state to an index number
            utilities = dataframe of cost data to be applied to the population array
    Output: results_utilities = [num_iterations x num_states x num_cycles] representing the utility at each state, 
                                in each cycle, relative to the popoulation in a given state.
    '''
    # Checks: population_array.shape[1] == number of states in mapping 
    #         all state mappings exist in utilities
    
    results_utilities = np.zeros(population_array.shape)

    iteration_utils = np.zeros((population_array.shape[0], population_array.shape[1]))
    for t in utilities.itertuples():
        state_index = state_mapping[t[1]]
        u_type = t[2]
        if u_type == 'beta':
            for i in range(0,iteration_utils.shape[0]):
                iteration_utils[i,state_index] = get_beta(t[3], t[4])
        elif u_type == 'gamma':
            for i in range(0,iteration_utils.shape[0]):
                iteration_utils[i,state_index] = get_gamma(t[3], t[4])
        elif u_type == 'static':
            iteration_utils[:,state_index] = t[3]
        else:
            raise ValueError('Error: Bad utility type specification', u_type)
    
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