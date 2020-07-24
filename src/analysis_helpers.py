#---------------------------------------------------------------------------------------------------
# Helper functions for cost-effectivness analysis
# Brandon Chan - July 2020
#---------------------------------------------------------------------------------------------------
# Import packages
#---------------------------------------------------------------------------------------------------
import numpy as np 
import math 
import matplotlib.pyplot as plt
import seaborn as sns
from markov_modeling import *
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

    pop_base, cost_base, util_base = run_model('../model_specifications/test_parameters_base.xlsx', save=save, model_name='base_'+analysis_name)
    pop_treat, cost_treat, util_treat = run_model('../model_specifications/../model_specifications/test_parameters_treat.xlsx', save=save, model_name='treatment_'+analysis_name)

    # calculate sum at each time point (assign to bottom row)
    # result array shape: [num_iterations, num_timesteps]
    # each row is the sum of the costs/utilitys at time point t of an iteration
    cost_sum_per_cycle_base = np.sum(cost_base, axis=1)
    utility_sum_per_cycle_base = np.sum(util_base, axis=1)

    cost_sum_per_cycle_treat = np.sum(cost_treat, axis=1)
    utility_sum_per_cycle_treat = np.sum(util_treat, axis=1)

    # Adjusting utility to QALY 
    # = cycle_length /  days_in_a_year (where cycle_length is in days)
    utility_sum_per_cycle_base = utility_sum_per_cycle_base * (cycle_length/365)
    utility_sum_per_cycle_treat = utility_sum_per_cycle_treat * (cycle_length/365)

    # Calculate the total cost/utility of an iteration
    cost_sum_per_iteration_base = np.sum(cost_sum_per_cycle_base, axis=1)
    utility_sum_per_iteration_base = np.sum(utility_sum_per_cycle_base, axis=1)
    cost_sum_per_iteration_treat = np.sum(cost_sum_per_cycle_treat, axis=1)
    utility_sum_per_iteration_treat = np.sum(utility_sum_per_cycle_treat, axis=1)

    # Difference in costs/utilities when comparing each arm per iteration
    c = cost_sum_per_iteration_treat - cost_sum_per_iteration_base
    u = utility_sum_per_iteration_treat - utility_sum_per_iteration_base

    delta_mean_utility = utility_sum_per_iteration_treat.mean() - utility_sum_per_iteration_base.mean() #Average util of treat - average util of basecase
    delta_mean_cost = cost_sum_per_iteration_treat.mean() - cost_sum_per_iteration_base.mean() #Average cost of treat - average cost of basecase

    def func(x):
        '''
        Function to return the ICER of the average of the boot strap sample.
        Uses globally defined cost and utility arrays (c and u):
            c and u are [1 x number_of_iterations] numpy arrays that represent 
            the sum of a iterations costs and utility respectivley
        x = indicies of interations to use from bootstrap function
        '''
        return c[x].mean()/u[x].mean()

    bs = IIDBootstrap(np.array(list(range(0,c.shape[0])))) #use a "dummy" of array indicies to sample from. Needed to correctly calculate ICER of the average
    ci = bs.conf_int(func, 1000, method='bca') #bias-corrected and accelerated method

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

    #  Checks: population_array.shape[1] == number of states in mapping
    #          all state mappings exist in costs
    results_costs = np.zeros(population_array.shape)

    for state in state_mapping:
        idx = state_mapping[state]  
        results_costs[:,idx,:] = results_costs[:,idx,:] * costs.loc[costs.state==state].cost.values[0]

    # Apply discount rate
    # cost * (1 / ((1+discount_rate)**year))
    for i in range(0,results_costs.shape[2]):
        year = math.floor((i*cycle_length)/365)
        results_costs[:,:,i] = results_costs[:,:,i] * (1 / ((1+discount_rate)**year))

    return results_costs

def apply_new_utilities(population_array, discount_rate, cycle_length, state_mapping, utilities):
    
    # Checks: population_array.shape[1] == number of states in mapping 
    #         all state mappings exist in utilities
    
    results_utilities = np.zeros(population_array.shape)
    
    for state in state_mapping:
        idx = state_mapping[state]
        results_utilities[:,idx,:] = results_utilities[:,idx,:] * utilities.loc[utilities.state==state].cost.values[0]

    # Apply discount rate
    # cost * (1 / ((1+discount_rate)**year))
    for i in range(0,results_utilities.shape[2]):
        year = math.floor((i*cycle_length)/365)
        results_utilities[:,:,i] = results_utilities[:,:,i] * (1 / ((1+discount_rate)**year))
    
    return results_utilities