#---------------------------------------------------------------------------------------------------
# Helper functions for cost-effectivness analysis
# Brandon Chan - July 2020
#---------------------------------------------------------------------------------------------------
# Import packages
#---------------------------------------------------------------------------------------------------
import numpy as np 
import math 
from arch.bootstrap import IIDBootstrap # Bootstrap analysis

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