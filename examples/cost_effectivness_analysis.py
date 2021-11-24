#---------------------------------------------------------------------------------------------------
# Example code for how to conduct a cost effectivness analysis using
# the framework. Mirrors code detailed in "Cost-Effectivness Analysis.ipynb"
# Please see the jupyter notebook for additional code on calculating other output metrics
#
# Brandon Chan | November 2021
#---------------------------------------------------------------------------------------------------
# Import Packages
#---------------------------------------------------------------------------------------------------
import numpy as np # Scientific computing - used for n-dimensional arrays, numerical functions
import math
import pandas as pd # Dataframe structure and reading of excel files
import matplotlib.pyplot as plt # Plotting package
import seaborn as sns # Another plotting package
import os # for file path things

import sys
sys.path.insert(0,"../src") # Direct to the src directory depending on where you're developing from
from markov_modeling import *

# Optional: specify random seed
np.random.seed(123)

#---------------------------------------------------------------------------------------------------
# Load base and treatment arm specifications, run model, and calculate cost and utility outputs
#---------------------------------------------------------------------------------------------------
specification_base = ModelSpec('../model_specifications/test_parameters_base.xlsx', 
                               model_name='test_base')
pop_base = run_model(specification_base)
cost_base = calculate_costs(pop_base, specification_base)
util_base = calculate_utilities(pop_base, specification_base)

specification_treat = ModelSpec('../model_specifications/test_parameters_treat.xlsx',
                                model_name='test_treat')
pop_treat = run_model(specification_treat)
cost_treat = calculate_costs(pop_treat, specification_treat)
util_treat = calculate_utilities(pop_treat, specification_treat)

#----------------------------------------------
# Optional: Saving outputs of running model and cost/utility calculations
#----------------------------------------------
'''
name = specification_base.model_name + '_' + strftime("%d-%m-%Y")
np.save('../model_outputs/'+name+'_population.npy', pop_base)
np.save('../model_outputs/'+name+'_costs.npy', cost_base)
np.save('../model_outputs/'+name+'_utilities.npy', util_base)

name = specification_treat.model_name + '_' + strftime("%d-%m-%Y")
np.save('../model_outputs/'+name+'_population.npy', pop_treat)
np.save('../model_outputs/'+name+'_costs.npy', cost_treat)
np.save('../model_outputs/'+name+'_utilities.npy', util_treat)
'''
#----------------------------------------------
# Optional: Load existing outputs
#----------------------------------------------
'''
pop_base = np.load('../model_outputs/test_base_20-07-2020_population.npy')
cost_base = np.load('../model_outputs/test_base_20-07-2020_costs.npy')
util_base = np.load('../model_outputs/test_base_20-07-2020_utilities.npy')

pop_treat = np.load('../model_outputs/test_treat_20-07-2020_population.npy')
cost_treat = np.load('../model_outputs/test_treat_20-07-2020_costs.npy')
util_treat = np.load('../model_outputs/test_treat_20-07-2020_utilities.npy')
'''

#---------------------------------------------------------------------------------------------------
# Use instances of ModelData to store and automatically condense model outputs
#---------------------------------------------------------------------------------------------------
outputs_base = ModelData(pop_base, cost_base, util_base)
outputs_treat = ModelData(pop_treat, cost_treat, util_treat)

#---------------------------------------------------------------------------------------------------
# Print basic model outputs
#---------------------------------------------------------------------------------------------------
significant_digits = 2 #convenience variable to format print statement numbers 
n_sample = 1

#-----------------------------
# Cost stats:
#-----------------------------
mean_cost_treat = np.mean(outputs_treat.iteration_cost_data)
std_cost_treat = np.std(outputs_treat.iteration_cost_data)
up_CI_cost_treat = mean_cost_treat + (1.96*(std_cost_treat/math.sqrt(n_sample)))
low_CI_cost_treat = mean_cost_treat - (1.96*(std_cost_treat/math.sqrt(n_sample)))

mean_cost_base = np.mean(outputs_base.iteration_cost_data)
std_cost_base = np.std(outputs_base.iteration_cost_data) 
up_CI_cost_base = mean_cost_base + (1.96*(std_cost_base/math.sqrt(n_sample)))
low_CI_cost_base = mean_cost_base - (1.96*(std_cost_base/math.sqrt(n_sample)))

print('--------------------------------------------------------------')
print('Treatment arm costs: ',round(mean_cost_treat, 2),'[',round(low_CI_cost_treat, 2),',',round(up_CI_cost_treat,2),']', '| SD:',round(std_cost_treat,2))
print('Base arm costs:      ',round(mean_cost_base, 2),'[',round(low_CI_cost_base, 2),',',round(up_CI_cost_base,2),']', '| SD:',round(std_cost_base,2))
print('--------------------------------------------------------------')

#-----------------------------
# Utility stats:
#-----------------------------
mean_util_treat = np.mean(outputs_treat.iteration_util_data)
std_util_treat = np.std(outputs_treat.iteration_util_data)
up_CI_util_treat = mean_util_treat + (1.96*(std_util_treat/math.sqrt(n_sample)))
low_CI_util_treat = mean_util_treat - (1.96*(std_util_treat/math.sqrt(n_sample)))

mean_util_base = np.mean(outputs_base.iteration_util_data)
std_util_base = np.std(outputs_base.iteration_util_data) 
up_CI_util_base= mean_util_base + (1.96*(std_util_base/math.sqrt(n_sample)))
low_CI_util_base = mean_util_base - (1.96*(std_util_base/math.sqrt(n_sample)))

print('Treatment arm utility: ',round(mean_util_treat, 2),'[',round(low_CI_util_treat, 2),',',round(up_CI_util_treat,2),']', '| SD:',round(std_util_treat,2))
print('Base arm utility:      ',round(mean_util_base, 2),'[',round(low_CI_util_base, 2),',',round(up_CI_util_base,2),']', '| SD:',round(std_util_base,2))
print('--------------------------------------------------------------')

#-----------------------------
# Delta arms stats:
#-----------------------------
diff_costs = outputs_treat.iteration_cost_data - outputs_base.iteration_cost_data
mean_cost_diff = np.mean(diff_costs)
std_cost_diff = np.std(diff_costs)
up_CI_cost_diff = mean_cost_diff + (1.96*(std_cost_diff/math.sqrt(n_sample)))
low_CI_cost_diff = mean_cost_diff - (1.96*(std_cost_diff/math.sqrt(n_sample)))

diff_utils = outputs_treat.iteration_util_data - outputs_base.iteration_util_data
mean_util_diff = np.mean(diff_utils)
std_util_diff = np.std(diff_utils)
up_CI_util_diff = mean_util_diff + (1.96*(std_util_diff/math.sqrt(n_sample)))
low_CI_util_diff = mean_util_diff - (1.96*(std_util_diff/math.sqrt(n_sample)))

print('Difference in cost:    ',round(mean_cost_diff, significant_digits),'[',round(low_CI_cost_diff, significant_digits),',',round(up_CI_cost_diff,significant_digits),']', '| SD:',round(std_cost_diff,significant_digits))
print('Difference in utility: ',round(mean_util_diff, significant_digits),'[',round(low_CI_util_diff, significant_digits),',',round(up_CI_util_diff,significant_digits),']', '| SD:',round(std_util_diff,significant_digits))
print('--------------------------------------------------------------')

#---------------------------------------------------------------------------------------------------
# Calculate ICER
#---------------------------------------------------------------------------------------------------
icer, ci = calculate_icer(outputs_base, outputs_treat, calculate_ci=True)
print('\nICER:', round(icer, 2), '(', round(ci[0][0], 2), 'to', round(ci[1][0], 2), ')')

#---------------------------------------------------------------------------------------------------
# Plot on CE plane
#---------------------------------------------------------------------------------------------------
# First need to calculate the differences in cost and utility per-iteration 
# between the treatment arm and base arm
c = outputs_treat.iteration_cost_data - outputs_base.iteration_cost_data
u = outputs_treat.iteration_util_data - outputs_base.iteration_util_data

plt.axhline(0, color='black')
plt.axvline(0, color='black')
plt.ylabel(r'$\Delta$C', rotation=0, fontsize=14, labelpad=14)
plt.xlabel(r'$\Delta$E', fontsize=14)
sns.scatterplot(x=u, y=c)

# Plots some threshold line...
ce_threshold = 20000
x = np.linspace(min(u)-50, max(u)+50, 100)
y = ce_threshold*x
plt.plot(x, y, '-r')

plt.ylim(min(c)-50, max(c)+50) # Add padding to axis limits
plt.xlim(min(u)-0.05, max(u)+0.05)

plt.show()

#----------------------------
# Quadrant counting
# Uses same lambda threshold as detailed in prior step: variable "ce_threshold"
#----------------------------
NW = 0
NE = 0
SW = 0
SE = 0
count_under_threshold = 0
num_iterations = c.shape[0]
for i in range(0, num_iterations):
    if c[i] > 0: # North
        if u[i] > 0: # East
            NE += 1
            y = u[i] * ce_threshold # maximum acceptatble cost at a given value of utility 
            if c[i] < y:
                count_under_threshold += 1
        else: # West
            NW += 1
    else: #South
        if u[i] > 0: # East
            SE += 1
        else: # West
            SW += 1
print()
print('Number of points in each quadrant of CE plane:')
print('North-West:', NW)
print('North-East:', NE, '| Under-threhsold:', count_under_threshold, '| Over-threshold:', NE-count_under_threshold)
print('South-West:', SW)
print('South-East:', SE)
print(NW+NE+SW+SE)
