#---------------------------------------------------------------------------------------------------
# Example code for how to conduct a cost effectivness analysis using
# the framework. Mirrors code detailed in "Cost-Effectivness Analysis.ipynb"
# Please see the jupyter notebook for additional code on calculating other output metrics
#
# Brandon Chan | January 2021
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
print()
print('mean cost treatment arm:', round(np.mean(outputs_treat.iteration_cost_data), significant_digits))
print('mean cost base case arm:', round(np.mean(outputs_base.iteration_cost_data), significant_digits))
print('mean utility treatment arm:', round(np.mean(outputs_treat.iteration_util_data), significant_digits))
print('mean utility base case arm:', round(np.mean(outputs_base.iteration_util_data), significant_digits))
print('------------------------')
print('std cost treatment arm:', round(np.std(outputs_treat.iteration_cost_data), significant_digits))
print('std cost base case arm:', round(np.std(outputs_base.iteration_cost_data), significant_digits))
print('std utility treatment arm:', round(np.std(outputs_treat.iteration_util_data), significant_digits))
print('std utility base case arm:', round(np.std(outputs_base.iteration_util_data), significant_digits))
print('------------------------')
print('Difference in costs:', round((np.mean(outputs_treat.iteration_cost_data) - np.mean(outputs_base.iteration_cost_data)), significant_digits))
print('Difference in utility:', round((np.mean(outputs_treat.iteration_util_data) - np.mean(outputs_base.iteration_util_data)), significant_digits))

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
