#---------------------------------------------------------------------------------------------------
# Helper functions for cost-effectivness analysis modeling - Monte Carlo Markov Chain
# Brandon Chan - July 2020
#---------------------------------------------------------------------------------------------------
# Import packages
#---------------------------------------------------------------------------------------------------
import pandas as pd 
import numpy as np 
import math 

#---------------------------------------------------------------------------------------------------
# File IO
#---------------------------------------------------------------------------------------------------
def read_file(filepath):
    model_specification_file = pd.ExcelFile(filepath)
    check_excel_file(model_specification_file) # Run checks
    return model_specification_file

def check_excel_file(excel_book):
    '''
    Excel file Q/A - checking for specification errors
    If error is detected, an exception will be raised (printed to console) and the program will terminate
    '''
    # Separate each sheet in the excel book
    transitions_df = pd.read_excel(excel_book, 'transitions')
    costs_df = pd.read_excel(excel_book, 'costs')
    utilities_df = pd.read_excel(excel_book, 'utilities')
    specification_df = pd.read_excel(excel_book, 'specification', header=None, index_col=0) # TODO:// add specification checks

    # Specification of variables regarding states in model
    start_state_names = transitions_df['start_state'].unique().tolist()
    end_state_names = transitions_df['end_state'].unique().tolist()
    unique_states = list(set(start_state_names+end_state_names))
    
    # Check only one unique start-end state pair is defined. (ie. no repeats/multiples of the same transition)
    check_multiple = transitions_df.groupby(['start_state','end_state']).size().reset_index().rename(columns={0:'count'})
    if (check_multiple['count'] > 1).any():
        print(check_multiple.loc[check_multiple['count'] > 1])
        raise ValueError('Multiple identical defined transitions found. Please check transitions sheet in input excel document')

    # Checks that each state has an entry/exit
    # TODO: consider that its possible for the entry/start state to be non-reenterable
    #       with that in mind its more logical that all states must have an outbound rather than enforcing entry
    #if np.setdiff1d(start_state_names, end_state_names).shape[0] > 0:
    #    raise ValueError('State missing exit:',np.setdiff1d(start_state_names, end_state_names).tolist(),'please check transitions sheet in excel document')
    if np.setdiff1d(end_state_names, start_state_names).shape[0] > 0:
        raise ValueError('State missing entry:',np.setdiff1d(end_state_names, start_state_names).tolist(),'please check transitions sheet in excel document')
            
    # Check for missing costs 
    if np.setdiff1d(unique_states, costs_df['state'].tolist()).shape[0] > 0:
        raise ValueError('Missing cost value for state(s):',np.setdiff1d(unique_states, costs_df['state'].tolist()),'please check cost sheet in excel document')

    # Check for missing utilities
    if np.setdiff1d(unique_states, utilities_df['state'].tolist()).shape[0] > 0:
        raise ValueError('Missing utility value for state(s):',np.setdiff1d(unique_states, utilities_df['state'].tolist()),'please check utilities sheet in excel document')


#---------------------------------------------------------------------------------------------------
# Helper functions for selecting transition probabilities
#---------------------------------------------------------------------------------------------------
# Not sure how to use yet... but ok.
def get_dirchlet(parameters):
    #Samples a value from the dirchlet distribution based on input parameters
    return np.random.dirichlet(parameters,1)

def get_gamma(mean, variance):
    '''
    Samples a value from the gamma distribution based on an input mean and variance
    Inputs: mean = mean
            variance = variance
    Output: A float between 0 and 1 sampled from the beta distribution defined
            by the input parameters
    '''
    # TODO: Error checking (look at wikipedia too) - consider moving out when checking the spreasheet. (all should be specified correctly before running)
    # mean-var > 1?
    # Actually used for sampling utility or cost?
    alpha = (mean**2)/(variance**2)  
    beta = (variance**2)/mean
    return np.random.gamma(alpha, beta)

def get_beta(mean, variance):
    '''
    Samples a value from the beta distribution based on an input mean and variance
    Inputs: mean = mean
            variance = variance
    Output: A float between 0 and 1 sampled from the beta distribution defined
            by the input parameters
    '''
    # TODO: Error checking (look at wikipedia too) - consider moving out when checking the spreasheet. (all should be specified correctly before runnin
    # if mean - close to 1 and 0 may return errors
    # if variance is > mean - maybe issues
    # Actually used for sampling utility or cost?
    alpha = mean*((mean*(1-mean)/variance**2) - 1) #(((1-mean)/variance) - (1/mean)) * mean**2  
    beta = (1-mean)*(mean/variance**2*(1-mean) - 1) #alpha * ((1/mean) - 1)
    return np.random.beta(alpha, beta)

def get_time_dependant(p, const, time, cycle_length):
    '''
    Obtains the transition probaility for a time-dependant transition by sampling the approximated function
    at the given time interval. 
    Inputs: time = time in model, so i in 1:ncycle
            cycle_length = cycle length, in days 
            p = weibull shape parameter from regression. 
            const = constant in regression
    Output: A float between 0 and 1 denoting the time-dependant transition probability from A to B 
            based on the input parameters
    '''
    # TODO:// if other survival curves are being used implment separate return functions? (or add to main "switch")
    return 1-math.exp((math.exp(const))*(((time*cycle_length)-cycle_length)**p)-((math.exp(const))*((time*cycle_length)**p)))

def set_transition(transition_type, params):
    '''
    Function that serves as a switch for a number of transition-probability retriever functions
    Inputs: transition_type = string that denotes which function to redirect params to
            params = a list of values that get passed to a secondary function and will dictate
                     the output
    Output: transition probability sampled/assigned according to the transition_type input
    '''
    #if transition_type == 'beta':
    #    return get_beta(params[0], params[1])
    if transition_type == 'beta':
        return np.random.beta(params[0], params[1])
    elif transition_type == 'gamma':
        return get_gamma(params[0], params[1]) 
    elif transition_type == 'time-dependant':
        return get_time_dependant(params[0], params[1], params[2], params[4])
    elif transition_type == 'constant':
        return params[0]
    else:
        raise ValueError('Invalid transition type provided:',str(transition_type))

#---------------------------------------------------------------------------------------------------
# Helper functions related to operations on the transition matrix
#---------------------------------------------------------------------------------------------------
def check_row_sums(matrix):
    '''
    Function that serves soley as a check that all outbound probabilities sum to 1
    Input: matrix = transition matrix of dimensions [num_states x num_states]
    Output: None. Will throw an error and terminate runtime if condition not met
    '''
    row_sums = np.round(matrix.sum(axis=1), 5) # rounding to 5 significant digits for tolerance for close to 1 values
    if np.any(row_sums != 1):
        print(matrix)
        print(row_sums)
        raise ValueError('Error: transitions do no add to 1. Suggest using normalize_transitions()...')

def calculate_residual(matrix, state_index):
    '''
    Function that calculates the residual of all states in the matrix passed
    Inputs: matrix = n x n numpy array representing the transition matrix of the model
            state_index = an integer representing the row of the state that we wish to calculate the residual of
    Output: residual of the outbound transitions
    '''
    row_sum = matrix[state_index, :].sum(axis=1)
    residual = 1 - row_sum
    print(residual)
    return residual
    

def normalize_transitions(matrix):
    '''
    Function that rescales/normalizes the outbound transtions (row-wise) to have a sum of 1
    Input: matrix = transition matrix of dimensions [num_states x num_states]
    Output: matrix rescaled to have each row sum to 1 and each outbound transition scaled appropriately 
    '''
    row_sums = matrix.sum(axis=1)
    return matrix / row_sums[:, np.newaxis] # need newaxis to "reshape" the array so its properly applied

#---------------------------------------------------------------------------------------------------
# Other helper functions
#---------------------------------------------------------------------------------------------------
def check_model_population(pop):
    '''
    Function to check that the entire popoulation is still in the model (ie. no disappearances)
    Input: pop = [1 x number_of_states] numpy array 
    Output: None - raises error if does not sum to 1
    '''
    if round(pop.sum(),5) != 1:
        raise ValueError('Error: proportions do not sum to 1. Loss or gain of population.', round(pop.sum(),5))