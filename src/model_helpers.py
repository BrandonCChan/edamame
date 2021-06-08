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

    # Check dirichlet parameters... All outbound for a given state must be dirichlet?
    dirichlet_transitions = transitions_df.loc[transitions_df.type == 'dirichlet']
    for start_state in dirichlet_transitions.start_state.unique():
        outbound_transitions = transitions_df.loc[transitions_df.start_state == start_state]
        if outbound_transitions.loc[outbound_transitions.type != 'dirichlet'].shape[0] > 0:
            raise ValueError('Not all outbound transitions for state:', start_state ,' are of type dirichlet. Please check transitions sheet in input excel document')

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

    # Check specification parameters
    if int(specification_df.loc['max_iterations'].values[0]) < 1:
        raise ValueError('Invalid number of model iterations. Needs to be greater than or equal to 1')
    if specification_df.loc['cycle_length'].values[0] > specification_df.loc['time_horizon'].values[0] * 365:
        raise ValueError('Invalid cycle length or time horizon. Cycle length must be smaller than horizon')
    if specification_df.loc['name_start_state'].values[0] not in unique_states:
        raise ValueError('Start state must be a valid defined state in "transitions" sheet. Please check. Invalid entry:',specification_df.loc['name_start_state'].values[0])
    if specification_df.loc['discount_rate'].values[0] < 0 or specification_df.loc['discount_rate'].values[0] > 1:
        raise ValueError('Invalid Discount Rate. Must be represented as a number between 0 and 1.')


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
    Output: A float between 0 and 1 sampled from the gamma distribution defined
            by the input parameters

    See: https://wiki.analytica.com/index.php?title=Gamma_distribution#:~:text=To%20estimate%20the%20parameters%20of,)%2FMean(X%2C%20I)&text=alpha%20%3A%3D%204%2FSkewness(X%2C%20I)%5E2
         https://www.itl.nist.gov/div898/handbook/eda/section3/eda366b.htm
    '''
    # TODO: Error checking (look at wikipedia too) - consider moving out when checking the spreasheet. (all should be specified correctly before running)
    # Currently based on criteria from David's old code
    # if variance > mean or mean < 10 or variance < 1.1:
    #    raise ValueError('Invalid mean and/or variance parameters for gamma prameter estimation')
    
    alpha = (mean**2)/variance  
    beta = variance/mean
    return np.random.gamma(alpha, beta)

def get_beta(mean, variance):
    '''
    Samples a value from the beta distribution based on an input mean and variance
    Inputs: mean = mean
            variance = variance
    Output: A float between 0 and 1 sampled from the beta distribution defined
            by the input parameters

    See: https://en.wikipedia.org/wiki/Beta_distribution#Mean_and_sample_size
    '''
    # TODO: Error checking (look at wikipedia too) - consider moving out when checking the spreasheet. (all should be specified correctly before running)
    # if mean - close to 1 and 0 may return errors
    if variance > mean*(1-mean):
        raise ValueError('Variance is too large to estimate parameters alpha and beta. Mean:', mean, 'Var:',variance)
    
    alpha = mean*(((mean*(1-mean))/variance) - 1)  
    beta = (1-mean)*(((mean*(1-mean))/variance) - 1)
    return np.random.beta(alpha, beta)

def get_beta_rate_to_prob(a, b, t, cycle_length):
    '''
    Sample a proportion from a beta distribution and convert to a rate then probability based on cycle length 
    using formulas outlined in the "green book" by Briggs et al.
    Inputs: a = the "numerator"
            b = the "denominator"
            t = the time span in years that the events a/b were observed
            cycle_length = cycle length of the model in days
    Ouputs: a float between 0 and 1 
    '''
    beta_p = np.random.beta(a, b-a)
    rate = -(math.log(1-beta_p))/t
    probability = 1 - math.exp(-rate * (cycle_length/365))
    if probability > 1 or probability < 0:
        raise ValueError('Invalid probability calculated >1 or <0. Please check params:', a, b, t, probability)
    return probability

def get_time_dependent_weibull(const, p, cycle, cycle_length):
    '''
    Obtains the transition probability for a time-dependent transition by sampling the approximated function 
    at the given time interval. 
    
    Implementation of formula: tp(t_u) = 1 - S(t)/S(t-u) from "Decision Modelling for Health Economic Evaluation" by Briggs et.al
    Based on survival functions outlined here: https://www.stata.com/manuals13/ststreg.pdf (Page 6.)
    
    Inputs: const = regression constant (as output from a PH fitted weibull in stata)
            p = paramter p (as output from a PH fitted weibull in stata)
            cycle = the current cycle of the model (ie. cycle i in 1:max_num_cycles)
            cycle_length = the length of a cycle in days. 
    Output: tdtp = A float between 0 and 1 denoting the time-dependent transition probability from A to B 
            based on the input parameters

    Note* : Assumes weibull fitted to survival curve at a time scale of YEARS on x-axis.
    '''
    lmbda = math.exp(const)

    # adjusts to yearly x-axis. subtract 1 from t1 as model "starts at time 1" however first transition calculation is based off of t0 and t1
    t1 = ((cycle-1)*cycle_length) / 365
    t2 = ((cycle)*cycle_length) / 365
    
    tdtp = 1 - ((math.exp(-lmbda*(t2**p))) / (math.exp(-lmbda*(t1**p))))

    if tdtp > 1 or tdtp < 0:
        raise ValueError('Transition sampled is greater than 1 or less than 0. Sampled value:',round(tdtp,4),'at cycle:',cycle)
    
    return tdtp

def get_time_dependent_gompertz(const, gamma, cycle, cycle_length):
    '''
    Obtains the transition probability for a time dependent transition represented as a gompertz
    
    Implementation of formula: tp(t_u) = 1 - S(t)/S(t-u) from "Decision Modelling for Health Economic Evaluation" by Briggs et.al
    Based on survival functions outlined here: https://www.stata.com/manuals13/ststreg.pdf (Page 6.)
    
    Inputs: const = the regression constant (as output from a PH fitted gompertz in stata)
            gamma = parameter gamma (as output from a PH gompertz fitted in stata)
            cycle = the current cycle of the model (ie. cycle i in 1:max_num_cycles)
            cycle_length = the length of a cycle in days.
    Output: tdtp = A float between 0 and 1 denoting the time-dependent transition probability from A to B 
            based on the input parameters

    Note** : Assumes gompertz fitted to survival curve at a time scale of YEARS on x-axis.
    '''
    lmbda = math.exp(const) 

    t1 = ((cycle-1)*cycle_length) / 365
    t2 = ((cycle)*cycle_length) / 365

    # adjusts to yearly x-axis. subtract 1 from t1 as model "starts at time 1" however first transition calculation is based off of t0 and t1
    tdtp = 1 - ((math.exp(-lmbda*gamma**-1*(math.exp(gamma*t2)-1))) / (math.exp(-lmbda*gamma**-1*(math.exp(gamma*t1)-1))))

    if tdtp > 1 or tdtp < 0:
        raise ValueError('Transition sampled is greater than 1 or less than 0. Sampled value:',round(tdtp,4),'at cycle:',cycle)
    
    return tdtp

def calculate_residual(matrix, state_index):
    '''
    Function that calculates the residual of all states in the matrix passed
    Inputs: matrix = n x n numpy array representing the transition matrix of the model
            state_index = an integer representing the row of the state that we wish to calculate the residual of
    Output: residual of the outbound transitions
    '''
    row_sum = matrix[state_index, :].sum()
    residual = 1 - row_sum
    
    # Error checking 
    # TODO: Think about how you want these errors to be raised.. Currently assumes if < 1, residual == 0
    if residual < 0:
        return 0
        raise ValueError('Error: residual is negative', round(residual,5), 'state_index:', state_index)
        
    return residual

def set_transition(transition_type, **kwargs):
    '''
    Function that serves as a switch for a number of transition-probability retriever functions
    Inputs: transition_type = string that denotes which function to redirect params to
            kwargs = a number of keyword arguments provided to the function. Different transition types require different
                     arguments. See code for details. Error is raised when incorrect arguments are provided for a transition type.
    Output: transition probability sampled/assigned according to the transition_type input
    '''
    if transition_type == 'beta':
        if not all (parameter in kwargs for parameter in ('a','b')):
            raise ValueError('Incorrect inputs specified for beta. Need a and b.')
        return np.random.beta(kwargs['a'], kwargs['b']-kwargs['a'])
    elif transition_type == 'beta_r2p':
        if not all (parameter in kwargs for parameter in ('a','b','t','cycle_length')):
            raise ValueError('Incorrect inputs specified for beta_rate. Need a and b and t and cycle_length')
        return get_beta_rate_to_prob(kwargs['a'], kwargs['b'], kwargs['t'], kwargs['cycle_length'])
    elif transition_type == 'gamma':
        if not all (parameter in kwargs for parameter in ('a','b')):
            raise ValueError('Incorrect inputs specified for gamma. Need a and b.')
        return get_gamma(kwargs['a'], kwargs['b']) 
    elif transition_type == 'time_dependent_weibull':
        if not all (parameter in kwargs for parameter in ('const','ancillary','cycle','cycle_length')):
            raise ValueError('Incorrect inputs specified for time dependent. Need const, p, cycle, and cycle_length.')
        return get_time_dependent_weibull(kwargs['const'],kwargs['ancillary'],kwargs['cycle'],kwargs['cycle_length'])
    elif transition_type == 'time_dependent_gompertz':
        if not all (parameter in kwargs for parameter in ('const','ancillary','cycle','cycle_length')):
            raise ValueError('Incorrect inputs specified for time dependent. Need const, gamma, cycle, and cycle_length.')
        return get_time_dependent_gompertz(kwargs['const'],kwargs['ancillary'],kwargs['cycle'],kwargs['cycle_length'])
    elif transition_type == 'constant':
        if 'transition' not in kwargs:
            raise ValueError('Incorrect inputs specified for constant:', kwargs,'Only need 1: transition')
        return kwargs['transition']
    elif transition_type == 'residual':
        if not all (parameter in kwargs for parameter in ('transition_matrix','update_index')):
            raise ValueError('Incorrect inputs specified for residual calculation. Need transition_matrix and update_index')
        return calculate_residual(kwargs['transition_matrix'], kwargs['update_index'])
    else:
        raise ValueError('Invalid transition type provided:',str(transition_type))

#---------------------------------------------------------------------------------------------------
# Helper functions related to operations on the transition matrix
#---------------------------------------------------------------------------------------------------
def check_row_sums(matrix):
    '''
    Function that serves solely as a check that all outbound probabilities sum to 1
    Input: matrix = transition matrix of dimensions [num_states x num_states]
    Output: None. Will throw an error and terminate runtime if condition not met
    '''
    row_sums = np.round(matrix.sum(axis=1), 5) # rounding to 5 significant digits for tolerance for close to 1 values
    if np.any(row_sums != 1):
        print(matrix)
        print(row_sums)
        raise ValueError('Error: transitions do no add to 1. Suggest using normalize_transitions()...')
    
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
    Function to check that the entire population is still in the model (ie. no disappearances)
    Input: pop = [1 x number_of_states] numpy array 
    Output: None - raises error if does not sum to 1
    '''
    if round(pop.sum(),5) != 1:
        raise ValueError('Error: proportions do not sum to 1. Loss or gain of population.', round(pop.sum(),5))