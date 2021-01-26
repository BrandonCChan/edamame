#---------------------------------------------------------------------------------------------------
# Dash app to plot model structure. Separated from notebook code in case jupyter dash is too
# much of an issue to set up compared to base Dash/Plotly
# Brandon Chan | January 2021
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# Boilerplate
#---------------------------------------------------------------------------------------------------
import pandas as pd
import numpy as np
import networkx as nx
import dash
import dash_core_components as dcc
import dash_cytoscape as cyto
import dash_html_components as html
import sys
sys.path.insert(0,"../src/")
from markov_modeling import *
from model_helpers import set_transition, calculate_residual, normalize_transitions, check_row_sums

#---------------------------------------------------------------------------------------------------
# Load model from specification
#---------------------------------------------------------------------------------------------------
model_specification = ModelSpec('../model_specifications/test_parameters_base.xlsx', 
                                model_name='test_base')

#---------------------------------------------------------------------------------------------------
# Generate matrix representation of model
#
# Dimensions of matrix = [num_states x num_states]
# Rows map to starting state, columns map to target state
#---------------------------------------------------------------------------------------------------
transitions_df = model_specification.structure
num_states = model_specification.num_states
state_mapping = model_specification.state_mapping
cycle_length = model_specification.cycle_length
cycle = 1

#specify empty transition matrix
transition_matrix = np.zeros((num_states,num_states)) 

# Keep track of specific indicies in the transition matrix that need to be updated "in-simulation"
# ie. time-dependant transitions and transitions that get resampled every iteration
# Intended to be stored as a list of dictionaries
resample_indicies = []
time_dependent_indicies = []
probabilistic_time_dependent_indicies = []
residual_indicies = []

# Iterate through specified transtions and initialize values of matrix 
# Log indicies of transitions that need to be updated to quickly index the correct position in the transition matrix
for t in transitions_df.itertuples():
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
        elif t_type in ['time_dependent_weibull', 'time_dependent_gompertz']:
            time_dependent_indicies += [{'start_state':t[1], 'end_state':t[2], 
                                         'i':start_state_index, 'j':end_state_index, 
                                         'type':t_type, 
                                         'const':params[0], 'ancillary':params[1]}]
        elif t_type in ['probabilistic_time_dependent_weibull', 'probabilistic_time_dependent_gompertz']:
            # Is treated the same as a time-dependent weibull or gompertz, hence why the type in the dictionary is adjusted to chop off the "probabilistic" part
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

#------------------------------------------------------------------------------------------------------
# Assign all transitions (initial)
#------------------------------------------------------------------------------------------------------
for t in resample_indicies:
    transition_matrix[t['i'],t['j']] = set_transition(t['type'], a=t['a'], b=t['b']) 

for t in time_dependent_indicies:
    transition_matrix[t['i'],t['j']] = set_transition(t['type'], const=t['const'], ancillary=t['ancillary'], cycle=cycle, cycle_length=cycle_length)

for t in probabilistic_time_dependent_indicies:
    transition_matrix[t['i'],t['j']] = set_transition(t['type'], const=t['const'], ancillary=t['ancillary'], cycle=cycle, cycle_length=cycle_length)
    
# Calculate residual transition probailities if needed - Update matrix as appropriate
for t in residual_indicies:
    transition_matrix[t['i'], t['j']] = 0 # AFTER RESAMPLING RESIDUAL IS RECALCULATED FROM ZERO (but only zero specific transition i to j)
    transition_matrix[t['i'], t['j']] = calculate_residual(transition_matrix, t['i'])

# Calculate residual transition probailities if needed - Update matrix as appropriate
for t in residual_indicies:
    transition_matrix[t['i'], t['j']] = 0 # AFTER RESAMPLING RESIDUAL IS RECALCULATED FROM ZERO (but only zero specific transition i to j)
    transition_matrix[t['i'], t['j']] = calculate_residual(transition_matrix, t['i'])

# rescale row-wise to ensure row sums (ie. all transitions out of a state sum to 1)
transition_matrix = normalize_transitions(transition_matrix) 
check_row_sums(transition_matrix)

#---------------------------------------------------------------------------------------------------
# Convert to graph and initialize app
#---------------------------------------------------------------------------------------------------
transition_matrix_pd = pd.DataFrame(transition_matrix, index=state_mapping, columns=state_mapping)
G = nx.from_pandas_adjacency(transition_matrix_pd, create_using=nx.DiGraph)
pos = nx.circular_layout(G)

cy = nx.readwrite.json_graph.cytoscape_data(G)
for node in G.nodes:
    G.nodes[node]['position'] = {'x': pos[node][0], 'y': pos[node][1]}

nodes = [{'data': {'id':node, 'label':node},
          'position': {'x': pos[node][0]*250, 'y':pos[node][1]*250}
         } for node in G.nodes]

edges = cy['elements']['edges']
converted_graph = nodes + edges

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div([
    cyto.Cytoscape(
        id = 'model-graph',
        layout = {'name':'preset'},
        style = {'width':'100%', 'height':'800px'},
        elements=converted_graph,
        stylesheet = [
            {
                'selector':'node',
                'style':{'content':'data(label)',
                         'text-halign':'center',
                         'text-valign':'center',
                         'width':'label',
                         'height':'label',
                         'shape':'ellipse',
                         'padding':'10px'}
            },
            {
                'selector':'edge',
                'style':{'curve-style':'bezier',
                         'target-arrow-shape': 'triangle',
                         'label':'data(weight)',
                         'loop-sweep':'-30deg'}
            }
        ])
])

if __name__ == '__main__':
    app.run_server(debug=True)
