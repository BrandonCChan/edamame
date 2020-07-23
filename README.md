# monte_carlo_markov
For decision analytic modeling and cost-effectiveness analysis

## Requirements (tested with):
Core:
* python 3.7
* numpy 1.18.5
* pandas 0.25.3
* openpyxl 3.0.1 
* arch 4.15 
* matplotlib 3.1.3 
* seaborn 0.9.0

For model results tracing:
* plotly 4.7.1
* chart-studio 1.1.0

For model visualiztion:
* dash 1.12.0 
* dash-cytoscape 0.1.1
* jupyter-dash 0.2.1.post1 
* networkx 2.4

* jupyter 1.0.0 (for running notebooks)

## 1) Model specification
Specification of model in excel file (examples in /model_specifications):
The excel workbook contains 3 sheets: 
1) transitions 
    Rows: each row represents a distinct transition between states in the described model
          the properties of the transtion are denoted by the following columns
    Columns: start_state :: the origin state (ie. the A in an A->B transition)
             end_state :: the destination state (ie. the B in an A->B transition)
             type :: specifies the 'type' of transition. Used by the code to handle 
                     the following parameter values when drawing a transtion probability
                     valid entries: time-dependant, beta, gamma, residual, constant

                     time-dependant | Based on using a weibull distribution. 
                                      Additional params (ie. time) are provided at runtime.
                                      - parameter_1 denotes the *weibull shape parameter 'p'
                                      - parameter_2 denotes the regression constant 'c'
                     beta and gamma | parameter_1 denotes mean, while paramter_2 denotes variance
                     residual | paramter_1 and parameter_2 hold no meaning, the transition is inputed
                                as the residual
                     constant | parameter_1 denotes a 'locked' transition probability
             
             parameter_1 :: meaning depends on the specified type (described above)
             parameter_2 :: meaning depends on the specified type (described above)
2) costs 
    Rows: Each row represents the cost associated with each state
    Columns: state :: name of the state with associated cost
             cost :: cost of being in the corresponding state (assume pre-adjusted to some year)
             cost_variance :: variance of the cost
3) utilities
    Rows: Each row represents the utility associated with each state
    Columns: state :: name of the state with associated cost
             utility :: cost of being in the corresponding state (assume pre-adjusted to some year)
             utility_variance :: variance of the cost
4) specification
    Rows: max_iterations :: number of iterations to run
          time_horizon :: maximum amount of time to model (in years)
          cycle_length :: "step size"/length of time of each cycle (days)
          discount_rate :: discount rate to be applied to costs and utilities 
          name_start_state :: which state to initialize the population in (ie. allocate all to initial state)

## 2) Running a model and anlaysis
See Cost-Effectivness Analysis.ipynb in the /notebooks folder for an end-to-end example 

## 3) Running checks
Model Diagnostics.ipynb contains two means of qualitativley checking models. First is the ability to visualize the model specified in an excel file as a graph. Second is to track and visualize the movement of the "population" throughout each model cycle.