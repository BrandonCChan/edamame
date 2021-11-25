# monte_carlo_markov
For decision analytic modeling and cost-effectiveness analysis. Allows users to flexibly define model schematics and parameters in excel workbooks before utilizing python packages to load, run, and generate model outputs.

## Requirements (tested with):
See requirements.txt

## 1) Model specification
Specification of model in an excel file (examples in /model_specifications):
The excel workbook contains 4 sheets: 
#### 1) transitions 
    Rows: each row represents a distinct transition between states in the described model
          the properties of the transtion are denoted by the following columns
    Columns: start_state :: the origin state (ie. the A in an A->B transition)
             end_state :: the destination state (ie. the B in an A->B transition)
             type :: specifies the 'type' of transition. Used by the code to handle 
                     the following parameter values when drawing a transtion probability
                     valid entries: time_dependent_weibull, time_dependent_gompertz, 
                                    probabilistic_time_dependent_weibull, probabilistic_time_dependent_gompertz,
                                    beta, gamma, residual, constant

                     time_dependent_weibull | Based on using a weibull distribution. 
                                              Additional params (ie. time) are provided at runtime.
                                              - parameter_1 denotes the regression constant 'const'
                                              - parameter_2 denotes the ancillary parameter 'p'
                     time_dependent_gompertz | Based on using a gompertz distribution. 
                                               Additional params (ie. time) are provided at runtime.
                                               - parameter_1 denotes the regression constant 'const'
                                               - parameter_2 denotes the ancillary parameter 'gamma'
                     probabilistic_time_dependent_weibull | Based on using a weibull distribution. 
                                                            Additional params (ie. time) are provided at runtime.
                                                            - parameter_1 denotes the regression constant 'const'
                                                            - parameter_2 denotes the ancillary parameter 'p'
                                                            - parameter_3 denotes the standard error of 'const'
                                                            - parameter_4 denotes the standard error of 'p'
                     probabilistic_time_dependent_gompertz | Based on using a gompertz distribution.
                                                             Additional params (ie. time) are provided at runtime.
                                                             - parameter_1 denotes the regression constant 'const'
                                                             - parameter_2 denotes the ancillary parameter 'p'
                                                             - parameter_3 denotes the standard error of 'const'
                                                             - parameter_4 denotes the standard error of 'p'
                     beta | parameter_1 denotes number of observed "sucesses", while paramter_2 denotes number of "non-sucesses"
                     gamma | parameter_1 denotes mean, while parameter_2 denotes variance 
                     residual | paramter_1 and parameter_2 hold no meaning, the transition is inputed
                                as the residual
                     constant | parameter_1 denotes a 'locked' transition probability
                     
             parameter_1 :: meaning depends on the specified type (described above)
             parameter_2 :: meaning depends on the specified type (described above)
             parameter_3 :: meaning depends on the specified type (described above)
             parameter_4 :: meaning depends on the specified type (described above)
             notes :: space for additional descriptions you may want to add           
#### 2) costs 
    Rows: Each row represents the cost associated with each state
    Columns: state :: name of the state with associated cost
             type :: method to apply costs. Options include 'static', 'beta' and 'gamma'.
                     beta and gamma utilize the variance and cost (mean) to sample a cost for a given state and iteration
             cost :: cost of being in the corresponding state (assume pre-adjusted to some year)
             cost_variance :: variance of the cost
             notes :: space for additional descriptions you may want to add
#### 3) utilities
    Rows: Each row represents the utility associated with each state
    Columns: state :: name of the state with associated cost
             type :: method to apply utilities. Options include 'static', 'beta' and 'gamma'.
                     beta and gamma utilize the variance and cost (mean) to sample a utility for a given state and iteration
             utility :: cost of being in the corresponding state (assume pre-adjusted to some year)
             utility_variance :: variance of the cost
             notes :: space for additional descriptions you may want to add
#### 4) specification
    Rows: max_iterations :: number of iterations to run
          time_horizon :: maximum amount of time to model (in years)
          cycle_length :: "step size"/length of time of each cycle (days)
          discount_rate :: discount rate to be applied to costs and utilities 
          name_start_state :: which state to initialize the population in (ie. allocate all to initial state)

## 2) Running a model and anlaysis
See Cost-Effectivness Analysis.ipynb in the /notebooks folder for an end-to-end example 

## 3) Running checks
Model Diagnostics.ipynb contains two means of qualitativley checking models. First is the ability to visualize the model specified in an excel file as a graph. Second is to track and visualize the movement of the "population" throughout each model cycle.
