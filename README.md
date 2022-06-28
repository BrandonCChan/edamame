[![edamame](https://github.com/BrandonCChan/edamame/blob/master/img/edamame.png)](https://github.com/BrandonCChan/edamame)


Economic Decision Analytic MArkov Model Evaluation (EDAMAME)

A framework and tools for decision analytic modeling and cost-effectiveness analysis written in Python. Mainly intended for use in the field of health economics to conduct economic evaluations using mote-carlo markov model simulations. Allows users to flexibly define model schematics and parameters in excel workbooks before utilizing edamame code to load, run, and generate model outputs.

## Requirements:
See requirements.txt for more details

Reflects development environment but is possible that older and newer versions of packages may still work.

## Running a cost-effectivness analysis
See [Cost-Effectivness Analysis.ipynb](https://github.com/BrandonCChan/edamame/blob/master/notebooks/Cost-Effectiveness%20Analysis.ipynb) for an end-to-end example with annotations and explainations!

Alternativley, [cost_effectivness_analysis.py](https://github.com/BrandonCChan/edamame/blob/master/examples/cost_effectivness_analysis.py) provides a script version of most of the code within the example jupyter notebook.

## Defining and specifying a model
The structure (i.e. states and transitions), associated state costs, associated state utilties, and simulation paramters (i.e. number of iterations, time-horizon) of a model are defined as a formatted excel document. Each model arm (comparator) is defined by it's own excel workbook. Templates of a simple 3 state model are provided in [here](https://github.com/BrandonCChan/edamame/tree/master/model_specifications).

See the [wiki page](https://github.com/BrandonCChan/edamame/wiki/Model-Specification) for more details.

## Other resources and examples

### Model disgnostics and debugging
[Model Diagnostics.ipynb](https://github.com/BrandonCChan/edamame/blob/master/notebooks/Model%20Diagnostics.ipynb) contains two means of qualitativley checking models. First is the ability to visualize the model specified in an excel file as a graph. Second is to track and visualize the movement of the "population" throughout each model cycle.

### Univariate sensitivity analysis
[Univariate Sensitivity Analysis.ipynb](https://github.com/BrandonCChan/edamame/blob/master/notebooks/Univariate%20Sensitivity%20Analysis.ipynb) provides an example of conducting two types of univariate sensitivity analysis using previously saved model outputs. 1) Is the adjustment of the cost of a state; 2) Is the adjustment of a transition probability.

### Example of time-varying transition probabilities
[Testing with two-state model](https://github.com/BrandonCChan/edamame/blob/master/notebooks/Testing%20with%20two-state%20model.ipynb) demonstrates the use of time-varying transition probabilities (i.e. a weibull) with a simple 2-state model. This example is useful in testing model behavior and validating the parameterization of a survivial regression you may have run with external data. 

## Documentation
See the relevant [wiki pages](https://github.com/BrandonCChan/edamame/wiki/Technical-Documentation)
