# CAIT
Code for Covariate Adjusted Interaction Trees


The file CodeUsedMainSimulations.R implments the simulations presented in the manuscript

**`Main.R`**

Implements the simulations presented in the manuscript using the following functions.

Functions | Description
----------|------------
**`CvMethod1.R`** | Contains functions implementing final tree selection method 1 for different estimators. `.EstNs` in the function names stands for node specific estimators; `.Est1` stands for model standardization estimators; `.Est2` stands for data adaptive estimators; `.TruePaper` stands for true model is fitted in the paper’s simulation setting (same below).
**`CvMethod2.R`** | Contains functions implementing final tree selection method 2 for different estimators. 
**`EstCondEff.R`** | Contains different functions to estimate the conditional expectation of the outcome `Y` for a given treatment adjusted for covariates, used for the data adaptive (DA) estimator.
**`EstDaTempFunc.R`** | Contains different functions to implement splitting/evaluation functions that uses data adaptive (DA) estimators. `.da` in the function names stands for data adaptive; `.b` is added to function names when the outcome is binary (same below).
**`EstMsTempFunc.R`** | Contains different functions to implement splitting/evaluation functions that uses model standardization (MS) estimators. `.ms` in the function names stands for model standardization.
**`EstNsTempFunc.R`** | Contains different functions to implement splitting/evaluation functions that uses node specific (NS) means. `.ns` in the function names stands for node specific means.
**`EvalMeas.R`** | Contains function to evaluate the final tree through multiple criteria.
**`iTemp.R`** | Initialization function for `rpart()`.
**`library.R`** | Contains code to install uninstalled and load libraries. 
**`MobVt.R`** | Implements MOB and VT algorithm simulations in the manuscript.
**`pruning.R`** | Contains functions to calculate the sequence of candidate trees.
**`SimData.R`** | Simulate data under different settings. 
