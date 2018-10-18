# Code for: How robusts are Structural Equation Models to model mis-specification

[![DOI](https://zenodo.org/badge/153457961.svg)](https://zenodo.org/badge/latestdoi/153457961)

An arXiv version of the manuscript can be found [here](https://arxiv.org/abs/1803.06186). 

This manuscript explore how different model mis-specification affect various SEM fitness metrics.

In this repo is stored the code used to generate the figures shown in the manuscript. There are two R files:

* simSEM.R: contain the code to run the simulations and generate the figures

* simSEM_function.R: contain all the functions used to run the simulations

The parameters that can be varied are:

* Type: the type of data-generation, one character from the following, random, exact, complex, simple or shuffled.
* X: the number of covariates in the models
* N: the sample size
* C: the graph connectance, how many links are drawn between the variables
* lv: whether to fit the model via lavaan (TRUE) or piecewiseSEM (FALSE)
* sd_res: the standard deviation of the residuals, ie the amount of white noise
* sd_eff: the standard deviation of the path coefficients, ie the amount of signal, larger values means that larger path coefficients are more likely to be generated.
