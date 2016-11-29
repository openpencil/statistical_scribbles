# Rethinking concepts

Concepts for recapitulating

* Regularizing or weakly informative priors --> Penalized likelihood
* A "subjective Bayesian"
* Choosing priors --> Choosing estimators,likelihood penalties and loss functions
* Combinations of data, likelihood, parameters and priors determine the estimates of the posterior distribution
* Expectation of data over the prior is called the marginal likelihood

Integration
* If the parameters are discrete, the marginal is summed over the values of the parameters. On the other hand, if parameters are continuous, the marginal is integrated over the values of the parameter.
* To compute the marginal one must average over the uncertainty in all other parameters, when describing the uncertainty in the value of the parameter of interest. e.g. if the parameter of interest is the result of a test for HIV status, one must average the positive tests over HIV+ as well as HIV-, where HIV status is a parameter with uncertainty. This is a tough integral, but with sampling, it is only a frequency format or natural frequencies.
* Samples can be used to summarize and simulate model output.
