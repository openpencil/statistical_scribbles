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

Gaussians
* Adding finite fluctuations results in a distribution of sums that have little extraneous information about the underlying data-generating process (DGP) aside from mean and spread. Thus statistical models based on Gaussians cannot reliably identify micro-process of the DGP.
* When all we know or are willing to say about the distributions is their mean and variance, then Gaussians are the most consistent with our assumptions. The Gaussian is the most natural expression of our state of ignorance.
* Gaussians have the least surprising and least informative assumption to make.
* The (x-u)^2: quadratic shape. Exponentiating it gives the bell shape. The rest is scaling and standardizing so it sums to one.
* 1/(sigma)^2 = Precision = _tau_
