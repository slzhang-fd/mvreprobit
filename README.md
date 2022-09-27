# mvreprobit

Efficient Gibbs sampling procedures for multivariate random effect probit model estimation.

## Description

The `mvreprobit` R package contains code associated with the article 
> Steele, F., Zhang, S., Grundy, E., and Burchardt, T. (2022). Longitudinal analysis of exchanges of support between parents and children in the UK.

It contains tailored Gibbs sampling procedures for multivariate random effect probit modelling.
It currently includes the following models fitted in the paper: a 2-level model
(1 random effect) and a 3-level model. The latter is an extension of the standard
3-level hierarchical model to include a time-varying level 3 random effect, giving
3 random effects in total. The code is concise, self-explanatory, and can be extended 
to any number of response variables and random effects of arbitrary multilevel designs.

## Installation

See [Wiki](https://github.com/slzhang-fd/mvreprobit/wiki) page

## Demo examples

[1. Four-process 2-level (1 random effect) probit model](https://github.com/slzhang-fd/mvreprobit/wiki/1.-Four-processes-one-random-effect-probit-model)

[2. Four-process 3-level (3 random effects) probit model, with additional time-varying level 3 random effect](https://github.com/slzhang-fd/mvreprobit/wiki/2.-Four-processes-three-random-effects-probit-model)

### Libraries and dependencies used by the code

-   truncnorm
-   mvtnorm
-   MCMCpack
-   matrixcalc

### License

GPL v3.0
