# Compute proportion in the survey sample (standard estimator)

Proportion estimated using the survey sample and confidence intervals
based on the Clopper-Pearson and the standard asymptotic approach.

## Usage

``` r
survey_mle(R, n, pi0 = 0, alpha = 0, beta = 0, gamma = 0.05, V = NULL, ...)
```

## Arguments

- R:

  A `numeric` that provides the people of positive people in the sample.

- n:

  A `numeric` that provides the sample size.

- pi0:

  A `numeric` that provides the prevalence or proportion of people (in
  the whole population) who are positive, as measured through a
  non-random, but systematic sampling (e.g. based on medical selection).
  Default value is `0` and in this case this information is not used in
  the estimation procedure.

- alpha:

  A `numeric` that provides the False Negative (FN) rate for the
  sample R. Default value is `0`.

- beta:

  A `numeric` that provides the False Positive (FP) rate for the
  sample R. Default value is `0`.

- gamma:

  A `numeric` that is used to compute a (1 - gamma) confidence region
  for the proportion. Default value is `0.05`.

- V:

  A `numeric` that corresponds to the average of squared sampling
  weights. Default value is `NULL`.

- ...:

  Additional arguments.

## Value

A `cpreval` object with the structure:

- estimate: Estimated proportion.

- sd: Estimated standard error of the estimator.

- ci_asym: Asymptotic confidence interval at the 1 - gamma confidence
  level.

- gamma: Confidence level (i.e. 1 - gamma) for confidence intervals.

- method: Estimation method (in this case sample survey).

- measurement: A vector with (alpha0, alpha, beta).

- boundary: A boolean variable indicating if the estimates falls at the
  boundary of the parameter space.

- pi0: Value of pi0 (input value).

- sampling: Type of sampling considered ("random" or "weighted").

- V: Average sum of squared sampling weights if weighted/stratified is
  used (otherwise NULL).

- ...: Additional parameters.

## Author

Stephane Guerrier, Maria-Pia Victoria-Feser, Christoph Kuzmics

## Examples

``` r
# Samples without measurement error
X = sim_Rs(theta = 30/1000, pi0 = 10/1000, n = 1500, seed = 18)
survey_mle(R = X$R, n = X$n)
#> Method: Survey MLE
#> 
#> Estimated proportion: 3.4667%
#> Standard error      : 0.4723%
#> 
#> Confidence intervals at the 95% level:
#> Asymptotic Approach: 2.5409% - 4.3924%
#> Clopper-Pearson    : 2.5997% - 4.5214%
#> 
#> Assumed measurement error: alpha = 0%, beta = 0% 
#> Sampling: Random

# With measurement error
X = sim_Rs(theta = 30/1000, pi0 = 10/1000, n = 1500, alpha = 0.01, beta = 0.05, seed = 18)
survey_mle(R = X$R, n = X$n)
#> Method: Survey MLE
#> 
#> Estimated proportion: 3.7333%
#> Standard error      : 0.4895%
#> 
#> Confidence intervals at the 95% level:
#> Asymptotic Approach: 2.7740% - 4.6927%
#> Clopper-Pearson    : 2.8322% - 4.8209%
#> 
#> Assumed measurement error: alpha = 0%, beta = 0% 
#> Sampling: Random
survey_mle(R = X$R, n = X$n, alpha = 0.01, beta = 0.05)
#> Method: Survey MLE
#> 
#> Estimated proportion: 2.9078%
#> Standard error      : 0.5207%
#> 
#> Confidence intervals at the 95% level:
#> Asymptotic Approach: 1.8872% - 3.9284%
#> Clopper-Pearson    : 1.9492% - 4.0648%
#> 
#> Assumed measurement error: alpha = 1%, beta = 5% 
#> Sampling: Random
```
