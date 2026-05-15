# Compute moment-based estimator.

Proportion estimated using the moment-based estimator and confidence
intervals based the asymptotic distribution of the estimator as well as
the Clopper-Pearson approach.

## Usage

``` r
moment_estimator(
  R3,
  n,
  pi0,
  gamma = 0.05,
  alpha = 0,
  beta = 0,
  alpha0 = 0,
  V = NULL,
  ...
)
```

## Arguments

- R3:

  A `numeric` that provides the number of participants in the survey
  sample that are tested positive only with the second testing device.

- n:

  A `numeric` that provides the sample size.

- pi0:

  A `numeric` that provides the prevalence or proportion of people (in
  the whole population) who are positive, as measured through a
  non-random, but systematic sampling (e.g. based on medical selection).

- gamma:

  A `numeric` that is used to compute a (1 - gamma) confidence region
  for the proportion. Default value is `0.05`.

- alpha:

  A `numeric` that provides the False Negative (FN) rate for the
  sample R. Default value is `0`.

- beta:

  A `numeric` that provides the False Positive (FP) rate for the
  sample R. Default value is `0`.

- alpha0:

  A `numeric` that corresponds to the probability that a random
  participant has been incorrectly declared positive through the
  nontransparent procedure. In most applications, this probability is
  likely very close to zero. Default value is `0`.

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

- ci_cp: Confidence interval (1 - gamma confidence level) based on the
  Clopper-Pearson approach.

- gamma: Confidence level (i.e. 1 - gamma) for confidence intervals.

- method: Estimation method (in this case moment estimator).

- measurement: A vector with (alpha0, alpha, beta).

- beta0: Estimated false negative rate of the official procedure.

- ci_beta0: Asymptotic confidence interval (1 - gamma confidence level)
  for beta0.

- boundary: A boolean variable indicating if the estimates falls at the
  boundary of the parameter space.

- pi0: Value of pi0 (input value).

- sampling: Type of sampling considered ("random" or "weighted").

- V: Average sum of squared sampling weights if weighted/stratified is
  used (otherwise NULL).

- n: Sample size.

- avar_beta0: Estimated asymptotic variance of beta0

- ...: Additional parameters.

## Author

Stephane Guerrier, Maria-Pia Victoria-Feser, Christoph Kuzmics

## Examples

``` r
# Samples without measurement error
X = sim_Rs(theta = 3/100, pi0 = 1/100, n = 1500, seed = 18)
moment_estimator(R3 = X$R3, n = X$n, pi0 = X$pi0)
#> Method: Moment Estimator
#> 
#> Estimated proportion: 3.2000%
#> Standard error      : 0.3787%
#> 
#> Confidence intervals at the 95% level:
#> Asymptotic Approach: 2.4577% - 3.9423%
#> Clopper-Pearson    : 2.5191% - 4.0758%
#> 
#> Assumed measurement error: alpha  = 0%, beta = 0%,
#>                            alpha0 = 0% 
#> 
#> Estimated false negative rate of the
#> official procedure: beta0 = 68.75%
#> CI at the 95% level: 61.50% - 76.00%
#> 
#> Estimated ascertainment rate: 
#> pi0/pi = 31.25%
#> CI at the 95% level: 24.00% - 38.50%
#> 
#> Sampling: Random

# With measurement error
X = sim_Rs(theta = 3/100, pi0 = 1/100, n = 1500, alpha0 = 0.001,
alpha = 0.01, beta = 0.05, seed = 18)
moment_estimator(R3 = X$R3, n = X$n, pi0 = X$pi0)
#> Method: Moment Estimator
#> 
#> Estimated proportion: 3.6667%
#> Standard error      : 0.4160%
#> 
#> Confidence intervals at the 95% level:
#> Asymptotic Approach: 2.8514% - 4.4820%
#> Clopper-Pearson    : 2.9118% - 4.6137%
#> 
#> Assumed measurement error: alpha  = 0%, beta = 0%,
#>                            alpha0 = 0% 
#> 
#> Estimated false negative rate of the
#> official procedure: beta0 = 72.73%
#> CI at the 95% level: 66.66% - 78.79%
#> 
#> Estimated ascertainment rate: 
#> pi0/pi = 27.27%
#> CI at the 95% level: 21.21% - 33.34%
#> 
#> Sampling: Random
moment_estimator(R3 = X$R3, n = X$n, pi0 = X$pi0, alpha0 = 0.001,
alpha = 0.01, beta = 0.05)
#> Method: Moment Estimator
#> 
#> Estimated proportion: 2.6864%
#> Standard error      : 0.4430%
#> 
#> Confidence intervals at the 95% level:
#> Asymptotic Approach: 1.8182% - 3.5546%
#> Clopper-Pearson    : 1.8825% - 3.6948%
#> 
#> Assumed measurement error: alpha  = 1%, beta = 5%,
#>                            alpha0 = 0.1% 
#> 
#> Estimated false negative rate of the
#> official procedure: beta0 = 66.40%
#> CI at the 95% level: 55.57% - 77.23%
#> 
#> Estimated ascertainment rate: 
#> pi0/pi = 33.60%
#> CI at the 95% level: 22.77% - 44.43%
#> 
#> Sampling: Random
```
