# Compute MLE based on the full information R1, R2, R3 and R4.

Proportion estimated using the MLE and confidence intervals based the
asymptotic distribution of the estimator.

## Usage

``` r
conditional_mle(
  R1 = NULL,
  R2 = NULL,
  R3 = NULL,
  R4 = NULL,
  n = R1 + R2 + R3 + R4,
  pi0,
  gamma = 0.05,
  alpha0 = 0,
  alpha = 0,
  beta = 0,
  V = NULL,
  ...
)
```

## Arguments

- R1:

  A `numeric` that provides the number of participants in the survey
  sample that were tested positive with both (medical) testing devices
  (and are, thus, members of the sub-population).

- R2:

  A `numeric` that provides the number of participants in the survey
  sample that are tested positive only with the first testing device
  (and are, thus, members of the sub-population).

- R3:

  A `numeric` that provides the number of participants in the survey
  sample that are tested positive only with the second testing device.

- R4:

  A `numeric` that provides the number of participants that are tested
  negative with the second testing device (and are either members of the
  sub-population and have tested negative with the first testing device
  or are not members of the sub-population).

- n:

  A `numeric` that provides the sample size. Default value R1 + R2 +
  R3 + R4. If this value is provided it is used to verify that R1 + R2 +
  R3 + R4 = n.

- pi0:

  A `numeric` that provides the prevalence or proportion of people (in
  the whole population) who are positive, as measured through a
  non-random, but systematic sampling (e.g. based on medical selection).

- gamma:

  A `numeric` that is used to compute a (1 - gamma) confidence region
  for the proportion. Default value is `0.05`.

- alpha0:

  A `numeric` that corresponds to the probability that a random
  participant has been incorrectly declared positive through the
  nontransparent procedure. In most applications, this probability is
  likely very close to zero. Default value is `0`.

- alpha:

  A `numeric` that provides the False Negative (FN) rate for the
  sample R. Default value is `0`.

- beta:

  A `numeric` that provides the False Positive (FP) rate for the
  sample R. Default value is `0`.

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

- method: Estimation method (in this case mle).

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
conditional_mle(R1 = X$R1, R2 = X$R2, R3 = X$R3, R4 = X$R4, pi0 = X$pi0)
#> Method: Conditional MLE
#> 
#> Estimated proportion: 3.2061%
#> Standard error      : 0.3792%
#> 
#> Confidence interval at the 95% level:
#> Asymptotic Approach: 2.4629% - 3.9494%
#> 
#> Assumed measurement error: alpha  = 0%, beta = 0%,
#>                            alpha0 = 0% 
#> 
#> Estimated false negative rate of the
#> official procedure: beta0 = 68.81%
#> CI at the 95% level: 61.58% - 76.04%
#> 
#> Estimated ascertainment rate: 
#> pi0/pi = 31.19%
#> CI at the 95% level: 23.96% - 38.42%
#> 
#> Sampling: Random

# With measurement error
X = sim_Rs(theta = 30/1000, pi0 = 10/1000, n = 1500, alpha0 = 0.001,
alpha = 0.01, beta0 = 0.05, beta = 0.05, seed = 18)
conditional_mle(R1 = X$R1, R2 = X$R2, R3 = X$R3, R4 = X$R4, pi0 = X$pi0)
#> Method: Conditional MLE
#> 
#> Estimated proportion: 3.6736%
#> Standard error      : 0.4164%
#> 
#> Confidence interval at the 95% level:
#> Asymptotic Approach: 2.8574% - 4.4899%
#> 
#> Assumed measurement error: alpha  = 0%, beta = 0%,
#>                            alpha0 = 0% 
#> 
#> Estimated false negative rate of the
#> official procedure: beta0 = 72.78%
#> CI at the 95% level: 66.73% - 78.83%
#> 
#> Estimated ascertainment rate: 
#> pi0/pi = 27.22%
#> CI at the 95% level: 21.17% - 33.27%
#> 
#> Sampling: Random
conditional_mle(R1 = X$R1, R2 = X$R2, R3 = X$R3, R4 = X$R4, pi0 = X$pi0,
alpha0 = 0.001, alpha = 0.01, beta = 0.05)
#> Method: Conditional MLE
#> 
#> Estimated proportion: 2.6929%
#> Standard error      : 0.4434%
#> 
#> Confidence interval at the 95% level:
#> Asymptotic Approach: 1.8238% - 3.5619%
#> 
#> Assumed measurement error: alpha  = 1%, beta = 5%,
#>                            alpha0 = 0.1% 
#> 
#> Estimated false negative rate of the
#> official procedure: beta0 = 66.48%
#> CI at the 95% level: 55.69% - 77.26%
#> 
#> Estimated ascertainment rate: 
#> pi0/pi = 33.52%
#> CI at the 95% level: 22.74% - 44.31%
#> 
#> Sampling: Random
```
