# Simulate data (R, R0, R1, R2, R3 and R4)

Simulation function for random variables of interest.

## Usage

``` r
sim_Rs(theta, pi0, n, alpha0 = 0, alpha = 0, beta = 0, seed = NULL, ...)
```

## Arguments

- theta:

  A `numeric` that provides the true prevalence of a given disease.

- pi0:

  A `numeric` that provides the prevalence or proportion of people (in
  the whole population) who are positive, as measured through a
  non-random, but systematic sampling (e.g. based on medical selection).

- n:

  A `numeric` that corresponds to the sample size.

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

- seed:

  A `numeric` that provides the simulation seed. Default value is
  `NULL`.

- ...:

  Additional arguments.

## Value

A `cpreval_sim` object (`list`) with the structure:

- R: the number of participants in the survey sample that were tested
  positive.

- R0: the number of participants in the survey sample that were tested
  positive with the first testing device (and are, thus, members of the
  sub-population).

- R1: the number of participants in the survey sample that were tested
  positive with both (medical) testing devices (and are, thus, members
  of the sub-population).

- R2: the number of participants in the survey sample that are tested
  positive only with the first testing device (and are, thus, members of
  the sub-population).

- R3: the number of participants in the survey sample that are tested
  positive only with the second testing device.

- R4: the number of participants that are tested negative with the
  second testing device (and are either members of the sub-population
  and have tested negative with the first testing device or are not
  members of the sub-population).

- n: the sample size.

- alpha: the False Negative (FN) rate for the sample R.

- beta: the False Positive (FP) rate for the sample R.

- alpha0: the alpha0 probability (as defined above).

- ...: additional arguments.

## Author

Stephane Guerrier

## Examples

``` r
# Samples without measurement error
sim_Rs(theta = 3/100, pi0 = 1/100, n = 1500, seed = 18)
#> Data: R = 52, R0 = 19, n = 1500
#>       R1 = 19, R2 = 0, R3 = 33, R4 = 1448
#> 
#> Assumed measurement error: alpha = 0%, alpha0 = 0%, beta  = 0% 
#> 
#> False negative rate of the official procedure: beta0 = 66.67%

# With measurement error
sim_Rs(theta = 3/100, pi0 = 1/100, n = 1500, alpha0 = 0,
alpha = 0.01, beta = 0.05, seed = 18)
#> Data: R = 56, R0 = 19, n = 1500
#>       R1 = 18, R2 = 1, R3 = 38, R4 = 1443
#> 
#> Assumed measurement error: alpha = 1%, alpha0 = 0%, beta  = 5% 
#> 
#> False negative rate of the official procedure: beta0 = 66.67%
```
