# Negative marginalized log-likelihood function based on R0 and R

Log-Likelihood function based on R1 and R3 multiplied by -1.

## Usage

``` r
neg_log_lik_integrated(theta, Rvect, n, pi0, alpha0, alpha, beta, ...)
```

## Arguments

- theta:

  A `numeric` value for the paramerer of interest.

- Rvect:

  A `vector` of observations, i.e. (R1, R2, R3, R4), but R2 and R4 are
  NOT used.

- n:

  A `numeric` that provides the sample size.

- pi0:

  A `numeric` that provides the prevalence or proportion of people (in
  the whole population) who are positive, as measured through a
  non-random, but systematic sampling (e.g. based on medical selection).

- alpha0:

  A `numeric` that corresponds to the probability that a random
  participant has been incorrectly declared positive through the
  nontransparent procedure. In most applications, this probability is
  likely very close to zero.

- alpha:

  A `numeric` that provides the False Negative (FN) rate for the sample
  R.

- beta:

  A `numeric` that provides the False Positive (FP) rate for the sample
  R.

- ...:

  Additional arguments.

## Value

Negative marginalized log-likelihood.

## Author

Stephane Guerrier
