# Negative Weighted Log-Likelihood function

Log-Likelihood function based on R1, R2, R3 and R4 with sampling weights
multiplied by -1.

## Usage

``` r
neg_log_wlik(theta, Rvect, n, pi0, alpha, beta, alpha0, ...)
```

## Arguments

- theta:

  A `numeric` value for the paramerer of interest.

- Rvect:

  A `vector` of observations, i.e. (R1, R2, R3, R4).

- n:

  A `numeric` that provides the sample size.

- pi0:

  A `numeric` that provides the prevalence or proportion of people (in
  the whole population) who are positive, as measured through a
  non-random, but systematic sampling (e.g. based on medical selection).

- alpha:

  A `numeric` that provides the False Negative (FN) rate for the sample
  R.

- beta:

  A `numeric` that provides the False Positive (FP) rate for the sample
  R.

- alpha0:

  A `numeric` that corresponds to the probability that a random
  participant has been incorrectly declared positive through the
  nontransparent procedure. In most applications, this probability is
  likely very close to zero.

- ...:

  Additional arguments.

## Value

Negative log-likelihood.

## Author

Stephane Guerrier
