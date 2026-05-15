# Compute sucess probabilities (tau_j's)

Compute joint probabilities of P(W = j, Y = k) for j, k = 0, 1.

## Usage

``` r
get_prob(theta, pi0, alpha, beta, alpha0)
```

## Arguments

- theta:

  A `numeric` that provides the true prevalence of a given disease.

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

## Value

A `vector` containing tau1, tau2, tau3 and tau4.

## Author

Stephane Guerrier

## Examples

``` r
prob1 = get_prob(theta = 0.02, pi0 = 0.01, alpha = 0, beta = 0, alpha0 = 0)
prob1
#> [1] 0.01 0.00 0.01 0.98
sum(prob1)
#> [1] 1

prob2 = get_prob(theta = 0.02, pi0 = 0.01, alpha = 0.001, beta = 0, alpha0 = 0.001)
prob2
#> [1] 0.00902098 0.00097902 0.01195902 0.97804098
sum(prob2)
#> [1] 1
```
