# Package index

## Estimation

Prevalence estimators with measurement-error correction.

- [`conditional_mle()`](https://stephaneguerrier.github.io/pempi/reference/conditional_mle.md)
  : Compute MLE based on the full information R1, R2, R3 and R4.
- [`marginal_mle()`](https://stephaneguerrier.github.io/pempi/reference/marginal_mle.md)
  : Compute (marginalized) MLE based on the partial information R1 and
  R3.
- [`moment_estimator()`](https://stephaneguerrier.github.io/pempi/reference/moment_estimator.md)
  : Compute moment-based estimator.
- [`survey_mle()`](https://stephaneguerrier.github.io/pempi/reference/survey_mle.md)
  : Compute proportion in the survey sample (standard estimator)
- [`update_prevalence()`](https://stephaneguerrier.github.io/pempi/reference/update_prevalence.md)
  : Update prevalence using new case prevalence rates

## Data

Austrian COVID-19 prevalence survey (November 2020).

- [`covid19_austria`](https://stephaneguerrier.github.io/pempi/reference/covid19_austria.md)
  : COVID-19 Survey Data from Statistics Austria (November 2020)

## Simulation

Generate synthetic surveys for sensitivity studies.

- [`sim_Rs()`](https://stephaneguerrier.github.io/pempi/reference/sim_Rs.md)
  : Simulate data (R, R0, R1, R2, R3 and R4)
- [`get_prob()`](https://stephaneguerrier.github.io/pempi/reference/get_prob.md)
  : Compute sucess probabilities (tau_j's)

## Print Methods

S3 methods for `cpreval` and `cpreval_sim` objects.

- [`print(`*`<cpreval>`*`)`](https://stephaneguerrier.github.io/pempi/reference/print.cpreval.md)
  : Print estimation results
- [`print(`*`<cpreval_sim>`*`)`](https://stephaneguerrier.github.io/pempi/reference/print.cpreval_sim.md)
  : Print (simulated) sample
