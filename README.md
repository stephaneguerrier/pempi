
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pempi — Proportion Estimation with Marginal Proxy Information

> Combine a small survey sample with population-level case-count data to
> estimate disease prevalence — with full handling of measurement error.
> The R companion to Guerrier, Kuzmics & Victoria-Feser (2024, *JASA*),
> winner of the **2024 JASA Reproducibility Award**.

<!-- badges: start -->

[![JASA Reproducibility
Award](https://img.shields.io/badge/JASA-2024%20Reproducibility%20Award-d97706?logo=academia&logoColor=white)](https://jasa-acs.github.io/repro-award/awardees.html)
[![DOI](https://img.shields.io/badge/DOI-10.1080%2F01621459.2024.2313790-4f46e5)](https://doi.org/10.1080/01621459.2024.2313790)
[![Licence](https://img.shields.io/badge/licence-AGPL--3.0-blue.svg)](https://opensource.org/licenses/AGPL-3.0)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-4.0.0-6666ff.svg)](https://cran.r-project.org/)
[![Last-changedate](https://img.shields.io/badge/last%20change-2026--05--15-green.svg)](https://github.com/stephaneguerrier/pempi)
[![R-CMD-check](https://github.com/stephaneguerrier/pempi/workflows/R-CMD-check/badge.svg)](https://github.com/stephaneguerrier/pempi/actions)
[![CRAN
downloads](http://cranlogs.r-pkg.org/badges/pempi)](https://www.r-pkg.org/pkg/pempi)
[![CRAN
total](https://cranlogs.r-pkg.org/badges/grand-total/pempi)](https://www.r-pkg.org/pkg/pempi)
<!-- badges: end -->

> 🏆 **Winner — 2024 JASA Reproducibility Award.** Awarded by the
> *Journal of the American Statistical Association* to Guerrier, Kuzmics
> & Victoria-Feser for the companion paper *Assessing COVID-19
> Prevalence in Austria with Infection Surveys and Case Count Data as
> Auxiliary Information* — JASA **119**(547):1722–1735. Full
> reproducibility material lives in the [Reproducibility
> vignette](https://stephaneguerrier.github.io/pempi/articles/reproducibility.html);
> see the [JASA awardees page](https://jasa-acs.github.io/repro-award/awardees.html)
> for the citation.

## Install

``` r
# CRAN (stable)
install.packages("pempi")

# GitHub (development)
# install.packages("devtools")
devtools::install_github("stephaneguerrier/pempi")
```

Windows users will need
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) installed.

## A Quick Taste — Austrian COVID-19 survey

In November 2020, Statistics Austria collected a survey sample of *n* =
2,287 to test for COVID-19 by PCR. 72 participants tested positive in
the survey (`R1 + R3 = 72`); 35 of them (`R1 = 35`) had also been
recorded positive by the official procedure that month. With 93,914
declared cases among ~7,166,167 inhabitants over 16, the official rate
is π₀ ≈ 1.31%.

``` r
# Load pempi
library(pempi)

# Austrian official rate (November 2020)
pi0 = 93914/7166167

# Bundled real-world dataset
data("covid19_austria")

# Random sampling
n  = nrow(covid19_austria)
R1 = sum(covid19_austria$Y == 1 & covid19_austria$Z == 1)
R2 = sum(covid19_austria$Y == 0 & covid19_austria$Z == 1)
R3 = sum(covid19_austria$Y == 1 & covid19_austria$Z == 0)
R4 = sum(covid19_austria$Y == 0 & covid19_austria$Z == 0)

# Print table
data_mat = c(R1, R2, R3, R4)
names(data_mat) = c("R1", "R2", "R3", "R4")
data_mat
#>   R1   R2   R3   R4 
#>   35    0   37 2218
```

The survey MLE alongside the conditional MLE and moment estimator from
Guerrier et al. (2024):

``` r
survey_mle(R = R1 + R3, n = n)
#> Method: Survey MLE
#> 
#> Estimated proportion: 3.1441%
#> Standard error      : 0.3647%
#> 
#> Confidence intervals at the 95% level:
#> Asymptotic Approach: 2.4294% - 3.8588%
#> Clopper-Pearson    : 2.4680% - 3.9433%
#> 
#> Assumed measurement error: alpha = 0%, beta = 0% 
#> Sampling: Random

conditional_mle(R1 = R1, R2 = R2, R3 = R3, R4 = R4, pi0 = pi0)
#> Method: Conditional MLE
#> 
#> Estimated proportion: 2.9317%
#> Standard error      : 0.2639%
#> 
#> Confidence interval at the 95% level:
#> Asymptotic Approach: 2.4145% - 3.4489%
#> 
#> Assumed measurement error: alpha  = 0%, beta = 0%,
#>                            alpha0 = 0% 
#> 
#> Estimated false negative rate of the
#> official procedure: beta0 = 55.30%
#> CI at the 95% level: 47.41% - 63.18%
#> 
#> Estimated ascertainment rate: 
#> pi0/pi = 44.70%
#> CI at the 95% level: 36.82% - 52.59%
#> 
#> Sampling: Random

moment_estimator(R3 = R3, n = n, pi0 = pi0)
#> Method: Moment Estimator
#> 
#> Estimated proportion: 2.9262%
#> Standard error      : 0.2635%
#> 
#> Confidence intervals at the 95% level:
#> Asymptotic Approach: 2.4099% - 3.4426%
#> Clopper-Pearson    : 2.4506% - 3.5308%
#> 
#> Assumed measurement error: alpha  = 0%, beta = 0%,
#>                            alpha0 = 0% 
#> 
#> Estimated false negative rate of the
#> official procedure: beta0 = 55.21%
#> CI at the 95% level: 47.31% - 63.12%
#> 
#> Estimated ascertainment rate: 
#> pi0/pi = 44.79%
#> CI at the 95% level: 36.88% - 52.69%
#> 
#> Sampling: Random
```

The conditional MLE and moment estimators correct the survey-only
estimate downward (from 3.14% to ~2.93%) by leveraging the auxiliary
case-count information — and shrink the standard error by roughly 30%.

## Where to Next

|     |                                                                                                                      |                                                                                                              |
|-----|----------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------|
| 📐  | [**Get Started — Methodology**](https://stephaneguerrier.github.io/pempi/articles/get_started.html)                  | The mathematical setup, the four R-counts, and how each estimator corrects the survey-only baseline.         |
| 🔬  | [**Reproducibility — JASA 2024**](https://stephaneguerrier.github.io/pempi/articles/reproducibility.html)            | Reproduce every table, figure, and simulation from the award-winning paper.                                  |
| 🇦🇹 | [**Austrian COVID-19 Survey**](https://stephaneguerrier.github.io/pempi/reference/covid19_austria.html)               | The bundled real-world dataset (*n* = 2,287, November 2020) used throughout the package.                     |
| 📑  | [**How to Cite**](https://stephaneguerrier.github.io/pempi/articles/citation.html)                                   | BibTeX, JASA citation, and DOI for the package and the underlying paper.                                     |

## Notation note

The notation in Guerrier, Kuzmics & Victoria-Feser (2024) is slightly
amended for convenience here. The package uses `R1`, `R2`, `R3`, `R4`
for the paper's $R_{11}$, $R_{10}$, $R_{01}$, $R_{00}$ respectively.

## Citation

If you use this package, please cite both the paper and the package.

**Paper**

    @article{guerrier2024prevalence,
      title   = {Assessing {COVID-19} Prevalence in {Austria} with Infection
                 Surveys and Case Count Data as Auxiliary Information},
      author  = {Guerrier, St\'ephane and Kuzmics, Christoph and
                 Victoria-Feser, Maria-Pia},
      journal = {Journal of the American Statistical Association},
      volume  = {119},
      number  = {547},
      pages   = {1722--1735},
      year    = {2024},
      doi     = {10.1080/01621459.2024.2313790}
    }

**Package**

    @Manual{guerrier2024pempi,
      title  = {{pempi}: Proportion Estimation with Marginal Proxy Information},
      author = {Guerrier, St\'ephane and Kuzmics, Christoph and
                Victoria-Feser, Maria-Pia},
      year   = {2024},
      note   = {R package},
      url    = {https://github.com/stephaneguerrier/pempi}
    }

## License

This source code is released under the GNU Affero General Public
License (AGPL) v3.0. See the [LICENSE](LICENSE) file or the
[GNU summary](https://www.gnu.org/licenses/agpl-3.0.en.html).

## References

Guerrier, S., Kuzmics, C., Victoria-Feser, M.-P. (2024). *Assessing
COVID-19 Prevalence in Austria with Infection Surveys and Case Count
Data as Auxiliary Information*. **Journal of the American Statistical
Association**, 119(547), 1722–1735.
[doi:10.1080/01621459.2024.2313790](https://doi.org/10.1080/01621459.2024.2313790)

Statistics Austria. (2020). *Prävalenz von SARS-CoV-2-Infektionen liegt
bei 3,1%*.
