
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![Licence](https://img.shields.io/badge/licence-AGPL--3.0-blue.svg)](https://opensource.org/licenses/AGPL-3.0)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-4.0.0-6666ff.svg)](https://cran.r-project.org/)
[![Last-changedate](https://img.shields.io/badge/last%20change-2024--01--31-green.svg)](https://github.com/stephaneguerrier/pempi)
[![R-CMD-check](https://github.com/stephaneguerrier/pempi/workflows/R-CMD-check/badge.svg)](https://github.com/stephaneguerrier/pempi/actions)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/pempi)](https://www.r-pkg.org/pkg/pempi)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/grand-total/pempi)](https://www.r-pkg.org/pkg/pempi)
<!-- badges: end -->

# `pempi` Overview

The proportion estimation with marginal proxy information (`pempi`)
package, allows to estimate and build confidence intervals for
proportions, from random or stratified samples and census data with
participation bias. Measurement errors in the form of false positive and
false negative are also included in the inferential procedure. The
`pempi` package also contains code for simulation studies and
sensitivity analysis reported in the companion paper Guerrier et
al. (2024), as well as the Austrian dataset on COVID-19 prevalence in
November 2020.

# Remark on notation

The notation and conventions used in Guerrier et al. (2024) are slightly
amended for convenience in this package. In particular, we use `R1`
instead $\textit{R}_{11}$, `R2` instead of $\textit{R}_{10}$, `R3`
instead of $\textit{R}_{01}$ and `R4` instead of $\textit{R}_{00}$.

# Package installation

The `pempi` package is available on both CRAN and GitHub. The CRAN
version is considered stable, whereas the GitHub version is subject to
modifications/updates which may lead to installation problems or broken
functions. You can install the stable version of the `pempi` package
with:

``` r
install.packages("pempi")
```

The latest version can install from GitHub as follows:

``` r
# Install devtools
install.packages("devtools")

# Install the package from GitHub
devtools::install_github("stephaneguerrier/pempi")
```

Note that Windows users are assumed that have Rtools installed (if this
is not the case, please visit this
[link](https://cran.r-project.org/bin/windows/Rtools/)).

# Example: Austrian COVID-19 survey

In November 2020, a survey sample of $\textit{n}=2,287$ was collected in
Statistics Austria (2020) to test for COVID-19 using PCR-tests. In this
study, seventy-two participants were tested positive (i.e.,
`R1 + R3 = 72`)), and among these ones, thirty-five (`R1 = 35`) had
reported being tested positive with the official procedure, during the
same month. In November, there were $93,914$ declared cases among the
official (approximately) $7,166,167$ inhabitants in Austria (above 16
years old), so that $\pi_0 \approx 1.3105\%$. For simplicity, we
consider a random (unweighted) sampling and assume that the PCR false
positive and negative rate as well as the false case positive rate are
equal to 0. The data from this study can be obtained as follows:

``` r
# Load pempi
library(pempi)

# Austrian data (November 2020)
pi0 = 93914/7166167

# Load data
data("covid19_austria")

# Random sampling
n = nrow(covid19_austria)
R1 = sum(covid19_austria$Y == 1 & covid19_austria$Z == 1)
R2 = sum(covid19_austria$Y == 0 & covid19_austria$Z == 1)
R3 = sum(covid19_austria$Y == 1 & covid19_austria$Z == 0)
R4 = sum(covid19_austria$Y == 0 & covid19_austria$Z == 0)

# Print table
data_mat =c(R1, R2, R3, R4)
names(data_mat) = c("R1", "R2", "R3", "R4")
data_mat
#>   R1   R2   R3   R4 
#>   35    0   37 2218
```

The survey MLE as well as the conditional MLE and moment estimator
proposed in Guerrier et al. (2024) can be computed as follows:

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

# Reproducibility

All results, including figures, tables, real data analysis, and
simulations from Guerrier et al. (2024), can be reproduced as detailed
here:
<https://stephaneguerrier.github.io/pempi/articles/reproducibility.html>.

# How to cite

    @Manual{guerrier2024pempi,
        title = {{pempi}: Proportion estimation with marginal proxy information},
        author = {Guerrier, S and Kuzmics, C and Victoria-Feser, M.-P.},
        year = {2024},
        note = {R package},
        url = {https://github.com/stephaneguerrier/pempi}
    }

# License

The license this source code is released under is the GNU AFFERO GENERAL
PUBLIC LICENSE (AGPL) v3.0. Please see the LICENSE file for full text.
Otherwise, please consult
[GNU](https://www.gnu.org/licenses/agpl-3.0.en.html) which will provide
a synopsis of the restrictions placed upon the code.

# References

Guerrier, Stéphane, Christoph Kuzmics, and Maria-Pia Victoria-Feser,
“Assessing COVID-19 Prevalence in Austria with Infection Surveys and
Case Count Data as Auxiliary Information”, Journal of the American
Statistical Association, in press, 2024.

Statistics Austria, “Pravalenz von SARS-CoV-2-Infektionen liegt bei
3,1%”, 2020.
