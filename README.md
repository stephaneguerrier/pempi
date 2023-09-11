
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![Licence](https://img.shields.io/badge/licence-AGPL--3.0-blue.svg)](https://opensource.org/licenses/AGPL-3.0)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-4.0.0-6666ff.svg)](https://cran.r-project.org/)
[![Last-changedate](https://img.shields.io/badge/last%20change-2023--09--11-green.svg)](https://github.com/stephaneguerrier/cape)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-blue.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![R-CMD-check](https://github.com/stephaneguerrier/cape/workflows/R-CMD-check/badge.svg)](https://github.com/stephaneguerrier/cape/actions)
<!-- badges: end -->

# `cape` Overview

The ContionAl Prevalence Estimation (cape) package, allows to estimate
and build confidence intervals for proportions, from random or
stratified samples and census data with participation bias. Measurement
errors in the form of false positive and false negative are also
included in the inferential procedure. The cape package also contains
code for simulation studies and sensitivity analysis reported in the
companion paper Guerrier et al. (2022), as well as the Austrian dataset
on COVID-19 prevalence in November 2020.

# Remark on notation

The notation and conventions used in Guerrier et al. (2022) are slightly
amended for convenience in this package. In particular, we use R1
instead R11, R2 instead of R10, R3 instead of R01 and R4 instead of R00.

# Package installation

The `cape` package can be installed from GitHub as follows:

``` r
# Install devtools
install.packages("devtools")

# Install the package from GitHub
devtools::install_github("stephaneguerrier/cape")
```

Note that Windows users are assumed that have Rtools installed (if this
is not the case, please visit this
[link](https://cran.r-project.org/bin/windows/Rtools/)).

# How to cite

    @Manual{guerrier2022cape,
        title = {{cape}: Conditional Prevalence Estimation using Random and Non-Random Sample Information},
        author = {Guerrier, S and Kuzmics, C and Victoria-Feser, M.-P.},
        year = {2022},
        note = {R package},
        url = {https://github.com/stephaneguerrier/cape}
    }

# License

The license this source code is released under is the GNU AFFERO GENERAL
PUBLIC LICENSE (AGPL) v3.0. Please see the LICENSE file for full text.
Otherwise, please consult [TLDR
Legal](https://tldrlegal.com/license/gnu-affero-general-public-license-v3-(agpl-3.0))
or [GNU](https://www.gnu.org/licenses/agpl-3.0.en.html) which will
provide a synopsis of the restrictions placed upon the code.

# References

Guerrier, Stéphane, Christoph Kuzmics, and Maria-Pia Victoria-Feser.
2022. “Assessing Coronavirus SARS-CoV-2 Prevalence in Austria in 2020,
with Sample Surveys and Census Data with Participation Bias”,
<http://arxiv.org/abs/2012.10745>.
