
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![Licence](https://img.shields.io/badge/licence-AGPL--3.0-blue.svg)](https://opensource.org/licenses/AGPL-3.0)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-4.0.0-6666ff.svg)](https://cran.r-project.org/)
[![Last-changedate](https://img.shields.io/badge/last%20change-2023--10--05-green.svg)](https://github.com/stephaneguerrier/pempi)
[![R-CMD-check](https://github.com/stephaneguerrier/pempi/workflows/R-CMD-check/badge.svg)](https://github.com/stephaneguerrier/pempi/actions)
<!-- badges: end -->

# `pempi` Overview

The proportion estimation with marginal proxy information (`pempi`)
package, allows to estimate and build confidence intervals for
proportions, from random or stratified samples and census data with
participation bias. Measurement errors in the form of false positive and
false negative are also included in the inferential procedure. The
`pempi` package also contains code for simulation studies and
sensitivity analysis reported in the companion paper Guerrier et
al. (2023), as well as the Austrian dataset on COVID-19 prevalence in
November 2020.

# Remark on notation

The notation and conventions used in Guerrier et al. (2023) are slightly
amended for convenience in this package. In particular, we use R1
instead R11, R2 instead of R10, R3 instead of R01 and R4 instead of R00.

# Package installation

The `pempi` package can be installed from GitHub as follows:

``` r
# Install devtools
install.packages("devtools")

# Install the package from GitHub
devtools::install_github("stephaneguerrier/pempi")
```

Note that Windows users are assumed that have Rtools installed (if this
is not the case, please visit this
[link](https://cran.r-project.org/bin/windows/Rtools/)).

# How to cite

    @Manual{guerrier2023cape,
        title = {{pempi}: Proportion estimation with marginal proxy information},
        author = {Guerrier, S and Kuzmics, C and Victoria-Feser, M.-P.},
        year = {2023},
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

Guerrier, Stéphane, Christoph Kuzmics, and Maria-Pia Victoria-Feser.
2023. “Assessing COVID-19 Prevalence in Austria with Infection Surveys
and Case Count Data as Auxiliary Information”,
<https://arxiv.org/abs/2012.10745>.
