
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/stephaneguerrier/cape.svg?branch=master)](https://travis-ci.com/stephaneguerrier/cape)
[![Licence](https://img.shields.io/badge/licence-AGPL--3.0-blue.svg)](https://opensource.org/licenses/AGPL-3.0)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.6.0-6666ff.svg)](https://cran.r-project.org/)
[![Last-changedate](https://img.shields.io/badge/last%20change-2021--01--19-green.svg)](https://github.com/stephaneguerrier/cape)
[![CRAN
status](https://www.r-pkg.org/badges/version/stacks)](https://CRAN.R-project.org/package=stacks)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-blue.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

# `cape` Overview

The ContionAl Prevalence Estimation (cape) package, allows to estimate
and build confidence intervals for proportions, from random samples and
census data with participation bias. Measurement errors in the form of
false positive and false negative are also included the the inferential
procedure. The cape package also contains code for simulation studies
and sensitivity analysis reported in the companion paper Guerrier,
Kuzmics, and Victoria-Feser (2020).

# Remark on notation

The notation and conventions used in Guerrier, Kuzmics, and
Victoria-Feser (2020) are slightly amended for convenience in this
package. In particular, we use R1 instead R11, R2 instead of R10, R3
instead of R01 and R4 instead of R00.

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

    @Manual{guerrier2020cape,
        title = {{cape}: Conditional Prevalence Estimation using Random and Non-Random Sample Information},
        author = {Guerrier, S and Kuzmics, C and Victoria-Feser, M.-P.},
        year = {2020},
        note = {R package},
        url = {https://github.com/stephaneguerrier/cape}
    }

# Graphical user interface

A Graphical User Interface (GUI) is available with the `cape` R package.
This GUI which can be launched as follows:

``` r
gui()
```

It allows to compare the methods considered in the package as
illustrated below:

![](https://i.imgur.com/k0UrJMC.gif)

# License

The license this source code is released under is the GNU AFFERO GENERAL
PUBLIC LICENSE (AGPL) v3.0. Please see the LICENSE file for full text.
Otherwise, please consult [TLDR
Legal](https://tldrlegal.com/license/gnu-affero-general-public-license-v3-\(agpl-3.0\))
or [GNU](https://www.gnu.org/licenses/agpl-3.0.en.html) which will
provide a synopsis of the restrictions placed upon the code.

# References

<div id="refs" class="references">

<div id="ref-guerrier2020prevalence">

Guerrier, Stéphane, Christoph Kuzmics, and Maria-Pia Victoria-Feser.
2020. “Prevalence Estimation from Random Samples and Census Data with
Participation Bias.” <http://arxiv.org/abs/2012.10745>.

</div>

</div>
