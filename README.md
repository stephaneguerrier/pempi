
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

Harrar and Spence (2013)

This R package provides an implementation of the ContionAl Prevalence
Estimation (or `cape`) approach proposed in *Prevalence Estimation using
Random and Non-Random Sample Information* by Stéphane Guerrier,
Christoph Kuzmics and Maria-Pia Victoria-Feser (submitted manuscript
available upon request).

# Package installation

The `cape` pacakge can be installed from GitHub as follows:

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
This GUI which can be lauched as follows:

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

<div id="refs" class="references">

<div id="ref-harrar2013taste">

Harrar, Vanessa, and Charles Spence. 2013. “The Taste of Cutlery: How
the Taste of Food Is Affected by the Weight, Size, Shape, and Colour of
the Cutlery Used to Eat It.” *Flavour* 2 (1): 21.

</div>

</div>
