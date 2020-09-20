
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Travis-CI Build
Status](https://travis-ci.org/SMAC-Group/simts.svg?branch=master)](https://travis-ci.org/SMAC-Group/simts)
[![Licence](https://img.shields.io/badge/licence-AGPL--3.0-blue.svg)](https://opensource.org/licenses/AGPL-3.0)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.4.0-6666ff.svg)](https://cran.r-project.org/)
[![Last-changedate](https://img.shields.io/badge/last%20change-2020--09--20-green.svg)](https://github.com/SMAC-Group/simts)

# `cape` Overview

This R package provides an implementation of the ContionAl Prevalence
Estimation (or `cape`) apporach proposed in *Accurate Prevalence
Estimation for Emerging or Rare Infectious Diseases* by Stéphane
Guerrier, Christoph Kuzmics and Maria-Pia Victoria-Feser (submitted
manuscript available upon request).

In epidemiological studies, an important quantity that needs to be
estimated with great accuracy is the prevalence rate of a disease in
order to learn more about central parameters such as case fatality rate
and/or to plan and refine decisions about measures with regards to an
epidemic or a pandemic. Traditionally, to measure the prevalence, a
survey sample (of randomly chosen subjects from the population) is
collected and the prevalence is estimated by the sample proportion. This
process involves some financial and logistic efforts, which increase
with the number of sampled participants, while at the same time,
increasing the accuracy of the estimator. Having a sufficiently large
sample is especially important when the (true) prevalence is very small,
for example at the beginning of a new pandemic such as the one of the
COVID-19. In this case, since the beginning of the outbreak, many
measurements have been taken on the number of infected people, but only
on the sub-population selected under some medical (and logistic)
criteria. It is obvious that using these measurements as if they would
be like a complete census, would lead to an underestimation of the
prevalence.

In Guerrier, Kuzmics, and Victoria-Feser (2020), we propose to
adequately use this information, together with data from a survey
sample, in order to improve the accuracy of the prevalence estimation.
In other terms, for a given (legal) statistical precision that might be
required by authorities that finance a survey, the sample size can be a
lot smaller if one adequately uses the information provided by the
information collected in the sub-population. The possible
misclassification errors of the (medical) testing devices used to
collect the data, are also taken into account. The misclassification
errors are actually induced by the sensitivity, i.e. the complement to
the false positive rate, and by the specificity, i.e. the complement to
the false negative rate, of the medical testing devices. The approach is
a frequentist one, i.e. using cutoff values for the sensitivity and
specificity, hence without the need to specify a (prior) distribution
for these quantities.

# Package installation

The `cape` pacakge can be installed from GitHub as follows:

``` r
# Install devtools
install.packages("devtools")

# Install the package from GitHub
devtools::install_github("stephaneguerrier/CPreval")
```

# Model

## Mathematical setup

Consider taking a (random) survey sample of \(n\) participants in some
population in order to estimate the prevalence \(\pi\) of, for example,
a given infectious disease. The framework also supposes that prior to
the collection of the survey sample, measurement have been taken on a
sub-population, not necessarily randomly chosen. Conditionally on this
information, one can get the prevalence computed as the number of
positive cases (in the sub-population) relative to the (whole)
population size, i.e. a prevalence parameter \(\pi_0\) which is such
that \(\pi_0 \leq \pi\). While \(\pi_0\) is known, it cannot be directly
used as a proxy for \(\pi\), but the information provided by \(\pi_0\),
can be used to improve the accuracy of the estimation of \(\pi\).

More precisely, for each participant in the survey sample, the following
information should be available:

Note that \(\Pr(X_i=1\vert V_i=1) = 1\) and
\(\Pr(X_i=0\vert V_i=1) = 0\). The objective is to provide an estimator
for the unknown prevalence

However, if we allow for possible misclassification error, the variables
in  are latent, and what is observed are the following variables

While \(Y_i\) is the minimal information collected in the survey sample,
\(W_i\) is an additional information that could be obtained by asking
the participants to provide it. It can be obtained while collecting the
information on \(Y_i\), or a posteriori using a followup procedure. In
that case, it is not necessary to proceed to the followup for all
sampled participants, but only for those with \(Y_i=1\).

We also suppose that the (medical) testing devices have known
sensitivity and specificity, or equivalently, False Positive (FP) rates
\(\alpha = 1-\mbox{specificity}\) and False Negative (FN) rates
\(\beta = 1-\mbox{sensitivity}\), that might be different for the two
populations. Therefore, we define the following quantities:

For proper modelling Guerrier, Kuzmics, and Victoria-Feser (2020), we
assume that \(\alpha+\beta<1\) and \(\alpha_0+\beta_0<1\). From the
above definitions, we can deduce

with \(\Delta_0=1-(\alpha_0+\beta_0)\). In the figure below the
probability tree associated to the variables presented previously:

<center>

<img src="./figures/prob.png" alt=" " width="600"/>

</center>

From these variables we deduce the following random variables that are
necessary for estimation and inference:

with:

  - \(R_1\): the number of participants in the survey sample that were
    tested positive with both (medical) testing devices (and are, thus,
    members of the sub-population);
  - \(R_3\) is the number of participants in the survey sample that are
    tested positive only with the second testing device.

If easily available, on can add the following information, that might
increase estimation and inference accuracy for survey sample sizes
smaller than \(n=1,500\):

so that \(R_0:= \sum_{i=1}^nW_i=R_1+R_2\), with

  - \(R_0\): the number of participants in the survey sample that were
    tested positive in the sub-population;
  - \(R_2\): the number of participants in the survey sample that are
    tested positive only with the first testing device (and are, thus,
    members of the sub-population).

The respective distributions of \(R_1,R_2,R_3\) are binomial
distributions with success probabilities given below in the next
section.

## Simulating data

Data can be simulated using the function `sim_Rs()` as follows:

``` r
# Measurement error
alpha = alpha0 = 0.01
beta = beta0 = 0.02

# Simulation "seed"
seed = 18

# True prevalence
theta = 4/100

# Prevalence in sub-population
pi0 = 2/100

# Sample size
n = 1500

# Simulation data
R = sim_Rs(theta = theta, pi0 = pi0, n = n, alpha0 = alpha0,
           alpha = alpha, beta0 = beta0, beta = beta, seed = seed)

# Print results
R
#> Data: R = 71, R0 = 36, n = 1500
#>       R1 = 19, R2 = 17, R3 = 52, R4 = 1412
#> 
#> Assumed measurement error: alpha0 = 1%, alpha = 1%, 
#>                            beta0  = 2%, beta  = 2%
```

# Estimators

The success probabilities associated to \(R_j,j=1,2,3\) are functions of
the (true) prevalence \(\pi\), and are given by

where \(\Delta :=1-(\alpha+\beta)\). A likelihood function can be built
based on the multinomial distribution (with categories
\(R_1,R_2,R_3,R_4\), with \(R_4=1-R_1-R_2-R_3\)) to obtain the Maximum
Likelihood Estimator (MLE) for \(\pi\) Guerrier, Kuzmics, and
Victoria-Feser (2020). The MLE has a closed form expressions if
\(\alpha_0=0\). If the information provided by \(R_2\) is not available,
the likelihood function is marginalized out on \(R_2\) and one can still
obtain an MLE for \(\pi\). Alternatively, the package also contains a
Method of Moments Estimator (MME) based on \(R_3\) (with expectation
\(\tau_3(\pi)\)), which is closed form for any
\(\alpha_0,\alpha,\beta_0,\beta\).

The confidence intervals for the MLE and the marginal MLE can be
computed using asymptotic theory. They provide good coverage, with small
interval length for survey sample sizes from \(n=1,000\) Guerrier,
Kuzmics, and Victoria-Feser (2020). For the MME, since its finite sample
distribution is known through \(R_3\sim\mathcal{B}(n,\tau_3(\pi))\),
exact confidence bounds can be computed using, for example, the
(fiducial) Clopper–Pearson (CP) approach method, see Clopper and Pearson
(1934), Fisher (1935) and Brown, Cai, and DasGupta (2001). In the
package are implemented, the CP method based on \(R_3\), and the
asymptotic theory method for the MLE (with \(R_2\) known or not). An
extensive simulation study is provided in Guerrier, Kuzmics, and
Victoria-Feser (2020) that compares coverage probabilities with the
different approaches.

To obtain these estimators, one can use the functions `mle()`,
`marginal_mle()` and `moment_estimator()`. Indeed, the MLE can be
estimate as follows on the data previously
simulated:

``` r
fit_mle = mle(R1 = R$R1, R2 = R$R2, R3 = R$R3, R4 = R$R4, n = R$n, pi0 = R$pi0,
              alpha0 = R$alpha0, alpha = R$alpha, 
              beta0 = R$beta0, beta = R$beta)
fit_mle
#> Method: MLE
#> 
#> Estimated proportion: 3.6168%
#> Standard error      : 0.4926%
#> 
#> Confidence intervals at the 95% level:
#> Asymptotic Approach: 2.6512% - 4.5823%
#> 
#> Assumed measurement error: alpha0 = 1%, alpha = 1%, 
#>                            beta0  = 2%, beta  = 2%
```

Note that the MLE is based on \(R_1\), \(R_2\), \(R_3\) and \(R_4\)
while the marginal MLE simply considered the information in \(R_1\) and
\(R_3\). This estimator can be obtained as follows:

``` r
fit_mmle = marginal_mle(R1 = R$R1, R3 = R$R3, n = R$n, pi0 = R$pi0,
                        alpha0 = R$alpha0, alpha = R$alpha, 
                        beta0 = R$beta0, beta = R$beta)
fit_mmle
#> Method: Marginal MLE
#> 
#> Estimated proportion: 3.6160%
#> Standard error      : 0.4927%
#> 
#> Confidence intervals at the 95% level:
#> Asymptotic Approach: 2.6504% - 4.5816%
#> 
#> Assumed measurement error: alpha0 = 1%, alpha = 1%, 
#>                            beta0  = 2%, beta  = 2%
```

Moreover, the moment estimator, which is simply based on \(R_3\), can be
computed as follows:

``` r
fit_mme = moment_estimator(R3 = R$R3, n = R$n, pi0 = R$pi0,
                           alpha0 = R$alpha0, alpha = R$alpha, 
                           beta0 = R$beta0, beta = R$beta)
fit_mme
#> Method: Moment Estimator
#> 
#> Estimated proportion: 3.5996%
#> Standard error      : 0.4919%
#> 
#> Confidence intervals at the 95% level:
#> Asymptotic Approach: 2.6355% - 4.5636%
#> Clopper-Pearson    : 2.6968% - 4.6979%
#> 
#> Assumed measurement error: alpha0 = 1%, alpha = 1%, 
#>                            beta0  = 2%, beta  = 2%
```

It can be observed that this estimator provides asymtoptic and CP
confidence intervals. Finally for comparison purpuses, we implemented
the standard sample average estimator which can be obtained as
follows:

``` r
fit_mean = survey_sample(R = R$R, n = R$n, alpha = R$alpha, beta = R$beta)
fit_mean
#> Method: Survey sample
#> 
#> Estimated proportion: 3.8488%
#> Standard error      : 0.5652%
#> 
#> Confidence intervals at the 95% level:
#> Asymptotic Approach: 2.7409% - 4.9567%
#> Clopper-Pearson    : 2.7990% - 5.0858%
#> 
#> Assumed measurement error: alpha = 1%, beta = 2%
```

# License

The license this source code is released under is the GNU AFFERO GENERAL
PUBLIC LICENSE (AGPL) v3.0. Please see the LICENSE file for full text.
Otherwise, please consult [TLDR
Legal](https://tldrlegal.com/license/gnu-affero-general-public-license-v3-\(agpl-3.0\))
or [GNU](https://www.gnu.org/licenses/agpl-3.0.en.html) which will
provide a synopsis of the restrictions placed upon the code.

# References

<div id="refs" class="references">

<div id="ref-brown2001">

Brown, L. D., T. Cai, and A. DasGupta. 2001. “Interval Estimation for a
Binomial Proportion.” *Statistical Science* 16: 101–33.

</div>

<div id="ref-ClPe34">

Clopper, C. J., and E. S. Pearson. 1934. “The Use of Confidence or
Fiducial Limits Illustrated in the Case of the Binomial.” *Biometrika*
26: 404–13.

</div>

<div id="ref-fisher35">

Fisher, R. A. 1935. “The Fiducial Argument in Statistical Inference.”
*Annals of Eugenics* 6: 391–98.

</div>

<div id="ref-guerrier2020accurate">

Guerrier, S, C Kuzmics, and M.-P. Victoria-Feser. 2020. “Accurate
Prevalence Estimation for Emerging or Rare Infectious Diseases.”

</div>

</div>
