# \`pempi\` Vignette

## Introduction

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

In Guerrier et al. (2024), we propose to combine this case-count
information with data from a survey sample to improve the accuracy of
the prevalence estimate. For a target statistical precision, the
required sample size shrinks substantially when case-count data are used
as auxiliary information. The possible misclassification errors of the
medical tests are taken into account through the sensitivity
($`1-\beta`$, the complement of the false negative rate) and the
specificity ($`1-\alpha`$, the complement of the false positive rate).
The approach is frequentist — it uses fixed values for the sensitivity
and specificity rather than priors.

## Mathematical setup

As in Guerrier et al. (2024), we define the following unobserved random
variable:

``` math
\begin{equation*}
         X_i=
         \left\{
    \begin{array}{ll}
        1  & \quad \mbox{if participant (in the survey sample) $i$ is positive,} \\
        0  & \quad \mbox{otherwise;}
    \end{array}
\right.
\end{equation*}
```

The objective is to provide an estimator for the unknown population
proportion, i.e. the prevalence, given by

``` math
\begin{equation*}
        \pi = \Pr\left(X_i=1\right). 
\end{equation*}
```

Next, we define the following quantities.

``` math
\begin{eqnarray}
         Y_i&=&
         \left\{
    \begin{array}{ll}
        1  & \quad \mbox{if participant $i$ is tested positive in the survey sample,} \\
        0  & \quad \mbox{otherwise;}
    \end{array}
\right.\nonumber \\
         Z_i&=& \left\{
    \begin{array}{ll}
        1  & \quad \mbox{if participant $i$ was declared positive with the official procedure,} \\
        0  & \quad \mbox{otherwise.}
    \end{array}
\right. \label{eqn:Y-Z}
\end{eqnarray}
```

and

``` math
\begin{alignat}{4}
    &R_{11}&&=\sum_{i=1}^n Y_i Z_i, \quad\quad\quad 
    &&R_{10}&&=\sum_{i=1}^n(1-Y_i)Z_i, \label{eqn:R-R4} \\
    &R_{01}&&= \sum_{i=1}^nY_i(1-Z_i), \quad\quad \quad 
    &&R_{00}&&= \sum_{i=1}^n(1-Y_i)(1-Z_i)=n-R_{11}-R_{10}+R_{01}.
    \nonumber
\end{alignat}
```

In words, $`R_{11}`$ is the number of participants in the survey sample
that are tested positive and have also been declared positive through
the official procedure; $`R_{10}`$ is the number of participants in the
survey sample that are tested negative but have been declared positive
through the official procedure; $`R_{01}`$ is the number of participants
in the survey sample that are tested positive but have been declared
negative through the official procedure; $`R_{00}`$ is the number of
participants in the survey sample that are tested negative and have been
declared negative through the official procedure. We also make use of
$`R_{\ast 1} = \sum_{i=1}^nY_i = R_{11} + R_{01}`$, the number of
participants that are tested positive in the survey sample. The
$`R_{jk}`$ are central quantities for computing the proposed estimators.
In the case of stratified sampling, these quantities can be modified
accordingly (see below in the section on stratified sampling).

## Random sampling

For simplicity of exposition, we first consider the unweighted (or
random) sampling case, the weighted one is treated below.

We also allow for the possibility that the outcome of (medical) tests
can be subject to misclassification error. Hence, we define
$`\alpha = \Pr(Y_i=1\vert X_i=0)`$ and
$`\beta = \Pr(Y_i=0\vert X_i=1)`$. The probabilities $`\alpha`$ and
$`\beta`$, are the (assumed known) FP rate
($`\alpha = 1-\mbox{specificity}`$) and FN rate
($`\beta = 1-\mbox{sensitivity}`$) of the particular medical test
employed in the survey. Moreover, we will make use of the (known)
prevalence $`\pi_0 =\Pr(Z_i=1)`$ from the official procedure. As
previously explained, $`\pi_0`$ is the joint probability of being
selected in the official procedure and declared positive, so that we
have $`\pi_0\leq \pi`$.

Since the objective is to take advantage of the information provided by
the official procedure in estimating the prevalence $`\pi`$, we also
need to take into account the possible biases in the official data.
Therefore, we define $`\alpha_0 = \Pr(Z_i=1\vert X_i=0)`$, the FP rate
of the official procedure and $`\beta_0 = \Pr(Z_i=0\vert X_i=1)`$ the FN
rate of the official procedure. It turns out (see Guerrier et al. (2024)
for details) that $`\alpha_0`$ is a negligible quantity and hence can be
set to $`0`$, and $`\beta_0`$ can be deduced from the other available or
estimable quantities as we have

``` math
\beta_0 = 1- \frac{\pi_0-\alpha_0(1-\pi)}{\pi}.
```

The success probabilities (see Guerrier et al. (2024) for their
derivation), denoted by $`\tau_{jk}(\pi)`$ associated to each
$`R_{jk}`$, $`j,k \in \{0,1\}`$ are given by %
``` math
\begin{eqnarray}
    \begin{aligned}
    \tau_{11}(\pi) &=  \Pr(Z_i=1, Y_i=1) =\pi\Delta\alpha_0+(\pi_0-\alpha_0)(1-\beta)+\alpha\alpha_0,   \\
    \tau_{10}(\pi) &=  \Pr(Z_i=1,Y_i=0) = -\pi\Delta\alpha_0+(\pi_0-\alpha_0)\beta+(1-\alpha)\alpha_0,  \\
    \tau_{01}(\pi) &=  \Pr(Z_i=0,Y_i=1) =  \pi\Delta(1-\alpha_0)-(\pi_0-\alpha_0)(1-\beta)+\alpha(1-\alpha_0), \\
    \tau_{00}(\pi) &=  \Pr(Z_i=0,Y_i=0) = -\pi\Delta(1-\alpha_0)-(\pi_0-\alpha_0)\beta+(1-\alpha)(1-\alpha_0),
    \end{aligned}
    \label{eqn:tau-pi}
\end{eqnarray}
```
% where $`\Delta =1-(\alpha+\beta)`$. Without misclassification error,
we would have $`\tau_{11}(\pi)=\pi_0`$, $`\tau_{10}(\pi)=0`$,
$`\tau_{01}(\pi)=\pi-\pi_0`$, $`\tau_{00}(\pi)=1-\pi`$.

Based on these definitions, we present a Conditional Maximum Likelihood
Estimator (CMLE) and a Method of Moments Estimator (MME). The formal
derivations and properties are provided in Guerrier et al. (2024).
There, we also consider a Marginal MLE (MMLE) when some data is missing,
and some Generalized Method of Moment (GMM) estimators. The likelihood
function for $`\pi`$ can be obtained from the multinomial distribution
with categories provided by $`R_{11},R_{10},R_{01},R_{00}`$ and their
associated success probabilities
$`\tau_{11}(\pi),\tau_{10}(\pi),\tau_{01}(\pi),\tau_{00}(\pi)`$. The
CMLE, conditional on the information provided by the official procedure,
$`\widehat{\pi}`$, generally, has no closed-form solution but can be
computed numerically. However, in the case when $`\alpha_0 = 0`$, we
obtain a closed-form solution given by

``` math
\begin{equation}
    \widehat{\pi} = \frac{\pi_0 R_{00} + R_{01}}{\Delta \left(R_{01}+R_{00}\right)} - \frac{\pi_0 \beta}{\Delta} - \frac{\alpha}{\Delta}.
    \label{eq:cmle-closedform}
\end{equation}
```

When $`\alpha_0=\alpha=\beta=0`$, this further reduces to

``` math
\begin{equation}
   \widehat{\pi} = \pi_0 \frac{n - R_{\ast 1}}{n - R_{11}} + \frac{R_{01}}{ \left(n - R_{11}\right)}.
   \label{eqn:MLE}
\end{equation}
```

Considering the data used in Guerrier et al. (2024), this estimator can
be computed as follows:

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

# Compute CMLE
conditional_mle(R1 = R1, R2 = R2, R3 = R3, R4 = R4, pi0 = pi0)
```

    ## Method: Conditional MLE
    ## 
    ## Estimated proportion: 2.9317%
    ## Standard error      : 0.2639%
    ## 
    ## Confidence interval at the 95% level:
    ## Asymptotic Approach: 2.4145% - 3.4489%
    ## 
    ## Assumed measurement error: alpha  = 0%, beta = 0%,
    ##                            alpha0 = 0% 
    ## 
    ## Estimated false negative rate of the
    ## official procedure: beta0 = 55.30%
    ## CI at the 95% level: 47.41% - 63.18%
    ## 
    ## Estimated ascertainment rate: 
    ## pi0/pi = 44.70%
    ## CI at the 95% level: 36.82% - 52.59%
    ## 
    ## Sampling: Random

Note that the notation in the paper is slightly amended for convenience
in this package. In particular, the package uses `R1` for $`R_{11}`$,
`R2` for $`R_{10}`$, `R3` for $`R_{01}`$ and `R4` for $`R_{00}`$.
Considering the following measurement error $`\alpha = 0.01`$,
$`\alpha_0 = 0`$ and $`\beta = 0.1`$, we obtain:

``` r

# Assumed measurement errors
alpha0 = 0
alpha  = 1/100
beta   = 10/100

# Compute CMLE with measurement error
conditional_mle(R1 = R1, R2 = R2, R3 = R3, R4 = R4, pi0 = pi0,
                alpha = alpha, alpha0 = alpha0, beta = beta)
```

    ## Method: Conditional MLE
    ## 
    ## Estimated proportion: 2.0200%
    ## Standard error      : 0.2962%
    ## 
    ## Confidence interval at the 95% level:
    ## Asymptotic Approach: 1.4394% - 2.6006%
    ## 
    ## Assumed measurement error: alpha  = 1%, beta = 10%,
    ##                            alpha0 = 0% 
    ## 
    ## Estimated false negative rate of the
    ## official procedure: beta0 = 35.12%
    ## CI at the 95% level: 16.48% - 53.77%
    ## 
    ## Estimated ascertainment rate: 
    ## pi0/pi = 64.88%
    ## CI at the 95% level: 46.23% - 83.52%
    ## 
    ## Sampling: Random

Alternatively, we can consider an estimator from the class of GMM
estimators based on the random variable
$`\mathbf{R} =[R_{11}/n, R_{10}/n, R_{01}/n, R_{00}/n]`$ with
expectation
$`\mathbb{E}[\mathbf{R}] = {\tau}(\pi)=[\tau_{11}(\pi), \tau_{10}(\pi), \tau_{01}(\pi), \tau_{00}(\pi)]`$.
A particular case is given by an MME based on $`R_{01}`$ (with
expectation $`\tau_{01}(\pi)`$), which, again assuming an interior
solution exists, is given by

``` math
\begin{equation}
    \widetilde{\pi} = \frac{1}{\Delta(1-\alpha_0)} \left(\frac{R_{01}}{n} + \pi_0 - \beta\pi_0 - \alpha_0 \Delta - \alpha\right).
    \label{eqn:MME-ME}
\end{equation}
```

When $`\alpha_0=\alpha=\beta=0`$, this reduces to

``` math
\begin{equation}
   \widetilde{\pi} = \pi_0 + \frac{R_{01}}{n}.
    \label{eqn:MME}
\end{equation}
```

Since we have that $`\mathbb{E}[\widetilde{\pi} ] = \pi`$, the MME is
unbiased. This estimator can be computed as follows:

``` r

# Without measurement error
moment_estimator(R3 = R3, n = n, pi0 = pi0)
```

    ## Method: Moment Estimator
    ## 
    ## Estimated proportion: 2.9262%
    ## Standard error      : 0.2635%
    ## 
    ## Confidence intervals at the 95% level:
    ## Asymptotic Approach: 2.4099% - 3.4426%
    ## Clopper-Pearson    : 2.4506% - 3.5308%
    ## 
    ## Assumed measurement error: alpha  = 0%, beta = 0%,
    ##                            alpha0 = 0% 
    ## 
    ## Estimated false negative rate of the
    ## official procedure: beta0 = 55.21%
    ## CI at the 95% level: 47.31% - 63.12%
    ## 
    ## Estimated ascertainment rate: 
    ## pi0/pi = 44.79%
    ## CI at the 95% level: 36.88% - 52.69%
    ## 
    ## Sampling: Random

``` r

# With measurement error
moment_estimator(R3 = R3, n = n, pi0 = pi0, alpha = alpha,
                 alpha0 = alpha0, beta = beta)
```

    ## Method: Moment Estimator
    ## 
    ## Estimated proportion: 2.0171%
    ## Standard error      : 0.2960%
    ## 
    ## Confidence intervals at the 95% level:
    ## Asymptotic Approach: 1.4369% - 2.5973%
    ## Clopper-Pearson    : 1.4827% - 2.6963%
    ## 
    ## Assumed measurement error: alpha  = 1%, beta = 10%,
    ##                            alpha0 = 0% 
    ## 
    ## Estimated false negative rate of the
    ## official procedure: beta0 = 35.03%
    ## CI at the 95% level: 16.34% - 53.72%
    ## 
    ## Estimated ascertainment rate: 
    ## pi0/pi = 64.97%
    ## CI at the 95% level: 46.28% - 83.66%
    ## 
    ## Sampling: Random

In the `pempi` package the marginal MLE is also implemented and can be
used as follows:

``` r

# Without measurement error
marginal_mle(R1 = R1, R3 = R3, n = n, pi0 = pi0)
```

    ## Method: Marginal MLE
    ## 
    ## Estimated proportion: 2.9317%
    ## Standard error      : 0.2639%
    ## 
    ## Confidence interval at the 95% level:
    ## Asymptotic Approach: 2.4145% - 3.4489%
    ## 
    ## Assumed measurement error: alpha  = 0%, beta = 0%,
    ##                            alpha0 = 0% 
    ## 
    ## Estimated false negative rate of the
    ## official procedure: beta0 = 55.30%
    ## CI at the 95% level: 47.41% - 63.18%
    ## 
    ## Estimated ascertainment rate: 
    ## pi0/pi = 44.70%
    ## CI at the 95% level: 36.82% - 52.59%
    ## 
    ## Sampling: Random

``` r

# With measurement error
marginal_mle(R1 = R1, R3 = R3, n = n, pi0 = pi0, 
             alpha = alpha, beta = beta, alpha0 = alpha0)
```

    ## Method: Marginal MLE
    ## 
    ## Estimated proportion: 2.0200%
    ## Standard error      : 0.2962%
    ## 
    ## Confidence interval at the 95% level:
    ## Asymptotic Approach: 1.4394% - 2.6006%
    ## 
    ## Assumed measurement error: alpha  = 1%, beta = 10%,
    ##                            alpha0 = 0% 
    ## 
    ## Estimated false negative rate of the
    ## official procedure: beta0 = 35.12%
    ## CI at the 95% level: 16.48% - 53.77%
    ## 
    ## Estimated ascertainment rate: 
    ## pi0/pi = 64.88%
    ## CI at the 95% level: 46.23% - 83.52%
    ## 
    ## Sampling: Random

These results can be compared to the standard survey MLE which can be
computed as follows:

``` r

# Without measurement error
survey_mle(R = R1 + R3, n = n)
```

    ## Method: Survey MLE
    ## 
    ## Estimated proportion: 3.1441%
    ## Standard error      : 0.3647%
    ## 
    ## Confidence intervals at the 95% level:
    ## Asymptotic Approach: 2.4294% - 3.8588%
    ## Clopper-Pearson    : 2.4680% - 3.9433%
    ## 
    ## Assumed measurement error: alpha = 0%, beta = 0% 
    ## Sampling: Random

``` r

# With measurement error
survey_mle(R = R1 + R3, n = n, alpha = alpha, beta = beta)
```

    ## Method: Survey MLE
    ## 
    ## Estimated proportion: 2.4091%
    ## Standard error      : 0.4097%
    ## 
    ## Confidence intervals at the 95% level:
    ## Asymptotic Approach: 1.6060% - 3.2122%
    ## Clopper-Pearson    : 1.6495% - 3.3070%
    ## 
    ## Assumed measurement error: alpha = 1%, beta = 10% 
    ## Sampling: Random

## Stratified sampling

In many applications sampling is not uniformly random, but stratified.
This is also the case for the COVID-19 data from the Austrian survey
sample that we apply our method to. For such cases, one can use
different approaches. In this section, we present a method that
considers a prevalence estimator formed as a weighted sum of prevalence
estimators $`\widetilde{\pi}^k`$ associated to each stratum $`k`$, i.e.,
a generalization of the MME, as well as a Weighted M-Estimator (WME). In
both cases, we have to rely on asymptotic theory for computing CIs.

It actually turns out that the resulting estimators are based on similar
quantities provided previously, but weighted ones. Let $`\gamma_i`$
denote the sampling weight associated to subject $`i=1,\ldots,n`$, which
is proportional to the reciprocal of the sampling probability for
subject $`i`$, and adjusted such that $`\sum_{i=1}^n \gamma_i = n`$. Let
also

``` math
\begin{eqnarray}
\overline{R}_{11}&=&\sum_{i=1}^n\gamma_iY_i Z_i=\sum_{i=1}^n\gamma_iR_{i11}, \nonumber \\
\overline{R}_{10}&=&\sum_{i=1}^n\gamma_i(1-Y_i)Z_i=\sum_{i=1}^n\gamma_iR_{i10}, \nonumber \\
\overline{R}_{01}&=& \sum_{i=1}^n\gamma_iY_i(1-Z_i)= \sum_{i=1}^n\gamma_iR_{i01}, \nonumber\\ 
\overline{R}_{00}&=& \sum_{i=1}^n\gamma_i(1-Y_i)(1-Z_i)=\sum_{i=1}^n\gamma_iR_{i00}, \nonumber \\
\overline{R}_{*1}&=&\sum_{i=1}^n\gamma_iY_i=\overline{R}_{11}+\overline{R}_{01}.     \nonumber
\end{eqnarray}
```

With the MME approach, we consider the possibility that different groups
of people (such as in different towns or provinces), or even each
participant $`i`$, are associated to different prevalence $`\pi_i`$. We
are, however, only interested in the overall prevalence
$`\pi=(1/n)\sum_{i=1}^n \gamma_i \pi_i`$, given the additional
information provided by the official procedure. Moreover, the (known but
biased) prevalence from the official procedure could also be different
for each (groups of) subject(s) $`i`$, with $`\pi_{i0}, i=1,\ldots,n`$.
Note, however, that $`\pi_0 =(1/n) \sum_{i=1}^n \gamma_i \pi_{i0}`$ and
that $`\pi_0`$ is known. A general approach then consists in obtaining
estimates for each $`\pi_i`$ using some method and then take their
weighted average as an estimator for $`\pi`$. Using the MME approach
based on the $`R_{i01}`$ and assuming, as done before, that the
parameters $`\alpha,\beta`$ are specific to the medical test used and,
thus, independent of subject $`i`$, we obtain a weighted MME for $`\pi`$
as

``` math
\begin{equation}
    \widetilde{\pi} = \frac{1}{\Delta(1-\alpha_0)} \left( \frac{\overline{R}_{01}}{n} + \pi_0 (1 - \beta) - \alpha_0 \Delta - \alpha\right).
    \label{eqn:CWMLE-alpha0}
\end{equation}
```

We also have that $`\mathbb{E}\left[\widetilde{\pi}\right]= \pi`$. When
$`\alpha_0=\alpha=\beta=0`$, we get

``` math
\begin{equation*}
   \widetilde{\pi} = \pi_0 + \frac{\overline{R}_{01}}{n}.
\end{equation*}
```

An alternative approach is to consider a WME $`\widehat{\pi}`$ as
proposed for example by Wooldridge (2001). Generally, it has no
closed-form solution but can be computed numerically. However, in the
case when $`\alpha_0 = 0`$, we obtain a closed-form solution given by

``` math
\begin{equation*}
    \widehat{\pi} = \frac{\pi_0 \overline{R}_{00} + \overline{R}_{01}}{\Delta \left(\overline{R}_{01}+\overline{R}_{00}\right)} - \frac{\pi_0 \beta}{\Delta} - \frac{\alpha}{\Delta}.
\end{equation*}
```

When $`\alpha_0=\alpha=\beta=0`$, this further reduces to

``` math
\begin{equation*}
   \widehat{\pi} = \pi_0 \frac{n - \overline{R}_{\ast 1}}{n - \overline{R}_{11}} + \frac{R_{01}}{ \left(n - R_{11}\right)}.
\end{equation*}
```

In other words, in the case of stratified sampling, one replaces the
$`R_{jk}`$ by their weighted counterparts $`\overline{R}_{jk}`$.

Considering the data used in Guerrier et al. (2024), these estimators
can be computed as follows:

``` r

# Load pempi
library(pempi)

# Austrian data (November 2020)
pi0 = 93914/7166167

# Weighted sampling
R1w = sum(covid19_austria$weights[covid19_austria$Y == 1 & covid19_austria$Z == 1])
R2w = sum(covid19_austria$weights[covid19_austria$Y == 0 & covid19_austria$Z == 1])
R3w = sum(covid19_austria$weights[covid19_austria$Y == 1 & covid19_austria$Z == 0])
R4w = sum(covid19_austria$weights[covid19_austria$Y == 0 & covid19_austria$Z == 0])

# Average of squared weights 
V = mean(covid19_austria$weights^2)

# Compute CMLE
conditional_mle(R1 = R1w, R2 = R2w, R3 = R3w, R4 = R4w, 
                pi0 = pi0, V = V)
```

    ## Method: Conditional MLE
    ## 
    ## Estimated proportion: 2.9841%
    ## Standard error      : 0.3294%
    ## 
    ## Confidence interval at the 95% level:
    ## Asymptotic Approach: 2.3385% - 3.6297%
    ## 
    ## Assumed measurement error: alpha  = 0%, beta = 0%,
    ##                            alpha0 = 0% 
    ## 
    ## Estimated false negative rate of the
    ## official procedure: beta0 = 56.08%
    ## CI at the 95% level: 46.58% - 65.59%
    ## 
    ## Estimated ascertainment rate: 
    ## pi0/pi = 43.92%
    ## CI at the 95% level: 34.41% - 53.42%
    ## 
    ## Sampling: Stratified with V = 1.51

``` r

# Compute MME
moment_estimator(R3 = R3w, pi0 = pi0, n = n, V = V)
```

    ## Method: Moment Estimator
    ## 
    ## Estimated proportion: 2.9818%
    ## Standard error      : 0.3292%
    ## 
    ## Confidence interval at the 95% level:
    ## Asymptotic Approach: 2.3365% - 3.6270%
    ## 
    ## Assumed measurement error: alpha  = 0%, beta = 0%,
    ##                            alpha0 = 0% 
    ## 
    ## Estimated false negative rate of the
    ## official procedure: beta0 = 56.05%
    ## CI at the 95% level: 46.54% - 65.56%
    ## 
    ## Estimated ascertainment rate: 
    ## pi0/pi = 43.95%
    ## CI at the 95% level: 34.44% - 53.46%
    ## 
    ## Sampling: Stratified with V = 1.51

``` r

# Survey MLE
survey_mle(R = R1w + R3w, pi0 = pi0, n = n, V = V)
```

    ## Method: Survey MLE
    ## 
    ## Estimated proportion: 3.1280%
    ## Standard error      : 0.4471%
    ## 
    ## Confidence interval at the 95% level:
    ## Asymptotic Approach: 2.2517% - 4.0042%
    ## 
    ## Assumed measurement error: alpha = 0%, beta = 0% 
    ## Sampling: Stratified with V = 1.51

These estimators are also defined in the stratified case with
measurement error. For example, for the MME:

``` r

# Compute MME
moment_estimator(R3 = R3w, pi0 = pi0, n = n, V = V, alpha = alpha,
                 alpha0 = alpha0, beta = beta)
```

    ## Method: Moment Estimator
    ## 
    ## Estimated proportion: 2.0794%
    ## Standard error      : 0.3699%
    ## 
    ## Confidence interval at the 95% level:
    ## Asymptotic Approach: 1.3544% - 2.8045%
    ## 
    ## Assumed measurement error: alpha  = 1%, beta = 10%,
    ##                            alpha0 = 0% 
    ## 
    ## Estimated false negative rate of the
    ## official procedure: beta0 = 36.98%
    ## CI at the 95% level: 15.00% - 58.95%
    ## 
    ## Estimated ascertainment rate: 
    ## pi0/pi = 63.02%
    ## CI at the 95% level: 41.05% - 85.00%
    ## 
    ## Sampling: Stratified with V = 1.51

## References

Guerrier, Stéphane, Christoph Kuzmics, and Maria-Pia Victoria-Feser.
2024. “Assessing COVID-19 Prevalence in Austria with Infection Surveys
and Case Count Data as Auxiliary Information.” *Journal of the American
Statistical Association* 119 (547): 1722–35.
<https://doi.org/10.1080/01621459.2024.2313790>.

Wooldridge, J. M. 2001. “Asymptotic Properties of Weighted M-Estimators
for Standard Stratified Samples.” *Econometric Theory* 17: 451–70.
