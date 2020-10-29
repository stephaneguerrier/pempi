#' @title Compute proportion in the survey sample (standard estimator)
#' @description Proportion estimated using the survey sample and confidence intervals based on the Clopper-Pearson and the standard
#' asymptotic approach.
#' @param R        A \code{numeric} that provides the people of positive people in the sample.
#' @param n        A \code{numeric} that provides the sample size.
#' @param pi0      A \code{numeric} that provides the prevalence or proportion of people (in the whole population) who are positive, as measured through a non-random,
#' but systematic sampling (e.g. based on medical selection). Default value is \code{0} and in this case this information is not used in the estimation procedure.
#' @param alpha    A \code{numeric} that provides the False Negative (FN) rate for the sample R. Default value is \code{0}.
#' @param beta     A \code{numeric} that provides the False Positive (FP) rate for the sample R. Default value is \code{0}.
#' @param gamma    A \code{numeric} that used to compute a (1 - gamma) confidence region for the proportion. Default value is \code{0.05}.
#' @param ...      Additional arguments.
#' @return A \code{cpreval} object with the structure:
#' \itemize{
#'  \item estimate:    Estimated proportion.
#'  \item sd:          Estimated standard error of the estimator.
#'  \item ci_asym:     Asymptotic confidence interval at the 1 - gamma confidence level.
#'  \item gamma:       Confidence level (i.e. 1 - gamma) for confidence intervals.
#'  \item method:      Estimation method (in this case sample survey).
#'  \item measurement: A vector with (alpha0, alpha, beta).
#'  \item beta0        Estimated false negative rate of the official procedure.
#'  \item ...:         Additional parameters.
#' }
#' @export
#' @author Stephane Guerrier, Maria-Pia Victoria-Feser, Christoph Kuzmics
#' @examples
#' # Samples without measurement error
#' X = sim_Rs(theta = 30/1000, pi0 = 10/1000, n = 1500, seed = 18)
#' survey_mle(R = X$R, n = X$n)
#'
#' # With measurement error
#' X = sim_Rs(theta = 30/1000, pi0 = 10/1000, n = 1500, alpha = 0.01, beta = 0.05, seed = 18)
#' survey_mle(R = X$R, n = X$n)
#' survey_mle(R = X$R, n = X$n, alpha = 0.01, beta = 0.05)
#' @importFrom stats qbeta qnorm
survey_mle = function(R, n, pi0 = 0, alpha = 0, beta = 0, gamma = 0.05, ...){
  # Check inputs
  if (R%%1!=0 || R < 0){
    stop("R should be non-negative integer.")
  }

  if (n%%1!=0 || n < 0){
    stop("n should be non-negative integer.")
  }

  if (R > n){
    stop("R and n should be such that R <= n.")
  }

  if (pi0 < 0 || pi0 > 1){
    stop("pi0 should be such that 0 <= pi0 <= 1.")
  }

  if (alpha < 0 || alpha >= 1 || beta < 0 || beta >= 1 || alpha + beta > 1){
    stop("alpha and beta should be such that: (i) 0 <= alpha < 1, (ii) 0 <= beta < 1, (iii) alpha + beta < 1.")
  }

  if (gamma < 0 || gamma > 1){
    stop("gamma should be between 0 and 1.")
  }

  # Compute survey proportion
  pi_bar = R/n

  # Adjust for false positive/negative
  if (max(alpha, beta) > 0){
    pi_bar = (pi_bar - alpha)/(1 - alpha - beta)
  }

  if (pi_bar < pi0 || pi_bar > 1){
    pi_bar = c(pi0, 1)[which.min(abs(pi_bar - c(pi0, 1)))]
    sd = NA
    ci_asym = NA

    if (pi_bar == pi0){
      upper = (qbeta(p = 1 - gamma/2, R + 1, n - R) - alpha)/(1 - alpha - beta)
      lower = pi0
    }else{
      upper = 1
      lower = (qbeta(p = gamma/2, R, n - R + 1) - alpha)/(1 - alpha - beta)
    }

    ci_cp = c(lower, upper)

    # Construct output
    out = list(estimate = pi_bar, sd = sd, ci_asym = ci_asym, ci_cp = ci_cp, gamma = gamma,
               method = "Survey MLE", measurement = c(NA, alpha, NA, beta),
               boundary = TRUE,...)
    class(out) = "cpreval"
    return(out)
  }

  # Estimated standard error
  Delta = 1 - alpha - beta
  sd = sqrt((pi_bar*Delta + alpha)*(1 - pi_bar*Delta - alpha)/(n*Delta^2))

  # Compute 1 - gamma asymptotic interval
  ci_asym = pi_bar + qnorm(1 - gamma/2)*c(-1,1)*sd

  # Compute 1 - gamma confidence interval - Clopper-Pearson approach
  upper = (qbeta(p = 1 - gamma/2, R + 1, n - R) - alpha)/(1 - alpha - beta)
  lower = (qbeta(p = gamma/2, R, n - R + 1) - alpha)/(1 - alpha - beta)

  # Adjust CI - CP
  if (lower < pi0){
    lower = pi0
  }

  if (upper > 1){
    upper = 1
  }

  ci_cp = c(lower, upper)

  # Construct output
  out = list(estimate = pi_bar, sd = sd, ci_asym = ci_asym, ci_cp = ci_cp, gamma = gamma,
             method = "Survey MLE", measurement = c(NA, alpha, beta), beta0 = NA)
  class(out) = "cpreval"
  out
}

#' @title Print estimation results
#' @description Simple print function for the four estimators (i.e. mle, conditional mle, moment based and survey sample).
#' @method print cpreval
#' @export
#' @keywords internal
#' @param x    A \code{cpreval} object
#' @param ...  Further arguments passed to or from other methods
#' @return Prints object
#' @author Stephane Guerrier
#' @examples
#' X = sim_Rs(theta = 3/100, pi0 = 1/100, n = 1500, seed = 18)
#' survey_mle(X$R, X$n)
print.cpreval = function(x, ...){
  cat("Method: ")
  cat(x$method)
  cat("\n\n")
  cat("Estimated proportion: ")
  cat(sprintf("%.4f", 100*x$estimate))
  cat("%\n")
  cat("Standard error      : ")
  cat(sprintf("%.4f", 100*x$sd))
  cat("%\n\n")
  cat("Confidence intervals at the ")
  cat(100*(1 - x$gamma))
  cat("% level:\n")
  cat("Asymptotic Approach: ")
  cat(sprintf("%.4f", 100*x$ci_asym[1]))
  cat("% - ")
  cat(sprintf("%.4f", 100*x$ci_asym[2]))
  if (x$method == "Survey MLE" || x$method == "Moment Estimator"){
    cat("%\n")
    cat("Clopper-Pearson    : ")
    cat(sprintf("%.4f", 100*x$ci_cp[1]))
    cat("% - ")
    cat(sprintf("%.4f", 100*x$ci_cp[2]))
  }
  cat("%\n\n")

  if (x$method == "Survey MLE"){
    cat("Assumed measurement error: alpha = ")
    cat(100*x$measurement[2])
    cat("%, beta = ")
    cat(100*x$measurement[3])
    cat("% \n")
  }else{
    cat("Assumed measurement error: alpha  = ")
    cat(100*x$measurement[2])
    cat("%, beta = ")
    cat(100*x$measurement[3])
    cat("%,\n")
    cat("                           alpha0 = ")
    cat(100*x$measurement[1])
    cat("% \n\n")

    cat("Estimated false negative rate of the\n")
    cat("official procedure: beta0 = ")
    cat(round(100*x$beta0,4))
    cat("%\n")
  }
}

#' @title Compute moment-based estimator.
#' @description Proportion estimated using the moment-based estimator and confidence intervals based the asymptotic distribution of the estimator as well as
#' the Clopper-Pearson approach.
#' @param R3        A \code{numeric} that provides the number of participants in the survey sample that are tested positive only with the second testing device.
#' @param n         A \code{numeric} that provides the sample size.
#' @param pi0       A \code{numeric} that provides the prevalence or proportion of people (in the whole population) who are positive, as measured through a non-random,
#' but systematic sampling (e.g. based on medical selection).
#' @param alpha0       A \code{numeric} that corresponds to the probability that a random participant
#' has been incorrectly declared positive through the nontransparent procedure. In most applications,
#' this probability is likely very close to zero. Default value is \code{0}.
#' @param alpha     A \code{numeric} that provides the False Negative (FN) rate for the sample R. Default value is \code{0}.
#' @param beta      A \code{numeric} that provides the False Positive (FP) rate for the sample R. Default value is \code{0}.
#' @param gamma     A \code{numeric} that used to compute a (1 - gamma) confidence region for the proportion. Default value is \code{0.05}.
#' @param ...       Additional arguments.
#' @return A \code{cpreval} object with the structure:
#' \itemize{
#'  \item estimate:    Estimated proportion.
#'  \item sd:          Estimated standard error of the estimator.
#'  \item ci_asym:     Asymptotic confidence interval at the 1 - gamma confidence level.
#'  \item ci_cp:       Confidence interval (1 - gamma confidence level) based on the Clopper-Pearson approach.
#'  \item gamma:       Confidence level (i.e. 1 - gamma) for confidence intervals.
#'  \item method:      Estimation method (in this case moment estimator).
#'  \item measurement: A vector with (alpha0, alpha, beta).
#'  \item beta0        Estimated false negative rate of the official procedure.
#'  \item ...:         Additional parameters.
#' }
#' @export
#' @author Stephane Guerrier, Maria-Pia Victoria-Feser, Christoph Kuzmics
#' @examples
#' # Samples without measurement error
#' X = sim_Rs(theta = 3/100, pi0 = 1/100, n = 1500, seed = 18)
#' moment_estimator(R3 = X$R3, n = X$n, pi0 = X$pi0)
#'
#' # With measurement error
#' X = sim_Rs(theta = 3/100, pi0 = 1/100, n = 1500, alpha0 = 0.001,
#' alpha = 0.01, beta = 0.05, seed = 18)
#' moment_estimator(R3 = X$R3, n = X$n, pi0 = X$pi0)
#' moment_estimator(R3 = X$R3, n = X$n, pi0 = X$pi0, alpha0 = 0.001,
#' alpha = 0.01, beta = 0.05)
#' @importFrom stats qbeta qnorm
moment_estimator = function(R3, n, pi0, gamma = 0.05, alpha = 0, beta = 0, alpha0 = 0, ...){
  # Check inputs
  if (R3%%1!=0 || R3 < 0){
    stop("R3 should be non-negative integer.")
  }

  if (alpha0 > pi0){
    stop("The inputs pi0 and alpha0 must be such that alpha0 <= pi0.")
  }

  if (n%%1!=0 || n < 0){
    stop("n should be non-negative integer.")
  }

  if (R3 > n){
    stop("R3 and n should be such that R3 <= n.")
  }

  if (pi0 < 0 || pi0 > 1){
    stop("pi0 should be such that 0 <= pi0 <= 1.")
  }

  if (alpha < 0 || alpha >= 1 || beta < 0 || beta >= 1 || alpha + beta > 1){
    stop("alpha and beta should be such that: (i) 0 <= alpha < 1, (ii) 0 <= beta < 1, (iii) alpha + beta < 1.")
  }

  if (alpha0 < 0 || alpha0 >= 1){
    stop("alpha0 should be close to 0 and such that: 0 <= alpha0 < 1.")
  }

  if (gamma < 0 || gamma > 1){
    stop("gamma should be between 0 and 1.")
  }

  # Compute point estimate
  beta0 = 0 # This is not used (TO DO remove beta0 here)
  Delta = 1 - alpha - beta
  Delta0 = 1 - alpha0 - beta0
  pi0_star = (pi0 - alpha0)/Delta0
  estimate = 1/(Delta*(1 - alpha0))*(R3/n + pi0 - beta*pi0 - alpha0*Delta - alpha)

  # Asymptotic CI (1-gamma)
  probs = get_prob(theta = estimate, pi0 = pi0, alpha = alpha,
                   beta = beta, alpha0 = alpha0)

  sd = sqrt((probs[3]*(1 - probs[3]))/(n*Delta^2*(1 - alpha0)^2))
  ci_asym = estimate + qnorm(1 - gamma/2)*c(-1,1)*sd

  # CP - CI (1 - gamma)
  upper = (qbeta(p = 1 - gamma/2, R3 + 1, n - R3) - alpha*(1 - alpha0) + (pi0 - alpha0)*(1 - beta))/(Delta*(1 - alpha0))
  lower = (qbeta(p = gamma/2, R3, n - R3 + 1) - alpha*(1 - alpha0) + (pi0 - alpha0)*(1 - beta))/(Delta*(1 - alpha0))

  # Adjust CI - CP
  if (lower < pi0){
    lower = pi0
  }

  if (upper > 1){
    upper = 1
  }
  ci_cp = c(lower, upper)

  # Compute estimated false negative rate of the official procedure
  beta0 = 1 - (pi0 - alpha0*(1 - estimate))/estimate

  # Construct output
  out = list(estimate = estimate, sd = sd, ci_asym = ci_asym, ci_cp = ci_cp, gamma = gamma,
             method = "Moment Estimator", measurement = c(alpha0, alpha, beta), beta0 = beta0)
  class(out) = "cpreval"
  out
}


#' @title Modified log function (internal function)
#' @description Simple function that return log(x) if x > 0 and 0 otherwise.
#' @param x     A \code{numeric} value assumed such that x >= 0.
#' @return Modified log function
#' @author Stephane Guerrier
log_modified = function(x){
  if (x == 0){
    return(0)
  }else{
    return(log(x))
  }
}

#' @title Negative Log-Likelihood function
#' @description Log-Likelihood function based on R1, R2, R3 and R4 multiplied by -1.
#' @param theta     A \code{numeric} value for the paramerer of interest.
#' @param Rvect     A \code{vector} of observations, i.e. (R1, R2, R3, R4).
#' @param n         A \code{numeric} that provides the sample size.
#' @param pi0       A \code{numeric} that provides the prevalence or proportion of people (in the whole population) who are positive, as measured through a non-random,
#' but systematic sampling (e.g. based on medical selection).
#' @param alpha     A \code{numeric} that provides the False Negative (FN) rate for the sample R.
#' @param beta      A \code{numeric} that provides the False Positive (FP) rate for the sample R.
#' @param alpha0       A \code{numeric} that corresponds to the probability that a random participant
#' has been incorrectly declared positive through the nontransparent procedure. In most applications,
#' this probability is likely very close to zero.
#' @param ...       Additional arguments.
#' @return Negative log-likelihood.
#' @author Stephane Guerrier
neg_log_lik = function(theta, Rvect, n, pi0, alpha, beta, alpha0, ...){
  probs = get_prob(theta = theta, pi0 = pi0, alpha = alpha, beta = beta, alpha0 = alpha0)
  (-1)*(Rvect[1]/n*log_modified(probs[1]) + Rvect[2]/n*log_modified(probs[2]) + Rvect[3]/n*log_modified(probs[3]) + Rvect[4]/n*log_modified(probs[4]))
}


#' @title Negative marginalized log-likelihood function based on R0 and R
#' @description Log-Likelihood function based on R1 and R3 multiplied by -1.
#' @param theta     A \code{numeric} value for the paramerer of interest.
#' @param Rvect     A \code{vector} of observations, i.e. (R1, R2, R3, R4), but R2 and R4 are NOT used.
#' @param n         A \code{numeric} that provides the sample size.
#' @param pi0       A \code{numeric} that provides the prevalence or proportion of people (in the whole population) who are positive, as measured through a non-random,
#' but systematic sampling (e.g. based on medical selection).
#' @param alpha0       A \code{numeric} that corresponds to the probability that a random participant
#' has been incorrectly declared positive through the nontransparent procedure. In most applications,
#' this probability is likely very close to zero.
#' @param alpha     A \code{numeric} that provides the False Negative (FN) rate for the sample R.
#' @param beta      A \code{numeric} that provides the False Positive (FP) rate for the sample R.
#' @param ...       Additional arguments.
#' @return Negative marginalized log-likelihood.
#' @author Stephane Guerrier
neg_log_lik_integrated = function(theta, Rvect, n, pi0, alpha0, alpha, beta, ...){
  probs = get_prob(theta = theta, pi0 = pi0, alpha = alpha, beta = beta, alpha0 = alpha0)
  (-1)*(Rvect[1]/n*log_modified(probs[1]) + probs[2]/n*log_modified(probs[2]) + Rvect[3]/n*log_modified(probs[3]) + (n - Rvect[1] - Rvect[3] - probs[2])/n*log_modified(probs[4]))
}


#' @title Compute MLE based on the full information R1, R2, R3 and R4.
#' @description Proportion estimated using the MLE and confidence intervals based the asymptotic distribution of the estimator.
#' @param R1        A \code{numeric} that provides the number of participants in the survey sample that were tested positive with both (medical) testing devices (and are, thus, members of the sub-population).
#' @param R2        A \code{numeric} that provides the number of participants in the survey sample that are tested positive only with the first testing device (and are, thus,  members of the sub-population).
#' @param R3        A \code{numeric} that provides the number of participants in the survey sample that are tested positive only with the second testing device.
#' @param R4        A \code{numeric} that provides the number of participants that are tested negative with the second testing device (and are either members of the sub-population and have tested negative with the first testing device or are not members of the sub-population).
#' @param n         A \code{numeric} that provides the sample size. Default value R1 + R2 + R3 + R4. If this value is provided it is used to verify that R1 + R2 + R3 + R4 = n.
#' @param pi0       A \code{numeric} that provides the prevalence or proportion of people (in the whole population) who are positive, as measured through a non-random,
#' but systematic sampling (e.g. based on medical selection).
#' @param gamma     A \code{numeric} that used to compute a (1 - gamma) confidence region for the proportion. Default value is \code{0.05}.
#' @param alpha0       A \code{numeric} that corresponds to the probability that a random participant
#' has been incorrectly declared positive through the nontransparent procedure. In most applications,
#' this probability is likely very close to zero. Default value is \code{0}.
#' @param alpha     A \code{numeric} that provides the False Negative (FN) rate for the sample R. Default value is \code{0}.
#' @param beta      A \code{numeric} that provides the False Positive (FP) rate for the sample R. Default value is \code{0}.
#' @param ...       Additional arguments.
#' @return A \code{cpreval} object with the structure:
#' \itemize{
#'  \item estimate:    Estimated proportion.
#'  \item sd:          Estimated standard error of the estimator.
#'  \item ci_asym:     Asymptotic confidence interval at the 1 - gamma confidence level.
#'  \item gamma:       Confidence level (i.e. 1 - gamma) for confidence intervals.
#'  \item method:      Estimation method (in this case mle).
#'  \item measurement: A vector with (alpha0, alpha, beta).
#'  \item beta0:       Estimated false negative rate of the official procedure.
#'  \item ...:         Additional parameters.
#' }
#' @export
#' @author Stephane Guerrier, Maria-Pia Victoria-Feser, Christoph Kuzmics
#' @examples
#' # Samples without measurement error
#' X = sim_Rs(theta = 3/100, pi0 = 1/100, n = 1500, seed = 18)
#' conditional_mle(R1 = X$R1, R2 = X$R2, R3 = X$R3, R4 = X$R4, pi0 = X$pi0)
#'
#' # With measurement error
#' X = sim_Rs(theta = 30/1000, pi0 = 10/1000, n = 1500, alpha0 = 0.001,
#' alpha = 0.01, beta0 = 0.05, beta = 0.05, seed = 18)
#' conditional_mle(R1 = X$R1, R2 = X$R2, R3 = X$R3, R4 = X$R4, pi0 = X$pi0)
#' conditional_mle(R1 = X$R1, R2 = X$R2, R3 = X$R3, R4 = X$R4, pi0 = X$pi0,
#' alpha0 = 0.001, alpha = 0.01, beta = 0.05)
#' @importFrom stats optimize qnorm
conditional_mle = function(R1 = NULL, R2 = NULL, R3 = NULL, R4 = NULL, n = R1 + R2 + R3 + R4, pi0, gamma = 0.05, alpha0 = 0, alpha = 0, beta = 0, ...){
  # Check inputs
  if (is.null(R1) + is.null(R2) + is.null(R3) + is.null(R4) == 1 && length(n) == 0){
    if (is.null(R1)){
      stop("R1 is missing.")
    }

    if (alpha0 > pi0){
      stop("The inputs pi0 and alpha0 must be such that alpha0 <= pi0.")
    }

    if (is.null(R2)){
      stop("R2 is missing.")
    }

    if (is.null(R3)){
      stop("R3 is missing.")
    }

    if (is.null(R4)){
      stop("R4 is missing.")
    }
  }

  if (is.null(R1) + is.null(R2) + is.null(R3) + is.null(R4) + is.null(n) == 1){
    if (is.null(R1)){
      R1 = n - R2 - R3 - R4
      warning(paste("R1 was computed using the other values: R1 = ", R1, sep = ""))
    }

    if (is.null(R2)){
      R2 = n - R1 - R3 - R4
      warning(paste("R2 was computed using the other values: R2 = ", R2, sep = ""))
    }

    if (is.null(R3)){
      R3 = n - R1 - R2 - R4
      warning(paste("R3 was computed using the other values: R3 = ", R3, sep = ""))
    }

    if (is.null(R4)){
      R4 = n - R1 - R2 - R3
      warning(paste("R4 was computed using the other values: R4 = ", R4, sep = ""))
    }
  }

  if (R1%%1!=0 || R1 < 0){
    stop("R1 should be non-negative integer.")
  }

  if (R2%%1!=0 || R2 < 0){
    stop("R2 should be non-negative integer.")
  }

  if (R3%%1!=0 || R3 < 0){
    stop("R3 should be non-negative integer.")
  }

  if (R4%%1!=0 || R4 < 0){
    stop("R4 should be non-negative integer.")
  }

  if (n%%1!=0 || n < 0){
    stop("n should be non-negative integer.")
  }

  if (R1 + R2 + R3 + R4 != n){
    stop("The sum of R_i should be equal to n.")
  }

  if (pi0 < 0 || pi0 > 1){
    stop("pi0 should be such that 0 <= pi0 <= 1.")
  }

  if (alpha < 0 || alpha >= 1 || beta < 0 || beta >= 1 || alpha + beta > 1){
    stop("alpha and beta should be such that: (i) 0 <= alpha < 1, (ii) 0 <= beta < 1, (iii) alpha + beta < 1.")
  }

  if (alpha0 < 0 || alpha0 >= 1){
    stop("alpha0 should be close to 0 and such that 0 <= alpha0 < 1.")
  }

  if (gamma < 0 || gamma > 1){
    stop("gamma should be between 0 and 1.")
  }

  # Find MLE (TODO use closed form when possible)
  R = c(R1, R2, R3, R4)
  eps = 10^(-5)
  estimate = optimize(neg_log_lik, interval = c(pi0 + eps, 1 - eps), R = R, n = n, pi0 = pi0, alpha = alpha, beta = beta, alpha0 = alpha0)$minimum

  # Compute implied probs
  probs = get_prob(theta = estimate, pi0 = pi0, alpha = alpha, beta = beta, alpha0 = alpha0)

  # Deltas
  beta0 = 0 # This is not used (TO DO remove beta0 here)
  Delta = 1 - alpha - beta
  Delta0 = 1 - alpha0 - beta0

  # Compute Fisher Info
  a1 = ((1 - alpha0)^2*Delta^2*probs[4])/((1 - alpha0)*(1 - estimate*Delta - alpha) - beta*(pi0 - alpha0))^2
  a2 = (alpha0^2*Delta^2*probs[1])/(alpha*alpha0 + (1-beta)*(pi0 - alpha0) + estimate*alpha0*Delta)^2
  a3 = (alpha0^2*Delta^2*probs[2])/(alpha0*(1 - estimate*Delta - alpha) + beta*(pi0 - alpha0))^2
  a4 = ((1 - alpha0)^2*Delta^2*probs[3])/((1-beta)*(pi0 - alpha0) - estimate*(1-alpha0)*Delta - alpha*(1-alpha0))^2

  fisher_info = 0

  if (!is.na(a1))
    fisher_info = fisher_info + a1

  if (!is.na(a2))
    fisher_info = fisher_info + a2

  if (!is.na(a3))
    fisher_info = fisher_info + a3

  if (!is.na(a4))
    fisher_info = fisher_info + a4

  # Asymptotic variance
  mle_var = 1/(n*fisher_info)

  # Compute Asymptotic CI
  sd = sqrt(mle_var)
  ci_asym = estimate + qnorm(1 - gamma/2)*c(-1,1)*sd

  # Compute estimated false negative rate of the official procedure
  beta0 = 1 - (pi0 - alpha0*(1 - estimate))/estimate

  # Construct output
  out = list(estimate = estimate, sd = sd, ci_asym = ci_asym, gamma = gamma,
             method = "Conditional MLE", measurement = c(alpha0, alpha, beta), beta0 = beta0, ...)
  class(out) = "cpreval"
  out
}

#' @title Compute (marginalized) MLE based on the partial information R1 and R3.
#' @description Proportion estimated using the MLE and confidence intervals based the asymptotic distribution of the estimator.
#' @param R1        A \code{numeric} that provides the number of participants in the survey sample that were tested positive with both (medical) testing devices (and are, thus, members of the sub-population).
#' @param R3        A \code{numeric} that provides the number of participants in the survey sample that are tested positive only with the second testing device.
#' @param n         A \code{numeric} that provides the sample size.
#' @param pi0       A \code{numeric} that provides the prevalence or proportion of people (in the whole population) who are positive, as measured through a non-random,
#' but systematic sampling (e.g. based on medical selection).
#' @param gamma     A \code{numeric} that used to compute a (1 - gamma) confidence region for the proportion. Default value is \code{0.05}.
#' @param alpha0       A \code{numeric} that corresponds to the probability that a random participant
#' has been incorrectly declared positive through the nontransparent procedure. In most applications,
#' this probability is likely very close to zero. Default value is \code{0}.
#' @param alpha     A \code{numeric} that provides the False Negative (FN) rate for the sample R. Default value is \code{0}.
#' @param beta      A \code{numeric} that provides the False Positive (FP) rate for the sample R. Default value is \code{0}.
#' @param ...       Additional arguments.
#' @return A \code{cpreval} object with the structure:
#' \itemize{
#'  \item estimate:    Estimated proportion.
#'  \item sd:          Estimated standard error of the estimator.
#'  \item ci_asym:     Asymptotic confidence interval at the 1 - gamma confidence level.
#'  \item gamma:       Confidence level (i.e. 1 - gamma) for confidence intervals.
#'  \item method:      Estimation method (in this case marginal mle).
#'  \item measurement: A vector with (alpha0, alpha, beta).
#'  \item beta0:       Estimated false negative rate of the official procedure.
#'  \item ...:         Additional parameters
#' }
#' @export
#' @author Stephane Guerrier, Maria-Pia Victoria-Feser, Christoph Kuzmics
#' @examples
#' # Samples without measurement error
#' X = sim_Rs(theta = 3/100, pi0 = 1/100, n = 1500, seed = 18)
#' conditional_mle(R1 = X$R1, R2 = X$R2, R3 = X$R3, R4 = X$R4, n = X$n, pi0 = X$pi0)
#'
#' # With measurement error
#' X = sim_Rs(theta = 30/1000, pi0 = 10/1000, n = 1500, alpha0 = 0.001,
#' alpha = 0.01, beta0 = 0.05, beta = 0.05, seed = 18)
#' marginal_mle(R1 = X$R1, R3 = X$R3, n = X$n, pi0 = X$pi0)
#' marginal_mle(R1 = X$R1, R3 = X$R3, n = X$n, pi0 = X$pi0,
#' alpha0 = 0.001, alpha = 0.01, beta0 = 0.05, beta = 0.05)
#' @importFrom stats optimize
marginal_mle = function(R1, R3, n, pi0, gamma = 0.05, alpha = 0, beta = 0, alpha0 = 0, ...){
  # Check inputs
  if (R1%%1!=0 || R1 < 0){
    stop("R1 should be non-negative integer.")
  }

  if (alpha0 > pi0){
    stop("The inputs pi0 and alpha0 must be such that alpha0 <= pi0.")
  }

  if (R3%%1!=0 || R3 < 0){
    stop("R3 should be non-negative integer.")
  }

  if (n%%1!=0 || n < 0){
    stop("n should be non-negative integer.")
  }

  if (R1 + R3 > n){
    stop("R1 + R3 can't be larger than n.")
  }

  if (pi0 < 0 || pi0 > 1){
    stop("pi0 should be such that 0 <= pi0 <= 1.")
  }

  if (alpha < 0 || alpha >= 1 || beta < 0 || beta >= 1 || alpha + beta > 1){
    stop("alpha and beta should be such that: (i) 0 <= alpha < 1, (ii) 0 <= beta < 1, (iii) alpha + beta < 1.")
  }

  if (alpha0 < 0 || alpha0 >= 1){
    stop("alpha0 should be close to 0 and such that 0 <= alpha0 < 1.")
  }

  if (gamma < 0 || gamma > 1){
    stop("gamma should be between 0 and 1.")
  }

  # Compute MLE (TODO use closed form when possible)
  R = c(R1, NA, R3, NA)
  eps = 10^(-5)
  estimate = optimize(neg_log_lik_integrated, interval = c(pi0 + eps, 1 - eps), R = R, n = n, pi0 = pi0, alpha = alpha, beta = beta, alpha0 = alpha0)$minimum

  # Compute probs
  probs = get_prob(theta = estimate, pi0 = pi0, alpha = alpha, beta = beta, alpha0 = alpha0)

  # Deltas
  beta0 = 0 # This is not used (TO DO remove beta0 here)
  Delta = 1 - alpha - beta
  Delta0 = 1 - alpha0 - beta0

  # Compute Fishier Info
  a1 =  -((alpha0 - 1)^2*(alpha + beta - 1)^2*(probs[1] + probs[3] - 1 - beta*(alpha0 - pi0) - alpha0*(alpha - 1) + alpha0*estimate*(alpha + beta - 1)))/(beta*(alpha0 - pi0) + (alpha - 1)*(alpha0 - 1) - estimate*(alpha0 - 1)*(alpha + beta - 1))^2
  a2 = (1/n)*(alpha0^2*(alpha + beta - 1)^2)/(beta*(alpha0 - pi0) + alpha0*(alpha - 1) - alpha0*estimate*(alpha + beta - 1))
  a3 = -(1/n)*(2*alpha0*(alpha0 - 1)*(alpha + beta - 1)^2)/(beta*(alpha0 - pi0) + (alpha - 1)*(alpha0 - 1) - estimate*(alpha0 - 1)*(alpha + beta - 1))
  a4 = (probs[1]*alpha0^2*(alpha + beta - 1)^2)/(alpha*alpha0 + (alpha0 - pi0)*(beta - 1) - alpha0*estimate*(alpha + beta - 1))^2
  a5 = (probs[3]*(alpha0 - 1)^2*(alpha + beta - 1)^2)/(alpha*(alpha0 - 1) + (alpha0 - pi0)*(beta - 1) - estimate*(alpha0 - 1)*(alpha + beta - 1))^2

  fisher_info = 0

  if (!is.na(a1))
    fisher_info = fisher_info + a1

  if (!is.na(a2))
    fisher_info = fisher_info + a2

  if (!is.na(a3))
    fisher_info = fisher_info + a3

  if (!is.na(a4))
    fisher_info = fisher_info + a4

  if (!is.na(a5))
    fisher_info = fisher_info + a5

  # Asymptotic variance
  mle_var = 1/(n*fisher_info)

  # Compute Asymptotic CI
  sd = sqrt(mle_var)
  ci_asym = estimate + qnorm(1 - gamma/2)*c(-1,1)*sd

  # Compute estimated false negative rate of the official procedure
  beta0 = 1 - (pi0 - alpha0*(1 - estimate))/estimate

  # Construct output
  out = list(estimate = estimate, sd = sd, ci_asym = ci_asym, gamma = gamma,
             method = "Marginal MLE", measurement = c(alpha0, alpha, beta), beta0 = beta0, ...)
  class(out) = "cpreval"
  out
}
