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
#' @param V        A \code{numeric} that corresponds to the average of squared sampling weights. Default value is \code{NULL}.
#' @param ...      Additional arguments.
#' @return A \code{cpreval} object with the structure:
#' \itemize{
#'  \item estimate:    Estimated proportion.
#'  \item sd:          Estimated standard error of the estimator.
#'  \item ci_asym:     Asymptotic confidence interval at the 1 - gamma confidence level.
#'  \item gamma:       Confidence level (i.e. 1 - gamma) for confidence intervals.
#'  \item method:      Estimation method (in this case sample survey).
#'  \item measurement: A vector with (alpha0, alpha, beta).
#'  \item boundary:    A boolean variable indicating if the estimates falls at the boundary of the parameter space.
#'  \item pi0:         Value of pi0 (input value).
#'  \item sampling:    Type of sampling considered ("random" or "weighted").
#'  \item V:           Average sum of squared sampling weights if weighted/stratified is used (otherwise NULL).
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
survey_mle = function(R, n, pi0 = 0, alpha = 0, beta = 0, gamma = 0.05, V = NULL, ...){
  # Check inputs
  if (R%%1!=0 || R < 0){
    if (is.null(V)){
      stop("R should be non-negative integer. To use sampling weights you need to provide a value for V.")
    }else{
      if (R < 0){
        stop("R should be non-negative.")
      }
    }
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

  # Check V
  if (!is.null(V)){
    if (V < 0){
      stop("V should be larger than 0")
    }
    if (V < 1){
      warning("V should be larger than 1.")
    }
  }

  # Compute survey proportion
  pi_bar = R/n

  # Adjust for false positive/negative
  if (max(alpha, beta) > 0){
    pi_bar = (pi_bar - alpha)/(1 - alpha - beta)
  }

  if (!is.null(V)){
    if ((pi_bar < pi0 || pi_bar > 1)){
      out = list(estimate = min(max(pi_bar, pi0), 1), sd = NA,
                 ci_asym = NA, ci_cp = NA, gamma = gamma,
                 method = "Survey MLE", measurement = c(NA, alpha, beta),
                 boundary = TRUE, pi0 = pi0, sampling = "weighted", V = V, ...)
      class(out) = "cpreval"
      return(out)
    }else{
      sd = sqrt((V/(1 - alpha - beta)^2*(alpha + pi_bar*(1 - alpha - beta))*(1 - alpha - pi_bar*(1 - alpha - beta)))/n)
      ci_asym = pi_bar + c(-1, 1)*qnorm(1 - gamma/2)*sd

      # Construct output
      out = list(estimate = pi_bar, sd = sd, ci_asym = ci_asym, ci_cp = NA, gamma = gamma,
                 method = "Survey MLE", measurement = c(NA, alpha, beta),
                 boundary = FALSE, pi0 = pi0, sampling = "weighted", V = V, ...)
      class(out) = "cpreval"
      return(out)
    }
  }

  if (pi_bar < pi0 || pi_bar > 1){
    pi_bar = c(pi0, 1)[which.min(abs(pi_bar - c(pi0, 1)))]
    sd = NA
    ci_asym = c(NA, NA)
    boundary = TRUE

    if (pi_bar == pi0){
      upper = (qbeta(p = 1 - gamma, R + 1, n - R) - alpha)/(1 - alpha - beta)
      lower = pi0
    }else{
      upper = 1
      lower = (qbeta(p = gamma, R, n - R + 1) - alpha)/(1 - alpha - beta)
    }

    ci_cp = c(lower, upper)

    # Construct output
    out = list(estimate = pi_bar, sd = sd, ci_asym = ci_asym, ci_cp = ci_cp, gamma = gamma,
               method = "Survey MLE", measurement = c(NA, alpha, beta),
               boundary = boundary, pi0 = pi0, sampling = "random", V = V, ...)
    class(out) = "cpreval"
    return(out)
  }else{
    boundary = FALSE
  }

  # Estimated standard error
  Delta = 1 - alpha - beta
  sd = sqrt((pi_bar*Delta + alpha)*(1 - pi_bar*Delta - alpha)/(n*Delta^2))

  # Compute 1 - gamma asymptotic interval
  ci_asym = pi_bar + qnorm(1 - gamma/2)*c(-1,1)*sd

  # Compute 1 - gamma confidence interval - Clopper-Pearson approach
  upper = (qbeta(p = 1 - gamma/2, R + 1, n - R) - alpha)/(1 - alpha - beta)
  lower = (qbeta(p = gamma/2, R, n - R + 1) - alpha)/(1 - alpha - beta)
  ci_cp = c(lower, upper)

  # Construct output
  out = list(estimate = pi_bar, sd = sd, ci_asym = ci_asym, ci_cp = ci_cp, gamma = gamma,
             method = "Survey MLE", measurement = c(NA, alpha, beta),
             boundary = boundary, pi0 = pi0, sampling = "random", V = V, ...)
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
  cat("Confidence interval")
  if (x$sampling == "random" && (x$method == "Survey MLE" || x$method == "Moment Estimator")){
    cat("s")
  }
  cat(" at the ")
  cat(100*(1 - x$gamma))
  cat("% level:\n")
  cat("Asymptotic Approach: ")
  cat(sprintf("%.4f", 100*x$ci_asym[1]))
  cat("% - ")
  cat(sprintf("%.4f", 100*x$ci_asym[2]))
  if (x$sampling == "random" && (x$method == "Survey MLE" || x$method == "Moment Estimator")){
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
    cat(sprintf("%.2f", 100*x$beta0))
    cat("%\n")
    cat("CI at the ")
    cat(100*(1 - x$gamma))
    cat("% level: ")
    cat(sprintf("%.2f", 100*x$ci_beta0[1]))
    cat("% - ")
    cat(sprintf("%.2f", 100*x$ci_beta0[2]))
    cat("%\n")
  }

  cat("Sampling: ")
  if (x$sampling == "random"){
    cat("Random")
  }else{
    cat("Stratified with V = ")
    cat(sprintf("%.2f", x$V))
    cat("\n")
  }

  if (x$boundary){
    warning("Parameter estimates is in the boundary of the space. Results may not be reliable.")
  }

  if(x$sampling == "random" && (x$method == "Survey MLE" || x$method == "Moment Estimator")){
    if (x$ci_cp[1] > x$ci_cp[2] || x$ci_cp[1] < x$pi0 || x$ci_cp[2] > 1){
      warning("The Clopper-Pearson confidence interval is not be reliable. Some of the values of alpha, beta and alpha0 may not be plausible given the data.")
    }
  }

  if (is.na(x$ci_asym[1]) || x$ci_asym[1] < x$pi0 || x$ci_asym[2] > 1){
    warning("The asymptotic confidence interval is not be reliable. The point estimate is either too close from the boundary of the space or some of the values of alpha, beta and alpha0 may not be plausible given the data.")
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
#' @param V         A \code{numeric} that corresponds to the average of squared sampling weights. Default value is \code{NULL}.
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
#'  \item beta0:       Estimated false negative rate of the official procedure.
#'  \item ci_beta0:    Asymptotic confidence interval (1 - gamma confidence level) for beta0.
#'  \item boundary:    A boolean variable indicating if the estimates falls at the boundary of the parameter space.
#'  \item pi0:         Value of pi0 (input value).
#'  \item sampling:    Type of sampling considered ("random" or "weighted").
#'  \item V:           Average sum of squared sampling weights if weighted/stratified is used (otherwise NULL).
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
moment_estimator = function(R3, n, pi0, gamma = 0.05, alpha = 0, beta = 0, alpha0 = 0, V = NULL, ...){
  # Check inputs
  if (R3%%1!=0 || R3 < 0){
    if (is.null(V)){
      stop("R3 should be non-negative integer. To use sampling weights you need to provide a value for V.")
    }else{
      if (R3 < 0){
        stop("R3 should be non-negative.")
      }
    }
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
  Delta = 1 - alpha - beta
  estimate = 1/(Delta*(1 - alpha0))*(R3/n + pi0 - beta*pi0 - alpha0*Delta - alpha)

  # Check if the estimate is in the range [pi0, 1]
  if (estimate < pi0 || estimate > 1){
    if (estimate < pi0){
      estimate = pi0
      boundary = "lower"
    }else{
      estimate = 1
      boundary = "upper"
    }
  }else{
    boundary = FALSE
  }

  # Asymptotic CI (1-gamma)
  probs = get_prob(theta = estimate, pi0 = pi0, alpha = alpha,
                   beta = beta, alpha0 = alpha0)

  if (is.null(V)){
    sd = sqrt((probs[3]*(1 - probs[3]))/(n*Delta^2*(1 - alpha0)^2))
    sampling = "random"
  }else{
    sd = sqrt((V*probs[3]*(1 - probs[3]))/(n*Delta^2*(1 - alpha0)^2))
    sampling = "weighted"
  }
  ci_asym = estimate + qnorm(1 - gamma/2)*c(-1,1)*sd

  if (is.null(V)){
    # CP - CI (1 - gamma)
    if (boundary == FALSE){
      upper = (qbeta(p = 1 - gamma/2, R3 + 1, n - R3) - alpha*(1 - alpha0) + (pi0 - alpha0)*(1 - beta))/(Delta*(1 - alpha0))
      lower = (qbeta(p = gamma/2, R3, n - R3 + 1) - alpha*(1 - alpha0) + (pi0 - alpha0)*(1 - beta))/(Delta*(1 - alpha0))
    }else{
      if (boundary == "lower"){
        upper = (1 - (gamma/2)^(1/n) - alpha*(1 - alpha0) + (pi0 - alpha0)*(1 - beta))/(Delta*(1 - alpha0))
        lower = pi0
      }else{
        upper = 1
        lower = ((gamma/2)^(1/n) - alpha*(1 - alpha0) + (pi0 - alpha0)*(1 - beta))/(Delta*(1 - alpha0))
      }
      boundary = TRUE
    }
    ci_cp = c(lower, upper)
  }else{
    ci_cp = NA
  }

  # Compute estimated false negative rate of the official procedure
  beta0 = 1 - (pi0 - alpha0*(1 - estimate))/estimate
  var_beta_asym = (sd^2*n*(pi0 - alpha0)^2)/(estimate^4)
  ci_beta0 = beta0 + qnorm(1 - gamma/2)*c(-1,1)*sqrt(var_beta_asym/n)

  # Construct output
  out = list(estimate = estimate, sd = sd, ci_asym = ci_asym, ci_cp = ci_cp,
             gamma = gamma, method = "Moment Estimator",
             measurement = c(alpha0, alpha, beta), beta0 = beta0,
             ci_beta0 = ci_beta0, boundary = boundary,
             pi0 = pi0, sampling = sampling, V = V, ...)
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
#' @param R1         A \code{numeric} that provides the number of participants in the survey sample that were tested positive with both (medical) testing devices (and are, thus, members of the sub-population).
#' @param R2         A \code{numeric} that provides the number of participants in the survey sample that are tested positive only with the first testing device (and are, thus,  members of the sub-population).
#' @param R3         A \code{numeric} that provides the number of participants in the survey sample that are tested positive only with the second testing device.
#' @param R4         A \code{numeric} that provides the number of participants that are tested negative with the second testing device (and are either members of the sub-population and have tested negative with the first testing device or are not members of the sub-population).
#' @param n          A \code{numeric} that provides the sample size. Default value R1 + R2 + R3 + R4. If this value is provided it is used to verify that R1 + R2 + R3 + R4 = n.
#' @param pi0        A \code{numeric} that provides the prevalence or proportion of people (in the whole population) who are positive, as measured through a non-random,
#' but systematic sampling (e.g. based on medical selection).
#' @param gamma      A \code{numeric} that used to compute a (1 - gamma) confidence region for the proportion. Default value is \code{0.05}.
#' @param alpha0     A \code{numeric} that corresponds to the probability that a random participant
#' has been incorrectly declared positive through the nontransparent procedure. In most applications,
#' this probability is likely very close to zero. Default value is \code{0}.
#' @param alpha      A \code{numeric} that provides the False Negative (FN) rate for the sample R. Default value is \code{0}.
#' @param beta       A \code{numeric} that provides the False Positive (FP) rate for the sample R. Default value is \code{0}.
#' @param V          A \code{numeric} that corresponds to the average of squared sampling weights. Default value is \code{NULL}.
#' @param ...        Additional arguments.
#' @return A \code{cpreval} object with the structure:
#' \itemize{
#'  \item estimate:    Estimated proportion.
#'  \item sd:          Estimated standard error of the estimator.
#'  \item ci_asym:     Asymptotic confidence interval at the 1 - gamma confidence level.
#'  \item gamma:       Confidence level (i.e. 1 - gamma) for confidence intervals.
#'  \item method:      Estimation method (in this case mle).
#'  \item measurement: A vector with (alpha0, alpha, beta).
#'  \item beta0:       Estimated false negative rate of the official procedure.
#'  \item ci_beta0:    Asymptotic confidence interval (1 - gamma confidence level) for beta0.
#'  \item boundary:    A boolean variable indicating if the estimates falls at the boundary of the parameter space.
#'  \item pi0:         Value of pi0 (input value).
#'  \item sampling:    Type of sampling considered ("random" or "weighted").
#'  \item V:           Average sum of squared sampling weights if weighted/stratified is used (otherwise NULL).
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
conditional_mle = function(R1 = NULL, R2 = NULL, R3 = NULL, R4 = NULL, n = R1 + R2 + R3 + R4, pi0, gamma = 0.05, alpha0 = 0, alpha = 0, beta = 0, V = NULL, ...){
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
    if (is.null(V)){
      stop("R1 should be non-negative integer. To use sampling weights you need to provide a value for V.")
    }else{
      if (R1 < 0){
        stop("R1 should be non-negative.")
      }
    }
  }

  if (R2%%1!=0 || R2 < 0){
    if (is.null(V)){
      stop("R2 should be non-negative integer. To use sampling weights you need to provide a value for V.")
    }else{
      if (R2 < 0){
        stop("R2 should be non-negative.")
      }
    }
  }

  if (R3%%1!=0 || R3 < 0){
    if (is.null(V)){
      stop("R3 should be non-negative integer. To use sampling weights you need to provide a value for V.")
    }else{
      if (R3 < 0){
        stop("R3 should be non-negative.")
      }
    }
  }

  if (R4%%1!=0 || R4 < 0){
    if (is.null(V)){
      stop("R4 should be non-negative integer. To use sampling weights you need to provide a value for V.")
    }else{
      if (R4 < 0){
        stop("R4 should be non-negative.")
      }
    }
  }

  # To avoid rounding error
  n = round(n, 6)
  if (n%%1!=0 || n < 0){
    stop("n should be non-negative integer.")
  }

  if (abs(R1 + R2 + R3 + R4 - n) > 10^(-6)){
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

  if (is.null(V)){
    sampling = "random"
    estimate = optimize(neg_log_lik, interval = c(pi0 + eps, 1 - eps), Rvect = R, n = n, pi0 = pi0, alpha = alpha, beta = beta, alpha0 = alpha0)$minimum

    # Check boundary
    LL_mle = (-1)*neg_log_lik(theta = estimate, Rvect = R, n = n, pi0 = pi0, alpha = alpha, beta = beta, alpha0 = alpha0)
    LL_pi0 = (-1)*neg_log_lik(theta = pi0 + 10^(-7), Rvect = R, n = n, pi0 = pi0, alpha = alpha, beta = beta, alpha0 = alpha0)
  }else{
    sampling = "weighted"
    estimate = optimize(neg_log_wlik, interval = c(pi0 + eps, 1 - eps), Rvect = R, n = n, pi0 = pi0, alpha = alpha, beta = beta, alpha0 = alpha0)$minimum

    # Check boundary
    LL_mle = (-1)*neg_log_lik(theta = estimate, Rvect = R, n = n, pi0 = pi0, alpha = alpha, beta = beta, alpha0 = alpha0)
    LL_pi0 = (-1)*neg_log_lik(theta = pi0 + 10^(-7), Rvect = R, n = n, pi0 = pi0, alpha = alpha, beta = beta, alpha0 = alpha0)
  }

  if (LL_mle < LL_pi0){
    boundary = TRUE
    estimate = pi0
  }else{
    boundary = FALSE
  }

  # Compute implied probs
  probs = get_prob(theta = estimate, pi0 = pi0, alpha = alpha, beta = beta, alpha0 = alpha0)

  # Delta
  Delta = 1 - alpha - beta

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
  if (is.null(V)){
    mle_var = 1/(n*fisher_info)
  }else{
    mle_var = V/(n*fisher_info)
  }

  # Compute Asymptotic CI
  sd = sqrt(mle_var)
  ci_asym = estimate + qnorm(1 - gamma/2)*c(-1,1)*sd

  # Compute estimated false negative rate of the official procedure
  beta0 = 1 - (pi0 - alpha0*(1 - estimate))/estimate
  var_beta_asym = (sd^2*n*(pi0 - alpha0)^2)/(estimate^4)
  ci_beta0 = beta0 + qnorm(1 - gamma/2)*c(-1,1)*sqrt(var_beta_asym/n)

  # Construct output
  out = list(estimate = estimate, sd = sd, ci_asym = ci_asym, gamma = gamma,
             method = "Conditional MLE",
             measurement = c(alpha0, alpha, beta),
             beta0 = beta0, ci_beta0 = ci_beta0,
             boundary = boundary, pi0 = pi0, sampling = sampling, V = V,...)
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
#' @param V        A \code{numeric} that corresponds to the average of squared sampling weights. Default value is \code{NULL} and
#' for the moment this method is currently only implemented for random sampling.
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
#'  \item ci_beta0:    Asymptotic confidence interval (1 - gamma confidence level) for beta0.
#'  \item boundary:    A boolean variable indicating if the estimates falls at the boundary of the parameter space.
#'  \item pi0:         Value of pi0 (input value).
#'  \item sampling:    Type of sampling considered ("random" or "weighted").
#'  \item V:           Average sum of squared sampling weights if weighted/stratified is used (otherwise NULL).
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
marginal_mle = function(R1, R3, n, pi0, gamma = 0.05, alpha = 0, beta = 0, alpha0 = 0, V = NULL,...){
  if (!is.null(V)){
    stop("The marginal MLE is currently only implemented for random sampling.")
  }

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
  estimate = optimize(neg_log_lik_integrated, interval = c(pi0 + eps, 1 - eps), Rvect = R, n = n, pi0 = pi0, alpha = alpha, beta = beta, alpha0 = alpha0)$minimum

  # Check boundary
  LL_mle = (-1)*neg_log_lik_integrated(theta = estimate, Rvect = R, n = n, pi0 = pi0, alpha = alpha, beta = beta, alpha0 = alpha0)
  LL_pi0 = (-1)*neg_log_lik_integrated(theta = pi0 + 10^(-7), Rvect = R, n = n, pi0 = pi0, alpha = alpha, beta = beta, alpha0 = alpha0)
  if (LL_mle < LL_pi0){
    boundary = TRUE
    estimate = pi0
  }else{
    boundary = FALSE
  }

  # Compute probs
  probs = get_prob(theta = estimate, pi0 = pi0, alpha = alpha, beta = beta, alpha0 = alpha0)

  # Delta
  Delta = 1 - alpha - beta

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
  var_beta_asym = (sd^2*n*(pi0 - alpha0)^2)/(estimate^4)
  ci_beta0 = beta0 + qnorm(1 - gamma/2)*c(-1,1)*sqrt(var_beta_asym/n)

  # Construct output
  out = list(estimate = estimate, sd = sd, ci_asym = ci_asym, gamma = gamma,
             method = "Marginal MLE", measurement = c(alpha0, alpha, beta),
             beta0 = beta0, ci_beta0 = ci_beta0, boundary = boundary, pi0 = pi0,
             sampling = "random", V = V, ...)
  class(out) = "cpreval"
  out
}

#' @title Negative Weighted Log-Likelihood function
#' @description Log-Likelihood function based on R1, R2, R3 and R4 with sampling weights multiplied by -1.
#' @param theta     A \code{numeric} value for the paramerer of interest.
#' @param Rvect     A \code{vector} of observations, i.e. (R1, R2, R3, R4).
#' @param n         A \code{numeric} that provides the sample size.
#' @param pi0       A \code{numeric} that provides the prevalence or proportion of people (in the whole population) who are positive, as measured through a non-random,
#' but systematic sampling (e.g. based on medical selection).
#' @param alpha     A \code{numeric} that provides the False Negative (FN) rate for the sample R.
#' @param beta      A \code{numeric} that provides the False Positive (FP) rate for the sample R.
#' @param alpha0    A \code{numeric} that corresponds to the probability that a random participant
#' has been incorrectly declared positive through the nontransparent procedure. In most applications,
#' this probability is likely very close to zero.
#' @param ...       Additional arguments.
#' @return Negative log-likelihood.
#' @author Stephane Guerrier
neg_log_wlik = function(theta, Rvect, n, pi0, alpha, beta, alpha0, ...){
  probs = get_prob(theta = theta, pi0 = pi0, alpha = alpha, beta = beta, alpha0 = alpha0)
  (-1)*(Rvect[1]/n*log_modified(probs[1]) + Rvect[2]/n*log_modified(probs[2]) + Rvect[3]/n*log_modified(probs[3]) + Rvect[4]/n*log_modified(probs[4]))
}
