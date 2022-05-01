#' @title Compute sucess probabilities (tau_j's)
#' @description Compute joint probabilities of P(W = j, Y = k) for j, k = 0, 1.
#' @param theta        A \code{numeric} that provides the true prevalence of a given disease.
#' @param pi0          A \code{numeric} that provides the prevalence or proportion of people (in the whole population) who are positive, as measured through a non-random,
#' but systematic sampling (e.g. based on medical selection).
#' @param alpha        A \code{numeric} that provides the False Negative (FN) rate for the sample R.
#' @param beta         A \code{numeric} that provides the False Positive (FP) rate for the sample R.
#' @param alpha0       A \code{numeric} that corresponds to the probability that a random participant
#' has been incorrectly declared positive through the nontransparent procedure. In most applications,
#' this probability is likely very close to zero.
#' @return A \code{vector} containing tau1, tau2, tau3 and tau4.
#' @export
#' @author Stephane Guerrier
#' @examples
#' prob1 = get_prob(theta = 0.02, pi0 = 0.01, alpha = 0, beta = 0, alpha0 = 0)
#' prob1
#' sum(prob1)
#'
#' prob2 = get_prob(theta = 0.02, pi0 = 0.01, alpha = 0.001, beta = 0, alpha0 = 0.001)
#' prob2
#' sum(prob2)
get_prob = function(theta, pi0, alpha, beta, alpha0){
  beta0 = 1 - (pi0 - alpha0*(1 - theta))/theta
  Delta = 1 - alpha - beta
  Delta0 = 1 - alpha0 - beta0
  pi0_star = (pi0 - alpha0)/Delta0

  tau1 = theta*Delta*alpha0 + pi0_star*Delta0*(1 - beta) + alpha*alpha0
  tau2 = -theta*Delta*alpha0 + pi0_star*Delta0*beta + (1 - alpha)*alpha0
  tau3 = theta*Delta*(1 - alpha0) - pi0_star*Delta0*(1 - beta) + alpha*(1 - alpha0)
  tau4 = -theta*Delta*(1 - alpha0) - pi0_star*Delta0*beta + (1 - alpha)*(1 - alpha0)

  c(tau1, tau2, tau3, tau4)
}


#' @title Simulate data (R, R0, R1, R2, R3 and R4)
#' @description Simulation function for random variables of interest.
#' @param theta        A \code{numeric} that provides the true prevalence of a given disease.
#' @param pi0          A \code{numeric} that provides the prevalence or proportion of people (in the whole population) who are positive, as measured through a non-random,
#' but systematic sampling (e.g. based on medical selection).
#' @param n            A \code{numeric} that corresponds to the sample size.
#' @param alpha        A \code{numeric} that provides the False Negative (FN) rate for the sample R. Default value is \code{0}.
#' @param beta         A \code{numeric} that provides the False Positive (FP) rate for the sample R. Default value is \code{0}.
#' @param alpha0       A \code{numeric} that corresponds to the probability that a random participant
#' has been incorrectly declared positive through the nontransparent procedure. In most applications,
#' this probability is likely very close to zero. Default value is \code{0}.
#' @param seed         A \code{numeric} that provides the simulation seed. Default value is \code{NULL}.
#' @param ...          Additional arguments.
#' @return A \code{cpreval_sim} object (\code{list}) with the structure:
#' \itemize{
#' \item R:      the number of participants in the survey sample that were tested positive.
#' \item R0:     the number of participants in the survey sample that were tested positive with the first testing device (and are, thus,  members of the sub-population).
#' \item R1:     the number of participants in the survey sample that were tested positive with both (medical) testing devices (and are, thus, members of the sub-population).
#' \item R2:     the number of participants in the survey sample that are tested positive only with the first testing device (and are, thus,  members of the sub-population).
#' \item R3:     the number of participants in the survey sample that are tested positive only with the second testing device.
#' \item R4:     the number of participants that are tested negative with the second testing device (and are either members of the sub-population and have tested negative with the first testing device or are not members of the sub-population).
#' \item n:      the sample size.
#' \item alpha:  the False Negative (FN) rate for the sample R.
#' \item beta:   the False Positive (FP) rate for the sample R.
#' \item alpha0: the alpha0 probability (as defined above).
#' \item ...:    additional arguments.
#' }
#' @export
#' @author Stephane Guerrier
#' @examples
#' # Samples without measurement error
#' sim_Rs(theta = 3/100, pi0 = 1/100, n = 1500, seed = 18)
#'
#' # With measurement error
#' sim_Rs(theta = 3/100, pi0 = 1/100, n = 1500, alpha0 = 0,
#' alpha = 0.01, beta = 0.05, seed = 18)
#' @importFrom stats rmultinom
sim_Rs = function(theta, pi0, n, alpha0 = 0, alpha = 0, beta = 0, seed = NULL, ...){
  # Simulation seed
  if (!is.null(seed)){
    set.seed(seed)
  }

  # Check inputs
  if (n %% 1 != 0){
    stop("The input n must be an integer.")
  }

  if (alpha0 > pi0){
    stop("The inputs pi0 and alpha0 must be such that alpha0 <= pi0.")
  }

  if (min(theta, pi0) <= 0 || max(theta, pi0) >= 1 || pi0 >= theta){
    stop("The inputs pi0 and theta must be such that 0 < pi0 < theta < 1.")
  }

  if (alpha0 < 0 || alpha0 > 1){
    stop("The input alpha0 must be such that 0 <= alpha0 < 1.")
  }

  if (alpha < 0 || alpha > 1){
    stop("The input alpha must be such that 0 <= alpha < 1.")
  }

  beta0 = 1 - (pi0 - alpha0*(1 - theta))/theta
  if (beta0 < 0 || beta0 > 1){
    stop("The input beta0 must be such that 0 <= beta0 < 1.")
  }

  if (beta < 0 || beta > 1){
    stop("The input beta must be such that 0 <= beta < 1.")
  }

  if (alpha + beta >= 1){
    stop("The inputs alpha and beta must be alpha + beta < 1 (see Assumption 1).")
  }

  if (alpha0 + beta0 >= 1){
    stop("The inputs alpha0 and beta0 must be alpha0 + beta0 < 1 (see Assumption 1).")
  }

  # Compute sucess probabilities
  probs = get_prob(theta = theta, pi0 = pi0, alpha = alpha, beta = beta, alpha0 = alpha0)

  # Draw data from a multinomial distribution
  X = rmultinom(n = 1, size = n, prob = probs)

  # Output list
  structure(list(R = X[1,1] + X[3,1], R0 = X[1,1] + X[2,1],
       R1 = X[1,1], R2 = X[2,1], R3 = X[3,1],
       R4 = X[4,1], n = n, pi0 = pi0,
       alpha = alpha, beta = beta,
       alpha0 = alpha0, beta0 = beta0, ...), class = "cpreval_sim")
}

#' @title Print (simulated) sample
#' @description Simple print function for the \code{cpreval_sim} objects
#' @method print cpreval_sim
#' @export
#' @keywords internal
#' @param x    A \code{cpreval_sim} object
#' @param ...  Further arguments passed to or from other methods
#' @return Prints object
#' @author Stephane Guerrier
#' @examples
#' # Samples without measurement error
#' sim_Rs(theta = 3/100, pi0 = 1/100, n = 1500, seed = 18)
#'
#' # With measurement error
#' sim_Rs(theta = 3/100, pi0 = 1/100, n = 1500, alpha0 = 0,
#' alpha = 0.01, beta = 0.05, seed = 18)
print.cpreval_sim = function(x, ...){
  cat("Data: R = ")
  cat(x$R)
  cat(", R0 = ")
  cat(x$R0)
  cat(", n = ")
  cat(x$n)
  cat("\n      R1 = ")
  cat(x$R1)
  cat(", R2 = ")
  cat(x$R2)
  cat(", R3 = ")
  cat(x$R3)
  cat(", R4 = ")
  cat(x$R4)
  cat("\n\n")

  cat("Assumed measurement error: alpha = ")
  cat(100*x$alpha)
  cat("%, alpha0 = ")
  cat(100*x$alpha0)
  cat("%, beta  = ")
  cat(100*x$beta)
  cat("% \n\n")

  cat("False negative rate of the official procedure: beta0 = ")
  cat(round(100*x$beta0,2))
  cat("%\n")
}


