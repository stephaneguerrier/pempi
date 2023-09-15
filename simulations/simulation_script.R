rm(list = ls())
library(pempi)

# Number of Monte-Carlo simulations
B = 5*10^4

# Initial simulation seed
seed = 1832

# Simulation settings
simu = data.frame(Simulation = 1:270,
                  p0 = 100*rep(c(rep(5/100, 30), rep(20/100, 30), rep(75/100, 30)), 3),
                  pi0 = 100*c(rep(c(seq(from = 1/1000, to = 0.9*(5/100), length.out = 30),
                                    seq(from = 1/1000, to = 0.9*(20/100), length.out = 30),
                                    seq(from = 1/1000, to = 0.9*(75/100), length.out = 30)), 2),
                              rep(c(seq(from = 1/100 + 1/1000, to = 0.9*(5/100), length.out = 30),
                                    seq(from = 1/100 + 1/1000, to = 0.9*(20/100), length.out = 30),
                                    seq(from = 1/100 + 1/1000, to = 0.9*(75/100), length.out = 30)), 1)),
                  n = rep(2000, 30*9),
                  alpha0 = 100*rep(0, 3*90),
                  alpha = 100*c(rep(0, 90), rep(0, 90), rep(0.01,90)),
                  beta = 100*c(rep(0, 90), rep(0.02,90), rep(0.02,90)))

m = tail(simu$Simulation, n=1)
error_survey = error_moment = error_mle = error_mmle = matrix(NA, B, m)
length_survey_cp = length_moment_cp = matrix(NA, B, m)
length_survey_asy = length_moment_asy = length_mle_asy = length_mmle_asy = matrix(NA, B, m)
coverage_survey_cp = coverage_moment_cp = matrix(NA, B, m)
coverage_survey_asy = coverage_moment_asy = coverage_mle_asy = coverage_mmle_asy = matrix(NA, B, m)
pb = txtProgressBar(min = 0, max = m, style = 3)

# Start Monte-Carlo
for (j in 1:m){
  for (i in 1:B){
    set.seed(i)
    # Simulate data
    X = sim_Rs(n = simu$n[j],
               theta = simu$p0[j]/100,
               pi0 = simu$pi0[j]/100,
               alpha = simu$alpha[j]/100,
               beta = simu$beta[j]/100,
               alpha0 = simu$alpha0[j]/100)

    # ----------------------------------------------------
    # Fit survey sample estimator
    # ----------------------------------------------------
    survey = survey_mle(R = X$R1 + X$R3, n = X$n,
                           alpha = X$alpha, beta = X$beta,
                           pi0 = X$pi0)

    # Estimation error
    error_survey[i,j] = survey$estimate - simu$p0[j]/100

    # Coverage - CP
    coverage_survey_cp[i,j] = survey$ci_cp[1] < simu$p0[j]/100 && survey$ci_cp[2] > simu$p0[j]/100
    length_survey_cp[i,j] = diff(survey$ci_cp)

    # Coverage - Asym
    coverage_survey_asy[i,j] = survey$ci_asym[1] < simu$p0[j]/100 && survey$ci_asym[2] > simu$p0[j]/100
    if (is.na(sum(survey$ci_asym))){
      length_survey_asy[i,j] = NA
    }else{
      length_survey_asy[i,j] = diff(survey$ci_asym)
    }

    # ----------------------------------------------------
    # Fit moment estimator
    # ----------------------------------------------------
    moment_estim = moment_estimator(R3 = X$R3, n = X$n, pi0 = X$pi0,
                                    alpha0 = X$alpha0, alpha = X$alpha,
                                    beta = X$beta)

    # Estimation error
    error_moment[i,j] = moment_estim$estimate - simu$p0[j]/100

    # Coverage - CP
    coverage_moment_cp[i,j] = moment_estim$ci_cp[1] < simu$p0[j]/100 && moment_estim$ci_cp[2] > simu$p0[j]/100
    length_moment_cp[i,j] = diff(moment_estim$ci_cp)

    # Coverage - Asym
    coverage_moment_asy[i,j] = moment_estim$ci_asym[1] < simu$p0[j]/100 && moment_estim$ci_asym[2] > simu$p0[j]/100
    length_moment_asy[i,j] = diff(moment_estim$ci_asym)


    # ----------------------------------------------------
    # Fit Condition MLE
    # ----------------------------------------------------
    mle_estim = conditional_mle(R1 = X$R1, R2 = X$R2, R3 = X$R3, R4 = X$R4,
                    n = X$n, pi0 = X$pi0, alpha = X$alpha,
                    beta = X$beta, alpha0 = X$alpha0)

    # Estimation error
    error_mle[i,j] = mle_estim$estimate - simu$p0[j]/100

    # Coverage - Asym
    coverage_mle_asy[i,j] = mle_estim$ci_asym[1] < simu$p0[j]/100 && mle_estim$ci_asym[2] > simu$p0[j]/100
    length_mle_asy[i,j] = diff(mle_estim$ci_asym)

    # ----------------------------------------------------
    # Fit Marginal MLE
    # ----------------------------------------------------
    marginal_mle_estim = marginal_mle(R1 = X$R1, R3 = X$R3, n = X$n, pi0 = X$pi0,
                                      alpha = X$alpha, beta = X$beta, alpha0 = X$alpha0)

    # Estimation error
    error_mmle[i,j] = marginal_mle_estim$estimate - simu$p0[j]/100

    # Coverage - Asym
    coverage_mmle_asy[i,j] = marginal_mle_estim$ci_asym[1] < simu$p0[j]/100 && marginal_mle_estim$ci_asym[2] > simu$p0[j]/100
    length_mmle_asy[i,j] = diff(marginal_mle_estim$ci_asym)

  }
  setTxtProgressBar(pb, j)
}
close(pb)

# Save results
coverage = list(mle = coverage_mle_asy,
                mmle = coverage_mmle_asy,
                mom = coverage_moment_asy,
                mom_cp = coverage_moment_cp,
                survey = coverage_survey_asy,
                survey_cp = coverage_survey_cp)

erreur = list(mle = error_mle,
              mmle = error_mmle,
              mom = error_moment,
              survey = error_survey)

len = list(mle = length_mle_asy,
           mmle = length_mmle_asy,
           mom = length_moment_asy,
           mom_cp = length_moment_cp,
           survey = length_survey_asy,
           survey_cp = length_survey_cp)

save(coverage, erreur, len, simu, file = "simulations/simulations.RData")

