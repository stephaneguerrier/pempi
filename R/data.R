#' @title COVID-19 Data from Statistics Austria
#'
#' @description Data collected in Austria in 2020 (see e.g. SORA, 2020; Kowarik et al., 2021, for more details),
#' allowing to COVID-19 prevalence.
#'
#' @format A \code{matrix} with 2290 rows and 3 variables:
#' \describe{
#'   \item{Y}{Binary variable, 1 if participant i is tested positive in the survey sample, 0 otherwise.}
#'   \item{Z}{Binary variable, 1 if participant i was declared positive with the official procedure, 0 otherwise.}
#'   \item{weights}{Sampling weights.}
#' }
#' @source Statistics Austria. 2020. "Prävalenz von SARS-CoV-2-Infektionen liegt bei 0.031." \url{https://www.statistik.at/web_de/presse/124846.html}
"covid19_austria"
