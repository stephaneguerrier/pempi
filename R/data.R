#' @title COVID-19 Survey Data from Statistics Austria (November 2020)
#'
#' @description Survey sample of \eqn{n = 2287} participants tested for
#' COVID-19 in November 2020 by Statistics Austria, used as the case study
#' in Guerrier, Kuzmics & Victoria-Feser (2024).
#'
#' @format A \code{data.frame} with 2287 rows and 3 variables:
#' \describe{
#'   \item{Y}{Binary; 1 if participant \eqn{i} tested positive in the survey sample, 0 otherwise.}
#'   \item{Z}{Binary; 1 if participant \eqn{i} was declared positive by the official procedure, 0 otherwise.}
#'   \item{weights}{Sampling weights.}
#' }
#' @source Statistics Austria. 2020. "Prävalenz von SARS-CoV-2-Infektionen liegt bei 3,1%."
#' @references Guerrier, S., Kuzmics, C., Victoria-Feser, M.-P. (2024). Assessing
#'   COVID-19 Prevalence in Austria with Infection Surveys and Case Count Data
#'   as Auxiliary Information. \emph{Journal of the American Statistical
#'   Association}, 119(547), 1722-1735.
#'   \doi{10.1080/01621459.2024.2313790}
"covid19_austria"
