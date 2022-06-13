#' Simulated data from 18 schools (from a cluster randomized controlled trial)
#'
#' Synthetic dataset used in the manuscript in the Journal of Research on Educational Effectiveness.
#'
#' @usage data(crct)
#' @format A data frame with 4233 rows and 12 variables:
#' \describe{
#'   \item{usid}{Unique school identifier (the grouping variable).}
#'   \item{stype}{School type (elementary, middle, or high school).}
#'   \item{trt}{Treatment indicator. 1 = intervention; 0 = control.}
#'   \item{odr_post}{Office disciplinary referral outcome.}
#'   \item{odr_pre}{Office disciplinary referral (baseline).}
#'   \item{size}{School enrollment size (to the nearest hundred).}
#'   \item{female}{Student is female: 1 = yes.}
#'   \item{stype_ms}{Dummy code for school type; middle school.}
#'   \item{stype_elem}{Dummy code for school type; elementary school.}
#'   \item{stype_high}{Dummy code for school type; high school.}
#'   \item{race_Black}{Dummy code for student race/ethnicity; Black student.}
#'   \item{race_Hispanic}{Dummy code for student race/ethnicity; Hispanic student.}
#' }
"crct"
