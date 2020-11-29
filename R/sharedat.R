#' Project SHARE
#'
#'Project SHARE (Sexual Health and Relationships) was a cluster randomized trial (CRT) in
#'Scotland carried out to measure the impact of a school-based
#'sexual health program (Wight et al., 2002).
#'
#' @docType data
#'
#' @usage data(sharedat)
#'
#' @format A data frame with 5399 observations and 7 variables.
#'
#' @param school The cluster variable.
#' @param sex factor indicating F or M
#' @param arm treatment arm = 1 vs control = 0
#' @param kscore Pupil knowledge of sexual health
#' @param idno student id number
#' @param sc factor showing the highest social class of the father or mother
#' based on occupation (coded 10: I (highest), 20: II,
#' 31: III non-manual, 32: III manual, 40: IV, 50: V (lowest), 99: not coded).
#' @param zscore standardized knowledge score
#'
#'
#' @keywords datasets
#'
#' @references
#'
#' Moulton, L. (2015). readme.txt contains an overall explanation of the data sets [Data set]. Harvard.
#' (\href{https://doi.org/10.7910/DVN/YXMQZM}{link})
#'
#' Wight, D., Raab, G. M., Henderson, M., Abraham, C., Buston, K., Hart, G., & Scott, S. (2002). Limits of teacher delivered sex education:
#' Interim behavioural outcomes from randomised trial. BMJ, 324, 1430.
#' (\href{https://doi.org/10.1136/bmj.324.7351.1430}{doi})
#'
#'
#' @source \href{https://doi.org/10.7910/DVN/YXMQZM}{Harvard dataverse}
#'
#' @examples
#' data(sharedat)
"sharedat"
