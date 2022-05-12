#' Grade point average (GPA) data of students from 25 schools
#'
#' For investigating heteroskedasticity.
#' @usage data(gpadat)
#' @format A data frame with 8,956 rows and 18 variables:
#' \describe{
#'   \item{gpa}{Grade point average. 1 = D ... 4 = A.}
#'   \item{female}{Gender. Female = 1.}
#'   \item{race}{Student race/ethnicity (factor).}
#'   \item{dis}{Disability status (1 = yes/0 = no).}
#'   \item{frpl}{Free/reduced price lunch status.}
#'   \item{race_w}{Dummy coded race (White).}
#'   \item{race_a}{Dummy coded race (Asian).}
#'   \item{race_b}{Dummy coded race (Black).}
#'   \item{race_h}{Dummy coded race (Hispanic).}
#'   \item{race_o}{Dummy coded race (Other).}
#'   \item{per_asian}{Group-aggregated Asian variable.}
#'   \item{per_black}{Group-aggregated Black variable.}
#'   \item{per_hisp}{Group-aggregated Hispanic variable.}
#'   \item{per_other}{Group-aggregated Other variable.}
#'   \item{per_fem}{Group-aggregated female variable.}
#'   \item{per_dis}{Group-aggregated disability variable.}
#'   \item{per_frpl}{Group-aggregated frpl variable.}
#'   \item{schoolid}{School identifier (cluster variable).}
#'   }
"gpadat"
