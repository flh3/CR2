#' Testing for nonconstant variance (ncv)
#'
#' @description{Function to detect heteroscedasticity in two-level random intercept models.
#' Uses a generalization of the Breusch-Pagan-type (using squared residuals)
#' and Levene-type test (using the absolute value of residuals). Note: this will
#' not tell you if including random slopes are warranted (for that, use the
#' \code{robust_mixed}) function and compare differences in model-based and
#' robust standard errors.}
#'
#'
#'
#' @importFrom stats resid anova update
#' @importFrom nlme ranef
#' @importFrom lme4 lmer
#' @param mx The \code{lme} or \code{merMod} model object.
#' @param bp Computes a Breusch-Pagan-type test (\code{TRUE}). If \code{FALSE} computes a Levene-type test.
#' @return A p-value (p < .05 suggests heteroskedasticity).
#'
#' @references
#' \cite{Huang, F., Wiedermann, W., & Zhang, B. (2022). Accounting for Heteroskedasticity Resulting from Between-group Differences in Multilevel Models. Multivariate Behavioral Research.
#' }
#'
#'
#' @examples
#' require(lme4)
#' data(sch25)
#' ncvMLM(lmer(math ~ byhomewk + male + ses + (1|schid), data = sch25)) #supported
#' ncvMLM(lmer(math ~ byhomewk + male + ses + minority + (1|schid), data = sch25)) #hetero
#' @export ncvMLM
ncvMLM <- function(mx, bp = TRUE){
  #if (class(mx) != 'lme') (stop("Only for lme objects"))
  if(is(mx, 'lme')){
  dat <- mx$data
  if (length(mx$dims$ngrps) > 3) {stop("Can only be used with two level data.")}
  if (ncol(ranef(mx)) > 1) {stop("Can only be used with random intercept models.")}

  if (bp == TRUE){
    dat$rr <- resid(mx)^2 #squared, these are conditional
  } else { #this is Levene's
    dat$rr <- abs(resid(mx)) #absolute value, these are conditional
  }
  tmp0 <- update(mx, rr ~ 1, method = 'ML', data = dat)
  tmp1 <- update(mx, rr ~ ., method = 'ML', data = dat)

  res <- anova(tmp1, tmp0)
  return(as.numeric(res$`p-value`[2]))
  } else if(is(mx, 'merMod')){

    dat <- mx@frame

    Gname <- names(getME(mx, 'l_i')) #name of clustering variable

    if (length(Gname) > 1) {
      stop("lmer: Can only be used with non cross-classified data.")
    }
    if (getME(mx, 'p_i') > 1) {
      stop("Can only be used with random intercept models.")
    }

    if (bp == TRUE){
      dat$rr <- resid(mx)^2 #squared, these are conditional
    } else { #this is Levene's
      dat$rr <- abs(resid(mx)) #absolute value, these are conditional
    }

    groups <- dat[,Gname]
    #print(groups)
    tmp0 <- lmer(rr ~ 1 + (1|groups), REML = FALSE, data = dat)
    tmp1 <- update(mx, rr ~ ., REML = FALSE, data = dat)

    res <- anova(tmp1, tmp0)
    return(as.numeric(res$`Pr(>Chisq)`[2]))

  } else {
    stop("Can only be used with lme or merMod objects.")
  }
}
