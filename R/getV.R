#' Get V matrix for merMod objects
#'
#' Function to extract V matrix.
#'
#'
#' @param x lme4 object
#' @return V matrix (weight) for multilevel models
#' @export
getV <- function(x) {
  lam <- data.matrix(getME(x, "Lambdat"))
  var.d <- crossprod(lam)
  Zt <- data.matrix(getME(x, "Zt"))
  vr <- sigma(x)^2
  var.b <- vr * (t(Zt) %*% var.d %*% Zt)
  sI <- vr * diag(nobs(x))
  var.y <- var.b + sI
}
