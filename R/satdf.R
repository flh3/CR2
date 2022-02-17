#' Compute Satterthwaite degrees of freedom
#'
#' Function to compute empirical degrees of freedom
#' based on Bell and McCaffrey (2002).
#'
#'
#' @importFrom stats nobs resid formula residuals var coef pt model.matrix family weights fitted.values
#' @param mod The \code{lmerMod} or \code{lme} model object.
#' @param Vinv2 Inverse of the variance matrix.
#' @param Vm2 The variance matrix.
#' @param br2 The bread component.
#' @param Gname The group (clustering variable) name'
#'
#' @export
## empirical DOF
satdf <- function(mod, Vinv2, Vm2, br2, Gname){

  Vinv2 <- as.matrix(Vinv2)
  Vm2 <- as.matrix(Vm2)

  if(class(mod) == 'lme') {
    dat <- mod$data
    fml <- formula(mod)
    X <- model.matrix(fml, data = dat)
    Gname <- names(mod$groups)
    gpsv <- dat[,Gname]

  } else if (class(mod) %in% c('lmerMod', 'lmerModLmerTest')){ #if lmer
    dat <- mod@frame
    X <- model.matrix(mod)

    if (is.null(Gname)){
      Gname <- names(getME(mod, 'l_i')) #name of clustering variable
      if (length(Gname) > 1) {
        stop("lmer: Can only be used with non cross-classified data. If more than two levels, specify highest level using Gname = 'clustername'")
      }
    }

    #print(Gname)
    #Gname <- names(getME(mod, 'l_i'))
    gpsv <- mod@frame[, Gname]

  } else {
    stop("Type of object is not an lmer or lme object.")
  }

  cnames <- names(table(gpsv))
  NG <- length(cnames)
  cdata <- data.frame(cluster = dat[,Gname])
  if(is.unsorted(gpsv)){
    cat("Data are not sorted by cluster. df will be wrong. Please sort your data first by cluster, run the analysis, and then use the function.\n")
  }

  if (NG > 50) cat("Computing Satterthwaite df. There are", NG, "clusters. This may take a while...\n")

  #print(cnames)
  cpx <- chol2inv(chol(t(X) %*% Vinv2 %*% X))

  # cpx <- br2
  Hm <- X %*% cpx %*% t(X) %*% Vinv2

  Id <- diag(nobs(mod))

  Hsel <- Id - Hm #this missing?
  ## STEP 1

  tHs <- function(s) {
    sel <- which(cdata$cluster == s)
    Hsel[,sel]
  }

  tH <- lapply(cnames, tHs) #per cluster

  ## STEP 2

  tXs <- function(s) {
    sel <- which(cdata$cluster == s)
    ss <- diag(length(sel))
    Hg <- Hm[sel, sel]

    #MatSqrtInverse(ss - Hg) %*% Vinv[sel, sel]

    V3 <- chol(Vm2[sel, sel]) #based on MBB
    Bi <- V3 %*% (ss - Hg) %*% Vm2[sel, sel] %*% t(V3)
    t(V3) %*% MatSqrtInverse(Bi) %*% V3


  } # A x Xs / Need this first

  tX <- lapply(cnames, tXs)

  ## step 3

  tFs <- function(s){
    sel <- which(cdata$cluster == s)
    Xs <- X[sel, , drop = FALSE]
    Xs %*% cpx
  }

  tF <- lapply(cnames, tFs)
  ## STEP 3
  k <- ncol(X)
  id <- diag(k) #number of coefficients // for different df
  degf <- numeric(k) #vector for df // container

  for (j in 1:k){ #using a loop since it's easier to see

    Gt <- sapply(seq(NG), function(i) tH[[i]] %*%
                   tX[[i]] %*% tF[[i]] %*% id[,j])
    #already transposed because of sapply: this is G'
    #ev <- eigen(Gt %*% t(Gt))$values #eigen values: n x n
    ev <- eigen(t(Gt) %*% Gt)$values #much quicker this way, same result: p x p:
    degf[j] <- (sum(ev)^2) / sum(ev^2) #final step to compute df
  }

  return(degf)
}
