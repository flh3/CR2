#' Cluster robust standard errors with degrees of freedom adjustments
#'
#' Function to compute the CR2 cluster
#' robust standard errors (SE) with Bell and McCaffrey (2002)
#' degrees of freedom (dof) adjustments.
#'
#' EXPERIMENTAL
#'
#' @importFrom stats nobs resid formula residuals var coef pt model.matrix family weights fitted.values
#' @param m1 The \code{merMod} model object.
#' @param digits Number of decimal places to display.
#' @param satt If Satterthwaite degrees of freedom are to be computed
#' @return A data frame with the CR adjustments with p-values.
#' \item{estimate}{The regression coefficient.}
#' \item{se.unadj}{The model-based (regular, unadjusted) SE.}
#'
#' @references
#' \cite{Bell, R., & McCaffrey, D. (2002). Bias reduction in standard errors for linear regression with multi-stage samples. Survey Methodology, 28, 169-182.
#' (\href{https://www150.statcan.gc.ca/n1/pub/12-001-x/2002002/article/9058-eng.pdf}{link})}
#'
#' Liang, K.Y., & Zeger, S. L. (1986). Longitudinal data analysis using generalized linear models. \emph{Biometrika, 73}(1), 13–22.
#' \doi{10.1093/biomet/73.1.13}
#'
#' @export
cr2_mixed <- function(m1, digits = 4, satt = FALSE){

  ### for lmer
  if(class(m1) %in%  c('lmerMod', 'lmerModLmerTest')){ #if lmer
    dat <- m1@frame
    X <- model.matrix(m1) #X matrix
    B <- fixef(m1) #coefficients
    y <- m1@resp$y #outcome
    Z <- getME(m1, 'Z') #sparse Z matrix
    b <- getME(m1, 'b') #random effects
    Gname <- names(getME(m1, 'l_i')) #name of clustering variable
    if (length(Gname) > 1) {
      stop("lmer: Can only be used with two level data.")
    }
    js <- table(dat[, Gname]) #how many observation in each cluster
    G <- bdiag(VarCorr(m1)) #G matrix

    #re <- as.numeric(y - (X %*% B + Z %*% b)) #not used, just checking
    #data.frame(re, resid(m1)) #the same
    #cor(re, resid(m1)) #1
    # qq <- getME(m1, 'q') #columns in RE matrix

    NG <- getME(m1, 'l_i') #number of groups :: ngrps(m1)
    NG <- NG[length(NG)]

    gpsv <- dat[, Gname] #data with groups

    { #done a bit later than necessary but that is fine
      if(is.unsorted(gpsv)){
       # stop("Data are not sorted by cluster. Please sort your data first by cluster, run the analysis, and then use the function.\n")
      }
    }

    getV <- function(x) {
      lam <- data.matrix(getME(x, "Lambdat"))
      var.d <- crossprod(lam)
      Zt <- data.matrix(getME(x, "Zt"))
      vr <- sigma(x)^2
      var.b <- vr * (t(Zt) %*% var.d %*% Zt)
      sI <- vr * diag(nobs(x))
      var.y <- var.b + sI
    }
    Vm <- getV(m1)
  }

  ## for nlme
  if(class(m1) == 'lme'){ #if nlme
    dat <- m1$data
    fml <- formula(m1)
    X <- model.matrix(fml, data = dat)
    B <- fixef(m1)
    NG <- m1$dims$ngrps[[1]]
    if (length(m1$dims$ngrps) > 3) {stop("Can only be used with two level data.")}
    Gname <- names(m1$groups)
    y <- dat[,as.character(m1$terms[[2]])]
    gpsv <- dat[,Gname]
    js <- table(gpsv)

    { #done a bit later than necessary but that is fine
      if(is.unsorted(gpsv)){
        stop("Data are not sorted by cluster. Please sort your data first by cluster, run the analysis, and then use the function.\n")
      }
    }

    ml <- list()
    for (j in 1:NG){
      test <- getVarCov(m1, individuals = j, type = 'marginal')
      ml[[j]] <- test[[1]]
    }

    Vm <- Matrix::bdiag(ml)
  }

  ### robust computation :: once all elements are extracted
  rr <- y - X %*% B #residuals with no random effects

  cdata <- data.frame(cluster = gpsv, r = rr)
  k <- ncol(X) #number of predictors (inc intercept)
  gs <- names(table(cdata$cluster)) #name of the clusters
  u <- matrix(NA, nrow = NG, ncol = k) #LZ
  uu <- matrix(NA, nrow = NG, ncol = k) #CR2

  #dat <- m1@frame
  cnames <- names(table(gpsv))

  #cpx <- solve(crossprod(X))
  #cpx <- chol2inv(qr.R(qr(X))) #using QR decomposition, faster, more stable?

  #cdata <- data.frame(cluster = dat[,Gname])
  #NG <- length(cnames)


  ### quicker way, doing the bread by cluster
  tmp <- split(X, cdata$cluster)
  XX <- lapply(tmp, function(x) matrix(x, ncol = k)) #X per clust

  # to get Vc
  aa <- function(x){
    sel <- which(cdata$cluster == x)
    chol2inv(chol(Vm[sel, sel]))
    #solve(Vm[sel, sel])
  }

  Vm2 <- lapply(cnames, aa) #Vc used
  names(Vm2) <- cnames #naming

  Vinv <- Matrix::bdiag(Vm2)
  # to get X V-1 X per cluster
  bb <- function(x){
    t(XX[[x]]) %*% Vm2[[x]] %*% XX[[x]]
  }

  dd <- lapply(cnames, bb)
  br <- solve(Reduce("+", dd)) #bread

  #Vinv <- solve(Vm) #slow
  #Vinv <- chol2inv(chol(Vm))

  #Q <- solve(t(X) %*% Vinv %*% X)
  #print(Q)
  ## mtsqrtinv is used below...


  tXs <- function(s) {
    # ind <- which(gpsv == s)
    # Xs <- X[ind, , drop = F]
    #
    # V2 <- solve(Vm[ind, ind])
    # #rr <- u[ind]
    # Ijj <- diag(length(ind))
    # Hjj <- Xs %*% Q %*% t(Xs) %*% V2

    ### new

    Ijj <- diag(nrow(XX[[s]]))
    Hjj <- XX[[s]] %*% br %*% t(XX[[s]]) %*% Vm2[[s]]

    MatSqrtInverse(Ijj - Hjj)

  } # A x Xs / Need this first

  tX <- lapply(cnames, tXs)

  # for(i in 1:NG){
  #   #print(i)
  #   #tmp <- js[i] #how many in group
  #   #LZ
  #   # u[i,] <- as.numeric(t(cdata$r[cdata$cluster == gs[i]]) %*% solve(ml[[i]]) %*% X[gpsv == gs[i], 1:k])
  #   #CR2
  #   # uu[i,] <- as.numeric(t(cdata$r[cdata$cluster == gs[i]]) %*% solve(ml[[i]]) %*% tX[[i]])
  #
  #   ind <- gpsv == gs[i]
  #   Xs <- X[ind, , drop = F]
  #   #V2 <- solve(Vm[ind, ind])
  #   #V2 <- chol2inv(chol(Vm[ind, ind])) #v^-1
  #   #CR0
  #   u[i,] <- as.numeric(t(cdata$r[ind]) %*% Vm2[[i]] %*% Xs) #V2
  #   #CR2
  #   uu[i,] <- as.numeric(t(cdata$r[ind]) %*% tX[[i]] %*% Vm2[[i]]  %*% XX[[i]]) #Xs
  #
  # }

  rrr <- split(rr, getME(m1, 'flist'))

  cc0 <- function(x){
    rrr[[x]] %*% Vm2[[x]] %*% XX[[x]]
  }

  u <- t(sapply(cnames, cc0))

  cc2 <- function(x){
    rrr[[x]] %*% tX[[x]] %*% Vm2[[x]] %*% XX[[x]]
  }

  uu <- t(sapply(1:NG, cc2)) #using 1:NG instead


  ## e'(Vg)-1 Xg ## CR0
  ## putting the pieces together

  #br2 <- solve(t(X) %*% Vinv %*% X) #bread
  mt <- t(u) %*% u #meat :: t(u) %*% u
  clvc2 <- br %*% mt %*% br
  rse <- sqrt(diag(clvc2))

  mt2 <- t(uu) %*% uu #meat :: t(u) %*% u
  clvc2a <- br %*% mt2 %*% br
  rse2 <- sqrt(diag(clvc2a))

  #print('done')

  ### HLM dof
  chk <- function(x){
    vrcheck <- sum(tapply(x, gpsv, var), na.rm = T) #L1,
    # na needed if only one observation with var = NA
    y <- 1 #assume lev1 by default
    if (vrcheck == 0) (y <- 2) #if variation, then L2
    return(y)
  }

  levs <- apply(X, 2, chk) #all except intercept
  # levs[1] <- 1 #intercept

  tt <- table(levs)
  l1v <- tt['1']
  l2v <- tt['2']

  l1v[is.na(l1v)] <- 0
  l2v[is.na(l2v)] <- 0

  ####
  n <- nobs(m1)
  #ns <- nobs(mod)
  df1 <- n - l1v - length(js)
  df2 <- NG - l2v

  dfn <- rep(df1, length(levs)) #naive
  dfn[levs == '2'] <- df2

  dfn.CR0 <- dfn

  if (satt == T){
    dfn <- satdf(m1)
  }

  robse <- as.numeric(rse)
  FE_auto <- fixef(m1)
  statistic.cr0 <- FE_auto / robse
  p.values.cr0 = round(2 * pt(-abs(statistic.cr0), df = dfn.CR0), digits) #using CR0

  stars.cr0 <- cut(p.values.cr0, breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   labels = c("***", "**", "*", ".", " "), include.lowest = TRUE)

  statistic.cr2 <- FE_auto / rse2
  p.values.cr2 = round(2 * pt(-abs(statistic.cr2), df = dfn), digits)
  ### changed for sim-- using HLM df for both
  ### change to df = dfn for satt dfn
  stars.cr2 <- cut(p.values.cr2, breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   labels = c("***", "**", "*", ".", " "), include.lowest = TRUE)

  ################# COMPARE RESULTS

  #gams <- solve(t(X) %*% solve(Vm) %*% X) %*% (t(X) %*% solve(Vm) %*% y)
  #SEm <- as.numeric(sqrt(diag(solve(t(X) %*% solve(Vm) %*% X)))) #X' Vm-1 X
  #SE <- as.numeric(sqrt(diag(vcov(m1)))) #compare standard errors
  SE <- as.numeric(sqrt(diag(br)))
  return(data.frame(
    #FE_manual = as.numeric(gams),
    FE_auto,
    #SE_manual = SEm,
    SE_auto = SE,
    dof = dfn,
    cr0 = robse,
    cr2 = rse2,
    # this is just to eval
    # robust_auto = sqrt(diag(clubSandwich::vcovCR(m1, cluster = cdata$cluster, type = 'CR2' ))),
    p.values.cr0,
    p.values.cr2,
    stars.cr0,
    stars.cr2
  )
  )

}

MatSqrtInverse <- function(A) {
  ##  Compute the inverse square root of a matrix
  ei <- eigen(A, symmetric = TRUE) #obtain eigenvalues and eigenvectors
  d <- pmax(ei$values, 10^-12) #set negatives values to zero
  #or near zero 10^-12
  d2 <- 1/sqrt(d) #get the inverse of the square root
  d2[d == 0] <- 0
  ei$vectors %*% (if (length(d2)==1) d2 else diag(d2)) %*% t(ei$vectors)
}


## empirical DOF

satdf <- function(mod, Vinv = Vinv){
  if(class(mod) == 'lme') {
    dat <- mod$data
    fml <- formula(mod)
    X <- model.matrix(fml, data = dat)
    Gname <- names(mod$groups)
    gpsv <- dat[,Gname]

  } else if (class(mod) %in% c('lmerMod', 'lmerModLmerTest')){ #if lmer
    dat <- mod@frame
    X <- model.matrix(mod)
    Gname <- names(getME(mod, 'l_i'))
    gpsv <- mod@frame[, Gname]

    getV <- function(x) {
      lam <- data.matrix(getME(x, "Lambdat"))
      var.d <- crossprod(lam)
      Zt <- data.matrix(getME(x, "Zt"))
      vr <- sigma(x)^2
      var.b <- vr * (t(Zt) %*% var.d %*% Zt)
      sI <- vr * diag(nobs(x))
      var.y <- var.b + sI
    }
    Vm <- getV(mod)
    #Vinv <- solve(Vm)
    Vinv <- chol2inv(chol(Vm))

  } else {
    stop("Type of object is not an lmer or lme object.")
  }

  cnames <- names(table(gpsv))
  #cpx <- solve(crossprod(X))
  cpx <- solve(t(X) %*% Vinv %*% X) #doing it the old way
  Hm <- X %*% cpx %*% t(X) %*% Vinv
  cdata <- data.frame(cluster = dat[,Gname])
  NG <- length(cnames)
  Id <- diag(nobs(mod))

  Hsel <- Id - Hm #this missing?
  ## STEP 1

  tHs <- function(s) {
    sel <- which(cdata$cluster == s)

    #Xs <- X[index, , drop = FALSE]
    #Xs <- X[cdata$cluster == s, , drop = F]

    #ss <- matrix(0, nrow = n, ncol = length(index)) #all 0, n x G
    #ss[cbind(index, 1:length(index))] <- 1 #indicator

    #ss - Hm[sel, sel] #overall X x crossprod x Xs'
    Hsel[,sel]
  }

  tH <- lapply(cnames, tHs) #per cluster


  ## STEP 2

  tXs <- function(s) {
    sel <- which(cdata$cluster == s)

    #Xs <- X[sel, , drop = FALSE]
    #Xs <- X[cdata$cluster == s, , drop = F]
    ss <- diag(length(sel))

    Hg <- Hm[sel, sel]

    MatSqrtInverse(ss - Hg) %*% Vinv[sel, sel] #solve(Vm[sel, sel])

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



