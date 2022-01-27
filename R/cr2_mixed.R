## FH robust_mixed function
## 2021.11.21

cr2_mixed <- function(m1, digits = 4, satt = F){
  require(Matrix)

  ### for lmer
  if(class(m1) %in%  c('lmerMod', 'lmerModLmerTest')){ #if lmer
    X <- model.matrix(m1) #X matrix
    B <- fixef(m1) #coefficients
    y <- m1@resp$y #outcome
    Z <- getME(m1, 'Z') #sparse Z matrix
    b <- getME(m1, 'b') #random effects
    Gname <- names(getME(m1, 'l_i')) #name of clustering variable
    if (length(Gname) > 1) {
      stop("lmer: Can only be used with two level data.")
    }
    js <- table(m1@frame[, Gname]) #how many observation in each cluster
    G <- bdiag(VarCorr(m1)) #G matrix

    #re <- as.numeric(y - (X %*% B + Z %*% b)) #not used, just checking
    #data.frame(re, resid(m1)) #the same
    #cor(re, resid(m1)) #1
    # qq <- getME(m1, 'q') #columns in RE matrix

    NG <- getME(m1, 'l_i') #number of groups :: ngrps(m1)
    NG <- NG[length(NG)]

    gpsv <- m1@frame[, Gname] #data with groups

    { #done a bit later than necessary but that is fine
      if(is.unsorted(gpsv)){
        stop("Data are not sorted by cluster. Please sort your data first by cluster, run the analysis, and then use the function.\n")
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

    Vm <- bdiag(ml)
  }

  ### robust computation :: once all elements are extracted
  rr <- y - X %*% B #residuals with no random effects

  cdata <- data.frame(cluster = gpsv, r = rr)
  k <- ncol(X) #
  gs <- names(table(cdata$cluster)) #name of the clusters
  u <- matrix(NA, nrow = NG, ncol = k) #LZ
  uu <- matrix(NA, nrow = NG, ncol = k) #CR2

  #dat <- m1@frame
  cnames <- names(table(gpsv))
  #cpx <- solve(crossprod(X))
  #cpx <- chol2inv(qr.R(qr(X))) #using QR decomposition, faster, more stable?

  #cdata <- data.frame(cluster = dat[,Gname])
  #NG <- length(cnames)
  Vinv <- solve(Vm) #slow

  Q <- solve(t(X) %*% Vinv %*% X)

  ## mtsqrtinv is used below...
  tXs <- function(s) {
    ind <- which(gpsv == s)
    Xs <- X[ind, , drop = F]

    #MatSqrtInverse(diag(NROW(Xs)) - Xs %*% cpx %*% t(Xs)) #%*%
    # Xs
    #MatSqrtInverse(diag(NROW(Xs)) - Xs %*% cpx %*% t(Xs) #%*% Vinv[ind, ind])
    V2 <- solve(Vm[ind, ind])
    #rr <- u[ind]
    Ijj <- diag(length(ind))
    Hjj <- Xs %*% Q %*% t(Xs) %*% V2
    MatSqrtInverse(Ijj - Hjj)

    #MatSqrtInverse(diag(length(ind)) - HAT[ind, ind]) #%*% Vinv[ind, ind])

  } # A x Xs / Need this first

  tX <- lapply(cnames, tXs)

  for(i in 1:NG){
    #tmp <- js[i] #how many in group
    #LZ
    # u[i,] <- as.numeric(t(cdata$r[cdata$cluster == gs[i]]) %*% solve(ml[[i]]) %*% X[gpsv == gs[i], 1:k])
    #CR2
    # uu[i,] <- as.numeric(t(cdata$r[cdata$cluster == gs[i]]) %*% solve(ml[[i]]) %*% tX[[i]])

    ind <- gpsv == gs[i]
    Xs <- X[ind, , drop = F]
    V2 <- solve(Vm[ind, ind])
    #CR0
    u[i,] <- as.numeric(t(cdata$r[ind]) %*% V2 %*% Xs)
    #CR2
    uu[i,] <- as.numeric(t(cdata$r[ind]) %*% tX[[i]] %*% V2  %*% Xs)

  }

  ## e'(Vg)-1 Xg ## CR0
  ## putting the pieces together

  br2 <- solve(t(X) %*% Vinv %*% X) #bread
  mt <- t(u) %*% u #meat :: t(u) %*% u
  clvc2 <- br2 %*% mt %*% br2
  rse <- sqrt(diag(clvc2))

  mt2 <- t(uu) %*% uu #meat :: t(u) %*% u
  clvc2a <- br2 %*% mt2 %*% br2
  rse2 <- sqrt(diag(clvc2a))

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

  gams <- solve(t(X) %*% solve(Vm) %*% X) %*% (t(X) %*% solve(Vm) %*% y)
  SEm <- as.numeric(sqrt(diag(solve(t(X) %*% solve(Vm) %*% X)))) #X' Vm-1 X
  SE <- as.numeric(sqrt(diag(vcov(m1)))) #compare standard errors
  return(data.frame(
    FE_manual = as.numeric(gams),
    FE_auto,
    SE_manual = SEm,
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
  ei <- eigen(A) #obtain eigenvalues and eigenvectors
  d <- pmax(ei$values, 10^-12) #set negatives values to zero
  #or near zero 10^-12
  d2 <- 1/sqrt(d) #get the inverse of the square root
  d2[d == 0] <- 0
  ei$vectors %*% (if (length(d2)==1) d2 else diag(d2)) %*% t(ei$vectors)
}


## empirical DOF

satdf <- function(mod){
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
    Vinv <- solve(Vm)

  } else {
    stop("Type of object is not an lmer or lme object.")
  }

  cnames <- names(table(gpsv))
  #cpx <- solve(crossprod(X))
  cpx <- solve(t(X) %*% Vinv %*% X)
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

    MatSqrtInverse(ss - Hg)

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


