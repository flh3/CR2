#' Cluster robust standard errors with degrees of freedom adjustments
#'
#' Function to compute the CR0, CR1, CR2 cluster
#' robust standard errors (SE) with Bell and McCaffrey (2002)
#' degrees of freedom (dof) adjustments. Useful when dealing with datasets with a few clusters.
#' Shows output using different CR types and dof choices (for comparative purposes only).
#'
#' @importFrom stats nobs resid residuals var coef pt model.matrix family weights fitted.values
#' @param mod The \code{lm} model object.
#' @param clust The cluster variable (with quotes).
#' @param digits Number of decimal places to display.
#' @param ztest If a normal approximation should be used as the naive degrees of freedom. If FALSE, the between-within degrees of freedom will be used.
#' @return A data frame with the CR adjustments with p-values.
#' \item{estimate}{The regression coefficient.}
#' \item{se.unadj}{The model-based (regular, unadjusted) SE.}
#' \item{CR0}{Cluster robust SE based on Liang & Zeger (1986).}
#' \item{CR1}{Cluster robust SE (using an adjustment based on number of clusters).}
#' \item{CR2}{Cluster robust SE based on Bell and McCaffrey (2002).}
#' \item{tCR2}{t statistic based on CR2.}
#' \item{dfn}{Degrees of freedom(naive): can be infinite (z) or between-within (default). User specified.}
#' \item{dfBM}{Degrees of freedom based on Bell and McCaffrey (2002).}
#' \item{pv.unadj}{p value based on model-based standard errors.}
#' \item{CR0pv}{p value based on CR0 SE with dfBM.}
#' \item{CR0pv.n}{p value  based on CR0 SE with naive df.}
#' \item{CR1pv}{p value based on CR1 SE with dfBM.}
#' \item{CR1pv.n}{p value  based on CR1 SE with naive df.}
#' \item{CR2pv}{p value based on CR2 SE with dfBM.}
#' \item{CR2pv.n}{p value  based on CR2 SE with naive df.}
#'
#' @examples
#' clustSE(lm(mpg ~ am + wt, data = mtcars), 'cyl')
#' data(sch25)
#' clustSE(lm(math ~ ses + minority + mses + mhmwk, data = sch25), 'schid')
#'
#' @references
#' \cite{Bell, R., & McCaffrey, D. (2002). Bias reduction in standard errors for linear regression with multi-stage samples. Survey Methodology, 28, 169-182.
#' (\href{https://www150.statcan.gc.ca/n1/pub/12-001-x/2002002/article/9058-eng.pdf}{link})}
#'
#' Liang, K.Y., & Zeger, S. L. (1986). Longitudinal data analysis using generalized linear models. \emph{Biometrika, 73}(1), 13–22.
#' \doi{10.1093/biomet/73.1.13}
#'
#' @export
clustSE <- function(mod, clust = NULL, digits = 4, ztest = FALSE){

  #if (is.null(data))
  data <- eval(mod$call$data) #OLD
  if (class(mod)[1] != 'glm' & class(mod)[1] != 'lm') stop("Must include a model object of class glm or lm.")
  if (is.null(clust)) stop("Must include a cluster name. clust = 'cluster'")

  tmp <- summary(mod)
  se.n <- tmp$coefficients[,2]
  pv.n <- tmp$coefficients[,4]

  #X <- model.matrix.lm(mod, data, na.action = "na.pass")
  X <- model.matrix(mod) #to keep NAs if
  #if (sum(class(mod) == 'glm')) X <- model.matrix(mod) * sqrt(weights(mod, "working"))
  if (family(mod)[[1]] != 'gaussian') X <- model.matrix(mod) * sqrt(weights(mod, "working"))

  if(nrow(X) != nrow(data)) {
    #warning("Just a note: Missing data in original data.")
    data <- data[names(fitted.values(mod)), ]
  }

  #data[,clust] <- as.character(data[,clust])

  NG <- length(table(data[,clust])) #how many clusters
  #cpx <- solve(crossprod(X)) #(X'X)-1 or inverse of the cp of X
  cpx <- chol2inv(qr.R(qr(X))) #using QR decomposition, faster, more stable?

  cnames <- names(table(data[,clust])) #names of the clusters
  js <- table(data[,clust]) #how many in each cluster
  n <- nobs(mod) #how many total observations
  k <- mod$rank #predictors + intercept

  Xj <- function(x){ #inverse of the symmetric square root (p. 709 IK)
    Xs <- X[data[,clust] == x, , drop = F] #X per cluster
    P <- Xs %*% cpx %*% t(Xs) # the Hat matrix
    return(MatSqrtInverse(diag(nrow(Xs)) - P)) #I - H
  }

  ml <- lapply(cnames, Xj) #need these matrices for CR2 computation

  #########
  #re <- resid(mod) #regular lm
  #if (family(mod)[[1]] != 'gaussian') re <- residuals(mod, "working")  * sqrt(weights(mod, "working")) #manual
  #if (family(mod)[[1]] != 'gaussian') re <- residuals(mod, "pearson")

  re <- resid(mod, 'pearson') #works for both lm and glm
  cdata <- data.frame(data[,clust], re ) #data with cluster and residuals
  names(cdata) <- c('cluster', 'r')
  gs <- names(table(cdata$cluster))

  u1 <- u2 <- matrix(NA, nrow = NG, ncol = k)
  dfa <- (NG / (NG - 1))  * ((n - 1)/(n - k)) #for HC1 / Stata
  #dfa <- NG / (NG - 1) #used by SAS

  ## function for Liang and Zeger SEs
  uu1 <- function(x){
    t(cdata$r[cdata$cluster == x]) %*%
      X[cdata$cluster == x, 1:k] #e'X #plain vanilla
  }

  u1 <- t(sapply(cnames, uu1)) #use as a list?
  #have to transpose to get into proper shape
  #because of sapply
  mt <- crossprod(u1)

  # br <- solve(crossprod(X)) #cpx
  br <- cpx #just copying, got this earlier
  clvc <- br %*% mt %*% br #LZ vcov matrix

  uu2 <- function(x){
    ind <- which(cnames == x)
    t(cdata$r[cdata$cluster == x]) %*% ml[[ind]] %*%
      X[cdata$cluster == x, 1:k]
  }

  u2 <- t(sapply(cnames, uu2)) #have to transpose because of sapply
  mt2 <- crossprod(u2)
  clvc2 <- br %*% mt2 %*% br  #BR LZ2 vcov matrix

  #### To compute empirically-based DF

  ## STEP 1
  tXs <- function(s) {
    Xs <- X[cdata$cluster == s, , drop = F]
    MatSqrtInverse(diag(NROW(Xs)) - Xs %*% cpx %*% t(Xs)) %*%
      Xs
  } # A x Xs / Need this first

  tX <- lapply(cnames, tXs)

  ## STEP 2
  tHs <- function(s) {
    Xs <- X[cdata$cluster == s, , drop = F]
    index <- which(cdata$cluster == s)
    ss <- matrix(0, nrow = n, ncol = length(index)) #all 0, n x G
    ss[cbind(index, 1:length(index))] <- 1 #indicator
    ss - X %*% cpx %*% t(Xs) #overall X x crossprod x Xs'
  }

  tH <- lapply(cnames, tHs) #per cluster

  ## STEP 3

  id <- diag(k) #number of coefficients // for different df
  degf <- numeric(k) #vector for df // container

  for (j in 1:k){ #using a loop since it's easier to see

    Gt <- sapply(seq(NG), function(i) tH[[i]] %*%
                   tX[[i]] %*% cpx %*% id[,j])
    #already transposed because of sapply: this is G'
    #ev <- eigen(Gt %*% t(Gt))$values #eigen values: n x n
    ev <- eigen(t(Gt) %*% Gt)$values #much quicker this way, same result: p x p:
    degf[j] <- (sum(ev)^2) / sum(ev^2) #final step to compute df
  }

  ### Computing the dofHLM

  if (ztest == FALSE){
    ### figuring out Number of L2 and L1 vars for dof

    chk <- function(x){
      vrcheck <- sum(tapply(x, data[,clust], var), na.rm = T) #L1,
      # na needed if only one observation with var = NA
      y <- 1 #assume lev1 by default
      if (vrcheck == 0) (y <- 2) #if variation, then L2
      return(y)
    }

    if (family(mod)[[1]] != 'gaussian') X <- model.matrix(mod) ## use original matrix to check for df
    ## don't need the X matrix after this

    levs <- apply(X, 2, chk) #all except intercept
    levs[1] <- 1 #intercept

    tt <- table(levs)
    l1v <- tt['1']
    l2v <- tt['2']

    l1v[is.na(l1v)] <- 0
    l2v[is.na(l2v)] <- 0

    ####

    #ns <- nobs(mod)
    df1 <- n - l1v - NG #l2v old HLM
    df2 <- NG - l2v - 1

    dfn <- rep(df1, length(levs)) #naive
    dfn[levs == '2'] <- df2
    dfn[1] <- df2 #intercept

      } else {

    dfn <- rep(Inf, k) #infinite
  }

  ### Putting it all together

  CR1 <- sqrt(diag(clvc * dfa))
  CR2 <- sqrt(diag(clvc2))
  CR0 = sqrt(diag(clvc))
  beta <- coef(mod)
  tCR0 <- beta / CR0
  tCR1 <- beta / CR1
  tCR2 <- beta / CR2
  CR0pv.n = round(2 * pt(-abs(tCR0), df = dfn), digits) #naive: either inf or HLM df
  CR0pv = round(2 * pt(-abs(tCR0), df = degf), digits) #using adjusted df
  CR1pv.n = round(2 * pt(-abs(tCR1), df = dfn), digits) #naive: either inf or HLM df
  CR1pv = round(2 * pt(-abs(tCR1), df = degf), digits) #using adjusted df
  CR2pv.n = round(2 * pt(-abs(tCR2), df = dfn), digits) #naive: either inf or HLM df
  CR2pv = round(2 * pt(-abs(tCR2), df = degf), digits) #using adjusted df

  ####
  res <- data.frame(estimate = round(beta, digits),
                    se.unadj = round(se.n, digits),
                    CR0 = round(CR0, digits), #LZeger
                    CR1 = round(CR1, digits), #stata
                    CR2 = round(CR2, digits),
                    tCR2 = round(tCR2, digits),
                    dfn = dfn,
                    dfBM = round(degf, 2),
                    pv.unadj = round(pv.n, digits),
                    CR0pv = round(CR0pv, digits),
                    CR0pv.n = round(CR0pv.n, digits),
                    CR1pv,
                    CR1pv.n,
                    CR2pv.n,
                    CR2pv) #BM
  return(res)
}


## Imbens and Kolesar function :: inverse of symmetric square root
## different from sandwich (if ev < 0, not PD)
MatSqrtInverse <- function(A) {
  ##  Compute the inverse square root of a matrix
  ei <- eigen(A, symmetric = TRUE) #obtain eigenvalues and eigenvectors
  d <- pmax(ei$values, 0) #set negatives values to zero
  d2 <- 1/sqrt(d) #get the inverse of the square root
  d2[d == 0] <- 0
  ei$vectors %*% (if (length(d2)==1) d2 else diag(d2)) %*% t(ei$vectors)
}


