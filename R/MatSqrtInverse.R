MatSqrtInverse <- function(A) {
  ##  Compute the inverse square root of a matrix
  ei <- eigen(A, symmetric = TRUE) #obtain eigenvalues and eigenvectors
  d <- pmax(ei$values, 10^-12) #set negatives values to zero
  #or near zero 10^-12
  d2 <- 1/sqrt(d) #get the inverse of the square root
  d2[d == 0] <- 0
  ei$vectors %*% (if (length(d2)==1) d2 else diag(d2)) %*% t(ei$vectors)
}

