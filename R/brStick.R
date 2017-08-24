brStick <- function (eigs) {
  if(max(Im(eigs)) > 1e-5) stop("broken-stick method does not work for complex eigenvalues")
  eigs <- Re(eigs)
  p <- length(eigs)
  a <- NULL
  r <- NULL

  for(j in 1:p){
    a[j] <- 1/p * sum(1/(j:p) )
    r[j] <- eigs[j]/(sum(eigs))
  }
  length(which(r > a))
}

