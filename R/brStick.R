brStick <- function (eigs) {
  p <- length(eigs)
  a <- NULL
  r <- NULL
  
  for(j in 1:p){
    a[j] <- 1/p * sum(1/(j:p) )
    r[j] <- eigs[j]/(sum(eigs))
  }
  length(which(r > a))
}

brStick(r)

