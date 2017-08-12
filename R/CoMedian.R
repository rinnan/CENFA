is.1num<-function (x) is.numeric(x) && length(x) == 1L

CorMed<-function (X) {   #correlation median equiv to Pearson's coefficient of correlation:
  stopifnot(is.1num(p <- ncol(X)), p >= 1)         #COM(X,Y)/(MAD(X)*MAD(Y))
  med <- colMedians(X)
  Y <- sweep(X, 2L, med, `-`)
  CM <- matrix(0, p, p)
  madY <- numeric(p)
  for (i in 1:p) {
    madY[i] <- madYi <- mad(Yi <- Y[, i])
    for (j in seq_len(i - 1)) {
      CM[j, i] <- CM[i, j] <- median(Yi * Y[, j])/(madYi * madY[j])
    }
    CM[i, i] <- median(Yi^2)/(madYi^2)
  }
  1.4826^2 * CM
}

CoMedian<-function (X, na.rm=TRUE) {
  stopifnot(is.1num(p <- ncol(X)), p >= 1)
  med <- colMedians(X, na.rm=na.rm)
  Y <- sweep(X, 2L, med, `-`)
  CM <- matrix(0, p, p)
  madY <- numeric(p)
  for (i in 1:p) {
    madY[i] <- madYi <- mad(Yi <- Y[, i])
    for (j in seq_len(i - 1)) {
      CM[j, i] <- CM[i, j] <- median(Yi * Y[, j])
    }
    CM[i, i] <- median(Yi^2)
  }
  1.4826^2 *CM
}

