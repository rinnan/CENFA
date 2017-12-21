#' @keywords internal
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach %dopar% foreach
#' @importFrom raster blockSize extension filename raster getValues ncell pbClose pbCreate pbStep setValues writeRaster writeStart writeStop writeValues
#' @importFrom snow makeCluster clusterExport stopCluster
#' @importFrom stats cov na.omit
#' @importFrom utils setTxtProgressBar txtProgressBar

.calc <- function(x, fun, filename='', na.rm, forcefun=FALSE, forceapply=FALSE, names, ...) {

  nl <- nlayers(x)

  test <- raster:::.calcTest(x[1:5], fun, na.rm, forcefun, forceapply)
  doapply <- test$doapply
  makemat <- test$makemat
  trans <- test$trans

  if (test$nlout == 1) {
    out <- raster(x)
  } else {
    out <- brick(x, values=FALSE)
    out@data@nlayers <- test$nlout
  }

  names(out) <- names

  fun <- raster:::.makeTextFun(fun)
  if (class(fun) == 'character') {
    doapply <- FALSE
    fun <- raster:::.getRowFun(fun)
  }

  filename <- trim(filename)

  if (canProcessInMemory(x, max(nlayers(x), nlayers(out)) * 2)) {
    x <- getValues(x)
    if (makemat) {
      x <- matrix(x, ncol=1)
    }
    if (missing(na.rm)) {
      if (! doapply ) {
        x <- fun(x )
      } else {
        x <- apply(x, 1, fun )
      }
    } else {
      if ( ! doapply ) {
        x <- fun(x, na.rm=na.rm )
      } else {
        x <- apply(x, 1, fun, na.rm=na.rm)
      }
    }
    if (trans) {
      x <- t(x)
    }
    x <- setValues(out, x)
    if (filename != '') {
      x <- writeRaster(x, filename, ...)
    }
    return(x)
  }

  # else

  out <- writeStart(out, filename=filename, ...)
  tr <- blockSize(out)
  pb <- pbCreate(tr$n, label='calc', ...)

  if (missing(na.rm)) {
    for (i in 1:tr$n) {
      v <- getValues(x, row=tr$row[i], nrows=tr$nrows[i])
      if ( ! doapply ) {
        v <- fun(v)
      } else {
        if (makemat) {
          v <- matrix(v, ncol=1)
        }
        v <- apply(v, 1, fun)
        if (trans) {
          v <- t(v)
        }
      }
      out <- writeValues(out, v, tr$row[i])
      pbStep(pb)
    }
  } else {
    for (i in 1:tr$n) {
      v <- getValues(x, row=tr$row[i], nrows=tr$nrows[i])
      if ( ! doapply ) {
        v <- fun(v, na.rm=na.rm)
      } else {
        if (makemat) {
          v <- matrix(v, ncol=1)
        }
        v <- apply(v, 1, fun, na.rm=na.rm)
        if (trans) {
          v <- t(v)
        }
      }
      out <- writeValues(out, v, tr$row[i])
      pbStep(pb)
    }
  }
  out <- writeStop(out)
  pbClose(pb)
  return(out)
}
