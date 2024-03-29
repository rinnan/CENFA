#' @keywords internal
#'
# @importFrom raster brick trim canProcessInMemory setValues writeRaster pbStep
#   nlayers getValues writeStart blockSize pbCreate writeValues writeStop pbClose
#   extract ncell
#' @importFrom foreach %dopar% foreach
#' @importFrom stats cov na.omit
#' @importFrom utils setTxtProgressBar txtProgressBar

# This function is the exact same as raster::overlay, with the addition of a
# 'names' argument, which allows the user to specify a character vector of
# layer names for the Raster* object that gets created.
#
# Internal functions from the raster package are included here, as
# necessary for raster::overlay to run.

.overlay <- function(x, y, ..., fun, filename="", recycle=TRUE, forcefun=FALSE, names){
  if (missing(fun)) {
    stop("you must supply a function 'fun'.\nE.g., 'fun=function(x,y){return(x+y)} or fun=sum'")
  }
  lst <- list(...)
  isRast <- sapply(lst, function(x) inherits(x, 'Raster'))
  if (sum(unlist(isRast)) > 0) {
    x <- c(x, y, lst[isRast])
    lst <- lst[! isRast ]
  } else {
    x <- list(x, y)
  }
  lst$fun <- fun
  lst$filename <- filename
  lst$recycle <- recycle
  lst$forcefun <- forcefun
  lst$x <- x
  lst$names <- names
  do.call(.overlayList, lst)
}



.overlayList <- function(x, fun, filename="", recycle=TRUE, forcefun=FALSE, ...){

  ln <- length(x)
  if (ln < 1) {
    stop('no Rasters')
  }
  if (ln > 1) {
    compareRaster(x)
  }

  nl <- sapply(x, nlayers)
  maxnl <- max(nl)

  filename <- trim(filename)

  testmat <- NULL
  testlst <- vector(length=length(x), mode='list')
  #w <- getOption('warn')
  suppressWarnings({
  for (i in 1:length(testlst)) {
    v <- extract(x[[i]], 1:5)
    testmat <- cbind(testmat, as.vector(v))
    testlst[[i]] <- v
  }
  })
  #options('warn'= w)

  test1 <- try ( apply(testmat, 1, fun) , silent=TRUE )

  if ( !("try-error" %in% class(test1)) & (!forcefun)) {
    doapply <- TRUE
    if (! is.null(dim(test1))) {
      test1 <- t(test1)
    } else {
      test1 <- matrix(test1, ncol=maxnl)
    }
    nlout <- NCOL(test1)
  } else {
    doapply <- FALSE
    dovec <- FALSE
    test2 <- try ( do.call(fun, testlst), silent=TRUE )
    nlout <- length(test2)/5
    if ( "try-error" %in% class(test2) | length(test2) < 5) {
      dovec <- TRUE
      testlst <- lapply(testlst, as.vector)
      test3 <- try ( do.call(fun, testlst), silent=TRUE )
      nlout <- length(test3)/5
      if ("try-error" %in% class(test3) | length(test3) < 5) {
        stop('cannot use this formula, probably because it is not vectorized')
      }
    }
  }

  if (nlout == 1) {
    out <- raster(x[[1]])
  } else {
    out <- brick(x[[1]], values=FALSE, nl=nlout)
  }

  if ( canProcessInMemory(out, sum(nl)+maxnl) ) {
    pb <- pbCreate(3, label='overlay', ...)
    pbStep(pb, 1)
    if (doapply) {
      valmat <- matrix(nrow=ncell(out)*maxnl, ncol=length(x))
      for (i in 1:length(x)) {
        if (ncell(x[[i]]) < nrow(valmat)) {
          #options('warn'=-1)
          suppressWarnings({
          valmat[,i] <- as.vector(getValues(x[[i]])) * rep(1, nrow(valmat))
          })
          #options('warn'= w)
        } else {
          valmat[,i] <- as.vector(getValues(x[[i]]))
        }
      }
      pbStep(pb, 2)

      vals <- apply(valmat, 1, fun)
      if (! is.null(dim(vals))) {
        vals <- t(vals)
      }
      vals <- matrix(vals, nrow=ncell(out))

    } else {
      for (i in 1:length(x)) {
        x[[i]] <- getValues(x[[i]])
      }
      if (dovec) {
        x <- lapply(x, as.vector)
      }
      pbStep(pb, 2)
      vals <- do.call(fun, x)
      vals <- matrix(vals, nrow=ncell(out))
    }
    pbStep(pb, 3)
    out <- setValues(out, vals)
    if (filename != "") {
      out <- writeRaster(out, filename=filename, ...)
    }
    pbClose(pb)
    return(out)

  } else {

    if (filename == "") {
      filename <- rasterTmpFile()
    }
    out <- writeStart(out, filename=filename, ...)


    tr <- blockSize(out, n=sum(nl)+maxnl)
    pb <- pbCreate(tr$n, label='overlay', ...)

    if (doapply) {
      valmat = matrix(nrow=tr$nrows[1]*ncol(out)*maxnl, ncol=length(x))
      for (i in 1:tr$n) {
        if (i == tr$n) {
          valmat = matrix(nrow=tr$nrows[i]*ncol(out)*maxnl , ncol=length(x))
        }
        for (j in 1:length(x)) {
          v <- as.vector(getValues(x[[j]], row=tr$row[i], nrows=tr$nrows[i]))
          if (length(v) < nrow(valmat)) {
            # options('warn'=-1)
            suppressWarnings({
            valmat[,j] <- v * rep(1, nrow(valmat))
            })
            # options('warn'=w)
          } else {
            valmat[,j] <- v
          }
        }

        vv <- apply(valmat, 1, fun)
        if (! is.null(dim(vv))) {
          vals <- t(vv)
        }
        vv <- matrix(vv, ncol=nlout)
        out <- writeValues(out, vv, tr$row[i])
        pbStep(pb, i)
      }

    } else {
      vallist <- list()
      for (i in 1:tr$n) {
        if (dovec) {
          for (j in 1:length(x)) {
            vallist[[j]] <- as.vector( getValues(x[[j]], row=tr$row[i], nrows=tr$nrows[i]) )
          }
        } else {
          for (j in 1:length(x)) {
            vallist[[j]] <- getValues(x[[j]], row=tr$row[i], nrows=tr$nrows[i])
          }
        }

        vv <- do.call(fun, vallist)
        vv <- matrix(vv, ncol=nlout)
        out <- writeValues(out, vv, tr$row[i])
        pbStep(pb, i)
      }
    }
    pbClose(pb)
    out <- writeStop(out)
  }
  return(out)
}
