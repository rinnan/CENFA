#' Efficient calculation of covariance matrices
#'
#' \code{covmat} efficiently calculates the covariance between different climate and ecological variables, taking advantage of parallel processing and pulling data into memory only as necessary. For large datasets with lots of variables, calculating the covariance matrix rapidly becomes unwieldy, as the number of calculations required grows quadratically with the number of variables.
#'
#' @param x Raster* object, typically a brick or stack of climate raster layers
#' @param cores numeric. Number of CPU cores you wish to utilize for parallel processing.
#' @param center logical. If \code{TRUE}, the Raster* object will get centered before calculating the covariance.
#' @param scale logical. If \code{TRUE}, the Raster* object will get scaled before calculating the covariance.
#' @param progress logical. If \code{TRUE}, a progress bar will be displayed.
#' @param sample logical. If \code{TRUE}, the sample covariance is calculated with a denominator of $n-1$.
#' @return Returns a matrix with the same row and column names as the layers of the Raster* object.
#' @export
#' @importFrom pbapply pbsapply pboptions


covmat <- function(x, cores = 1, center = FALSE, scale = FALSE, progress = TRUE, sample = TRUE) {
  stopifnot(is.numeric(cores) & cores >= 0)

  small <- canProcessInMemory(x)
  if(small){
    dat <- na.omit(values(x))
    mat <- cov(dat, method = "pearson")
    return(mat)
  }
  nl <- nlayers(x)
  mat <- matrix(NA, nrow=nl, ncol=nl)
  colnames(mat) <- rownames(mat) <- names(x)

  ii <- rep(1,nl)
  for(i in 2:nl) ii <- c(ii, rep(i, each = (nl - i + 1)))
  jj <- 1:nl
  for(i in 2:nl) jj <- c(jj, i:nl)
  s <- 1:length(ii)

  if(center){
    means <- cellStats(x, stat = 'mean', na.rm = T)
    x <- (x - means)
  }

  if(scale){
    sds <- cellStats(x, stat = 'sd', na.rm = T)
    x <- x/sds
  }

  if(cores == 1){
    if(progress){
      pbapply::pboptions(type = "txt", char = "=", txt.width = NA)
      result <- pbapply::pbsapply(s, function(p) do.call(.covij, list(x = x, i = ii[p], j = jj[p])))
    } else {result <- sapply(s, function(p) do.call(.covij, list(x = x, i = ii[p], j = jj[p])))}
  }

  if(cores > 1){
    cl <- makeCluster(getOption("cl.cores", cores))
    clusterExport(cl, c(".covij", "raster", "cellStats", "x", "ii", "jj", "s", "nl", "canProcessInMemory", "values"), envir = environment())
    registerDoSNOW(cl)
    if(progress){
      pb <- txtProgressBar(min = 0, max = length(s), style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      result <- foreach(p = s, .combine=c, .options.snow = opts) %dopar% {
        do.call(.covij, list(x = x, i = ii[p], j = jj[p]))
      }
      close(pb)
    }   else if(!progress){
      result <- foreach(p = s, .combine = c, .options.snow = opts) %dopar% {
        do.call(.covij, list(x = x, i = ii[p], j = jj[p]))
      }
    }
    stopCluster(cl)
  }

  for(p in s){
    mat[ii[p], jj[p]] <- mat[jj[p], ii[p]] <- result[p]
  }

  closeAllConnections()
  return(mat)
}
