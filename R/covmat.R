#' Efficient calculation of covariance matrices for Raster* objects
#'
#' \code{covmat} efficiently calculates the covariance between different climate
#' and ecological variables, taking advantage of parallel processing and pulling
#' data into memory only as necessary. For large datasets with lots of variables,
#' calculating the covariance matrix rapidly becomes unwieldy, as the number of
#' calculations required grows quadratically with the number of variables.
#'
#' @param x Raster* object, typically a brick or stack of climate raster layers
#' @param y NULL (default) or a Raster* object with the same extent and resolution
#'   as \code{x}
#' @param ... additional arguments, including any of the following:
#' @param center logical. If \code{TRUE}, the Raster* object will get centered
#'   before calculating the covariance
#' @param scale logical. If \code{TRUE}, the Raster* object will get scaled
#'   before calculating the covariance
#' @param w optional Raster* object of weights for a weighted covariance matrix
#' @param sample logical. If \code{TRUE}, the sample covariance is calculated
#'   with a denominator of $n-1$
#' @param parallel logical. If \code{TRUE} then multiple cores are utilized
#' @param n numeric. Optional number of CPU cores to utilize for parallel processing
#'
#' @examples
#' mat1 <- covmat(x = climdat.hist)
#' mat2 <- covmat(x = climdat.hist, center = TRUE, scale = TRUE)
#'
#' # covariance between two Raster* objects
#' mat3 <- covmat(x = climdat.hist, y = climdat.fut)
#' dat.h <- values(climdat.hist)
#' dat.f <- values(climdat.fut)
#' mat4 <- cov(dat.h, dat.f, use = "na.or.complete", method = "pearson")
#'
#' # same results either way
#' all.equal(mat3, mat4)
#'
#' @return Returns a matrix with the same row and column names as the layers of
#'   \code{x}. If \code{y} is supplied, then the covariances between the layers
#'   of \code{x} and the layers of code{y} are computed.
#'
#' @seealso \code{\link[stats]{cov}}
#'
#' @export
#' @importFrom pbapply pbsapply pboptions

setGeneric("covmat", function(x, y, ...){
  standardGeneric("covmat")})

#' @rdname covmat
setMethod("covmat",
          signature(x = "Raster", y = "missing"),
          function(x, center = FALSE, scale = FALSE, w = NULL, sample = TRUE, parallel = FALSE, n){

            small <- canProcessInMemory(x)
            if(small){
              dat <- na.omit(values(x))
              if(center) dat <- dat - colMeans(dat)
              if(scale) {
                sds <- apply(dat, 2, sd)
                dat <- dat/sds}
              mat <- cov(dat, method = "pearson")
              return(mat)
            }
            nl <- nlayers(x)
            mat <- matrix(NA, nrow = nl, ncol = nl)
            colnames(mat) <- rownames(mat) <- names(x)

            ii <- rep(1, nl)
            for(i in 2:nl) ii <- c(ii, rep(i, each = (nl - i + 1)))
            jj <- 1:nl
            for(i in 2:nl) jj <- c(jj, i:nl)
            #combs <- cbind(ii, jj)
            s <- 1:length(ii)

            if(center){
              means <- cellStats(x, stat = 'mean', na.rm = T)
              x <- (x - means)
            }

            if(scale){
              sds <- cellStats(x, stat = 'sd', na.rm = T)
              x <- x/sds
            }

            if(!parallel){
              pboptions(char = "-", txt.width = NA, type = "txt")
              result <- pbsapply(s, function(p) do.call(.covij, list(x = x[[ ii[p] ]], y = x[[ jj[p] ]], w = w, sample = sample)))
            }

            if(parallel){
              if (missing(n)) {
                n <- parallel::detectCores()
                message(n, ' cores detected, using ', n-1)
                n <- n-1
              }
              cl <- snow::makeCluster(getOption("cl.cores", n))
              snow::clusterExport(cl, c(".covij", "raster", "cellStats", "x", "ii", "jj", "s", "w", "canProcessInMemory", "values", "sample"),
                                  envir = environment())
              doSNOW::registerDoSNOW(cl)
              pb <- txtProgressBar(min = 0, max = length(s), style = 3, char = "-")
              progress <- function(n) setTxtProgressBar(pb, n)
              opts <- list(progress = progress)
              result <- foreach::foreach(p = s, .combine=c, .options.snow = opts) %dopar% {
                do.call(.covij, list(x = x[[ ii[p] ]], y = x[[ jj[p] ]], w = w, sample = sample))
              }
              close(pb)
              snow::stopCluster(cl)
            }

            for(p in s){
              mat[ii[p], jj[p]] <- mat[jj[p], ii[p]] <- result[p]
            }

            closeAllConnections()
            return(mat)
          }
)

#' @rdname covmat
setMethod("covmat",
          signature(x = "Raster", y = "Raster"),
          function(x, y, w = NULL, sample = TRUE, parallel = FALSE, n){

            small <- canProcessInMemory(x)
            if(small){
              x.dat <- values(x)
              y.dat <- values(y)
              mat <- cov(x.dat, y.dat, method = "pearson", use = "na.or.complete")
              return(mat)
            }

            nlx <- nlayers(x)
            nly <- nlayers(y)
            mat <- matrix(NA, nrow = nlx, ncol = nly)
            rownames(mat) <- names(x)
            colnames(mat) <- names(y)
            z <- .expand.grid.unique(1:nlx, 1:nly)
            s <- 1:nrow(z)

            if(!parallel){
              pboptions(char = "-", txt.width = NA, type = "txt")
              result <- pbsapply(s, function(p) do.call(.covij, list(x = x[[ z[p, 1] ]], y = y[[ z[p, 2] ]], w = w, sample = sample)))
            }

            if(parallel){
              if (missing(n)) {
                n <- parallel::detectCores()
                message(n, ' cores detected, using ', n-1)
                n <- n-1
              }
              cl <- snow::makeCluster(getOption("cl.cores", n))
              snow::clusterExport(cl, c(".covij", "raster", "cellStats", "x", "y", "z", "s", "w", "canProcessInMemory", "values", "sample"),
                                  envir = environment())
              doSNOW::registerDoSNOW(cl)
              pb <- txtProgressBar(min = 0, max = length(s), style = 3, char = "-")
              progress <- function(n) setTxtProgressBar(pb, n)
              opts <- list(progress = progress)
              result <- foreach::foreach(p = s, .combine=c, .options.snow = opts) %dopar% {
                do.call(.covij, list(x = x[[ z[p, 1] ]], y = y[[ z[p, 2] ]], w = w, sample = sample))
              }
              close(pb)
              snow::stopCluster(cl)
            }

            for(p in s){
              mat[ z[p,2], z[p,1] ] <- result[p]
              if(nlx > 1 & nly > 1) mat[ z[p,1], z[p,2] ] <- mat[ z[p,2], z[p,1] ]
            }

            closeAllConnections()
            return(mat)
          }
)


# covmat <- function(x, cores = 1, center = FALSE, scale = FALSE, progress = TRUE, sample = TRUE) {
#   stopifnot(is.numeric(cores) & cores >= 0)
#
#   small <- raster::canProcessInMemory(x)
#   if(small){
#     dat <- stats::na.omit(raster::values(x))
#     if(center) dat <- dat - colMeans(dat)
#     if(scale) {
#       sds <- apply(dat, 2, sd)
#       dat <- dat/sds}
#     mat <- cov(dat, method = "pearson")
#     return(mat)
#   }
#   nl <- raster::nlayers(x)
#   mat <- matrix(NA, nrow=nl, ncol=nl)
#   colnames(mat) <- rownames(mat) <- names(x)
#
#   ii <- rep(1,nl)
#   for(i in 2:nl) ii <- c(ii, rep(i, each = (nl - i + 1)))
#   jj <- 1:nl
#   for(i in 2:nl) jj <- c(jj, i:nl)
#   s <- 1:length(ii)
#
#   if(center){
#     means <- raster::cellStats(x, stat = 'mean', na.rm = T)
#     x <- (x - means)
#   }
#
#   if(scale){
#     sds <- raster::cellStats(x, stat = 'sd', na.rm = T)
#     x <- x/sds
#   }
#
#   if(cores == 1){
#     if(progress){
#       pbapply::pboptions(type = "txt", char = "=", txt.width = NA)
#       result <- pbapply::pbsapply(s, function(p) do.call(.covij, list(x = x, i = ii[p], j = jj[p])))
#     } else {result <- sapply(s, function(p) do.call(.covij, list(x = x, i = ii[p], j = jj[p])))}
#   }
#
#   if(cores > 1){
#     cl <- snow::makeCluster(getOption("cl.cores", cores))
#     snow::clusterExport(cl, c(".covij", "raster", "cellStats", "x", "ii", "jj", "s", "nl", "canProcessInMemory", "values"), envir = environment())
#     doSNOW::registerDoSNOW(cl)
#     if(progress){
#       pb <- utils::txtProgressBar(min = 0, max = length(s), style = 3)
#       progress <- function(n) utils::setTxtProgressBar(pb, n)
#       opts <- list(progress = progress)
#       result <- foreach::foreach(p = s, .combine=c, .options.snow = opts) %dopar% {
#         do.call(.covij, list(x = x, i = ii[p], j = jj[p]))
#       }
#       close(pb)
#     }   else if(!progress){
#       result <- foreach::foreach(p = s, .combine = c, .options.snow = opts) %dopar% {
#         do.call(.covij, list(x = x, i = ii[p], j = jj[p]))
#       }
#     }
#     snow::stopCluster(cl)
#   }
#
#   for(p in s){
#     mat[ii[p], jj[p]] <- mat[jj[p], ii[p]] <- result[p]
#   }
#
#   closeAllConnections()
#   return(mat)
# }

#' @keywords internal
.expand.grid.unique <- function(x, y){

  nx <- length(x)
  ny <- length(y)

  if(ny == 1){
    dat <- cbind(x, y)
  } else if(nx <= ny) {
    dat <- NULL
    for(i in 1:nx){
      for(j in i:ny){
        dat <- rbind(dat, c(x[i], y[j]))
      }
    }
  } else if(nx > ny) {
    dat <- NULL
    for(j in 1:ny){
      for(i in j:nx){
        dat <- rbind(dat, c(x[i], y[j]))
      }
    }
  }
  return(dat)
}
