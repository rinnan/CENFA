#' difRaster
#'
#' An object of class \code{enfa} is created from performing ecological-niche
#' factor analysis on species presence data using the \code{enfa} function.
#'
#' @param x.hist Raster* object of historical climate values
#' @param x.fut  Raster* object of future climate values
#' @param scale logical. If \code{TRUE} then the values of \code{x.hist} and
#'   \code{x.fut} will be centered and scaled by the means and sds of the historical
#'   climate data
#' @param parallel logical. If \code{TRUE} then multiple cores are utilized
#' @param n numeric. Optional number of CPU cores to utilize for parallel processing
#'
#' @return A Raster* object of class \code{difRasterBrick}


setGeneric("difRaster", function(x.hist, x.fut, scale = TRUE, parallel = FALSE, n){
  standardGeneric("difRaster")})

#' @rdname difRaster
setMethod("difRaster",
          signature(x.hist = "Raster", x.fut = "Raster"),
          function(x.hist, x.fut, scale = TRUE, parallel = FALSE, n){

            if(names(x.hist) != names(x.fut)) stop("historical and future raster layers do not match")
            if(!all.equal(x.hist, x.fut)) stop("historical and future rasters resolutions or extent do not match")

            if(parallel == T) {
              if(scale == T) {
                if (missing(n)) {
                  n <- parallel::detectCores()
                  message(n, ' cores detected, using ', n-1)
                  n <- n-1
                }
                means <- cellStats(x.hist, mean)
                sds <- cellStats(x.hist, sd)
                beginCluster(n = n)
                x.hist <- clusterR(x.hist, fun = scale, export = list("means", "sds"), args = list(center = means, scale = sds))
                x.fut <- clusterR(x.fut, scale, export = list("means", "sds"), args = list(center = means, scale = sds))
                endCluster()
              }
            } else {
              if(scale == T) {
                means <- cellStats(x.hist, mean)
                sds <- cellStats(x.hist, sd)
                x.hist <- scale(x.hist)
                x.fut <- scale(x.fut, center = means, scale = sds)
              }
            }

            x.dif <- abs(x.fut - x.hist)
            names(x.dif) <- names(x.hist)
            x.dif <- methods::new('difRasterBrick',
                                  file = x.dif@file,
                                  data = x.dif@data,
                                  legend = x.dif@legend,
                                  title = x.dif@title,
                                  extent = x.dif@extent,
                                  rotated = x.dif@rotated,
                                  rotation = x.dif@rotation,
                                  ncols = x.dif@ncols,
                                  nrows = x.dif@nrows,
                                  crs = x.dif@crs,
                                  history = x.dif@history,
                                  z = x.dif@z)
            return(x.dif)
          }
)
