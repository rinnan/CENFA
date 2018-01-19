#' difRaster
#'
#' Calculates the absolute differences between two Raster* objects.
#'
#' @param x Raster* object of historical climate values
#' @param y  Raster* object of future climate values
#' @param scale logical. If \code{TRUE} then the values of \code{x} and
#'   \code{y} will be centered and scaled by the means and sds of the historical
#'   climate data
#' @param parallel logical. If \code{TRUE} then multiple cores are utilized
#' @param n numeric. Optional number of CPU cores to utilize for parallel processing
#'
#' @examples
#' difRaster(x = climdat.hist, y = climdat.fut)
#'
#' @return Raster* object
#' @export


setGeneric("difRaster", function(x, y, scale = TRUE, parallel = FALSE, n){
  standardGeneric("difRaster")})

#' @rdname difRaster
setMethod("difRaster",
          signature(x = "Raster", y = "Raster"),
          function(x, y, scale = TRUE, parallel = FALSE, n){

            if(!all.equal(names(x), names(y))) stop("historical and future raster layers do not match")
            if(!compareRaster(x, y)) stop("historical and future rasters resolutions or extent do not match")

            if(parallel == T) {
              if(scale == T) {
                if (missing(n)) {
                  n <- parallel::detectCores()
                  message(n, ' cores detected, using ', n-1)
                  n <- n-1
                }
                means <- cellStats(x, mean)
                sds <- cellStats(x, sd)
                beginCluster(n = n)
                x <- clusterR(x, fun = scale, export = list("means", "sds"), args = list(center = means, scale = sds))
                y <- clusterR(y, scale, export = list("means", "sds"), args = list(center = means, scale = sds))
                endCluster()
              }
            } else {
              if(scale == T) {
                means <- cellStats(x, mean)
                sds <- cellStats(x, sd)
                x <- scale(x)
                y <- scale(y, center = means, scale = sds)
              }
            }

            x.dif <- abs(y - x)
            names(x.dif) <- names(x)
            # x.dif <- methods::new('difRasterBrick',
            #                       file = x.dif@file,
            #                       data = x.dif@data,
            #                       legend = x.dif@legend,
            #                       title = x.dif@title,
            #                       extent = x.dif@extent,
            #                       rotated = x.dif@rotated,
            #                       rotation = x.dif@rotation,
            #                       ncols = x.dif@ncols,
            #                       nrows = x.dif@nrows,
            #                       crs = x.dif@crs,
            #                       history = x.dif@history,
            #                       z = x.dif@z)
            return(x.dif)
          }
)
