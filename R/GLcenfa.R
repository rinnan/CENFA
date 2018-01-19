#' GLcenfa
#'
#' This function is used to facilitate comparisons between species in the same
#' study area. It speeds up the computation of multiple CNFAs or ENFAs by calculating
#' the global covariance matrix as a first step, which can then be fed into the
#' \code{\link{cnfa}} or \code{\link{enfa}} functions as their first argument.
#' This saves the user from having to calculate the global covariance matrix for
#' each species, which can take quite a bit of time.
#'
#' @aliases print.GLcenfa, show.GLcenfa
#'
#' @param x Raster* object, typically a brick or stack of p environmental raster
#'   layers
#' @param center logical or numeric. If \code{TRUE}, centering is done by
#'   subtracting the layer means (omitting NAs), and if \code{FALSE}, no centering
#'   is done. If \code{center} is a numeric vector with length equal to the
#'   \code{nlayers(x)}, then each layer of \code{x} has the corresponding value
#'   from center subtracted from it
#' @param scale logical or numeric. If \code{TRUE}, scaling is done by dividing
#'   the (centered) layers of \code{x} by their standard deviations if center is
#'   \code{TRUE}, and the root mean square otherwise. If scale is \code{FALSE},
#'   no scaling is done. If scale is a numeric vector with length equal to
#'   \code{nlayers(x)}, each layer of \code{x} is divided by the corresponding
#'   value. Scaling is done after centering
#' @param progress logical. If \code{TRUE} then progress updates are printed
#' @param parallel logical. If \code{TRUE} then multiple cores are utilized
#' @param n numeric. Optional number of CPU cores to utilize for parallel processing
#' @param filename character. Optional filename to save the RasterBrick output
#'   to file. If this is not provided, a temporary file will be created for large
#'   \code{x}
#' @param ... Additonal arguments for \code{\link[raster]{writeRaster}}
#'
#' @examples
#' glc <- GLcenfa(x = climdat.hist)
#'
#' @return Returns an S4 object of class \code{GLcenfa} with the following components:
#' \describe{
#'   \item{global_ras}{Raster* \code{x} of p layers, possibly centered and scaled}
#'   \item{cov}{Global p x p covariance matrix}
#'   }
#   \item{center}{Vector of layer means of \code{x}}
#   \item{sd}{Vector of layer standard deviations of \code{x}}
# }
#'
#' @seealso \code{\link{cnfa}}, \code{\link{enfa}}
#' @export

setGeneric("GLcenfa", function(x, center = TRUE, scale = TRUE, filename = '', progress = TRUE, parallel = FALSE, n, ...) {
  standardGeneric("GLcenfa")
})

#' @rdname GLcenfa
setMethod("GLcenfa",
          signature(x = "Raster"),
          function(x, center = TRUE, scale = TRUE, filename = '', progress = TRUE, parallel = FALSE, n, ...){

            filename <- trim(filename)
            if (!canProcessInMemory(x) && filename == '') {
              filename <- rasterTmpFile()
            }

            # if(progress) cat("\nCalculating layer means...")
            # if(center) cent <- cellStats(x, 'mean')
            # if(progress) cat("\nCalculating layer sds...")
            # sds <- cellStats(x, 'sd')

            x <- parScale(x, center = center, scale = scale, parallel = parallel, n = n)

            if(!(filename == "")){
              cat("Writing data to file...")
              writeRaster(x, filename = filename, ...)
            }
            cov.mat <- parCov(x, center = FALSE, scale = FALSE, parallel = parallel, n = n)

            GLcenfa <- methods::new("GLcenfa", global_ras = x, cov = cov.mat)
            return(GLcenfa)
          }
)

