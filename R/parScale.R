#' Efficient scaling of Raster* objects
#'
#' \code{parScale} expands the \code{raster::scale} function to allow for
#' faster parallel processing, scaling each layer of \code{x} in parallel.
#'
#' @param x Raster* object
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
#' @param filename character. Optional filename to save the Raster* output to
#'   file. If this is not provided, a temporary file will be created for large \code{x}
#' @param progress logical. If \code{TRUE}, messages and progress bar will be
#'   printed
#' @param parallel logical. If \code{TRUE} then multiple cores are utilized
#' @param n numeric. Number of CPU cores to utilize for parallel processing
#' @param ... Additional arguments for \code{\link[raster]{writeRaster}}
#'
#' @examples
#' ch.scale <- parScale(x = climdat.hist)
#'
#' @return Raster* object
#'
#' @seealso \code{\link[base]{scale}}, \code{\link[raster]{scale}}
#'
#' @export
#' @importFrom pbapply pbsapply pboptions
#' @importFrom foreach '%dopar%'
#' @importFrom raster subset

setGeneric("parScale", function(x, ...){
  standardGeneric("parScale")})

#' @rdname parScale
setMethod("parScale",
          signature(x = "Raster"),
          function(x, center = TRUE, scale = TRUE, filename = '', progress = FALSE, parallel = FALSE, n = 1, ...){

            if (is.logical(center) & is.logical(scale)) {
              if (!center & !scale) return(x) }

            if (canProcessInMemory(x) & !parallel) {
              v <- values(x)
              x <- setValues(x, scale(v, center = center, scale = scale))
              return(x)
            }

            filename <- trim(filename)
            if (filename == '') filename <- rasterTmpFile()

            if (!parallel) {
              x <- scale(x, center = center, scale = scale)
              writeRaster(x, filename = filename, ...)
              return(x)
            }

            on.exit(closeAllConnections())
            nl <- nlayers(x)
            s <- 1:nl
            if (is.logical(center)) center <- rep(center, nl)
            if (is.logical(scale)) scale <- rep(scale, nl)

            if (n < 1 | !is.numeric(n)) {
              n <- min(parallel::detectCores() - 1, nl)
              if (progress) message('incorrect number of cores specified, using ', n)
            } else if (n > parallel::detectCores()) {
              n <- min(parallel::detectCores() - 1, nl)
              if (progress) message('too many cores specified, using ', n)
            }
            cl <- snow::makeCluster(getOption("cl.cores", n))
            snow::clusterExport(cl, c("x", "s", "scale", "center", "subset"),
                                envir = environment())
            doSNOW::registerDoSNOW(cl)
            if (progress) {
              pb <- txtProgressBar(min = 0, max = length(s), style = 3, char = "-")
              progress <- function(n) setTxtProgressBar(pb, n)
              opts <- list(progress = progress)
              result <- foreach::foreach(i = s, .options.snow = opts) %dopar% {
                do.call(raster::scale, list(x = raster::subset(x, i), center = center[i], scale = scale[i]))
              }
              close(pb)
            } else if (!progress) {
              result <- foreach::foreach(i = s) %dopar% {
                do.call(raster::scale, list(x = raster::subset(x, i), center = center[i], scale = scale[i]))
              }
            }
            snow::stopCluster(cl)

            # resolves error message for global binding of i
            for(i in s){}

            x <- brick(result)
            writeRaster(x, filename = filename)#, ...)

            closeAllConnections()
            return(x)
          }
)
