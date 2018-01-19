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
#' @param parallel logical. If \code{TRUE} then multiple cores are utilized
#' @param n numeric. Optional number of CPU cores to utilize for parallel processing
#'
#' @examples
#' mat <- parScale(x = climdat.hist)
#'
#' @return Raster* object
#'
#' @seealso \code{\link[base]{scale}}, \code{\link[raster]{scale}}
#'
#' @export
#' @importFrom pbapply pbsapply pboptions

setGeneric("parScale", function(x, center = FALSE, scale = FALSE, parallel = FALSE, n){
  standardGeneric("parScale")})

#' @rdname parScale
setMethod("parScale",
          signature(x = "Raster"),
          function(x, center = FALSE, scale = FALSE, parallel = FALSE, n){

            if(canProcessInMemory(x) & !parallel){
              v <- values(x)
              x <- setValues(x, scale(v, center=center, scale=scale))
              return(x)
            }

            if(!parallel){
              if (!is.logical(center)) {

                stopifnot(length(center) == nlayers(x))
                x <- x - center

              } else if (center) {
                m <- cellStats(x, 'mean', na.rm=TRUE)
                x <- x - m
              }

              if (!is.logical(scale)) {
                stopifnot(length(scale) == nlayers(x))
                x <- x / scale

              } else if (scale) {
                if (center[1] & is.logical(center[1])) {
                  st <- cellStats(x, 'sd', na.rm=TRUE)
                } else {
                  st <- cellStats(x, 'rms', na.rm=TRUE)
                }
                x <- x / st
              }
              return(x)
            }

            nl <- nlayers(x)
            s <- 1:nl
            if(is.logical(center)) center <- rep(center, nl)
            if(is.logical(scale)) scale <- rep(scale, nl)

            if (missing(n)) {
              n <- parallel::detectCores()
              message(n, ' cores detected, using ', n-1)
              n <- n-1
            }
            cl <- snow::makeCluster(getOption("cl.cores", n))
            snow::clusterExport(cl, c("x", "s", "scale", "center"),
                                envir = environment())
            doSNOW::registerDoSNOW(cl)
            pb <- txtProgressBar(min = 0, max = length(s), style = 3, char = "-")
            progress <- function(n) setTxtProgressBar(pb, n)
            opts <- list(progress = progress)
            result <- foreach::foreach(i = s, .options.snow = opts) %dopar% {
              do.call(raster::scale, list(x = x[[ i ]], center = center[i], scale = scale[i]))
            }
            close(pb)
            snow::stopCluster(cl)

            # resolves error message for global binding of i
            for(i in s){}

            x <- brick(result)

            closeAllConnections()
            return(x)
          }
)
