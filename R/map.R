#' Contrast adjustments for RasterLayer plots
#'
#' A plotting function that provides methods for improving the contrast
#' between values.
#'
#'
#' @param x a RasterLayer
#' @param type character. Possible values are "linear", "stretch", and "sd"
#' @param n number of standard deviations to include if \code{type = "sd"}
#' @param legend.args a list of arguments to pass on to adjust the plot legend.
#'   See \code{\link[fields]{image.plot}} for possible options
#' @param ... Additional arguments for raster::plot
#'
#' @details
#' If \code{type = "stretch"}, a histogram equalization procedure will be
#' applied to the values of \code{x}. If \code{type = "sd"}, the values of
#' \code{x} will be scaled between values that fall between \code{n} standard
#' deviations of the mean.
#'
#' @examples
#' mdr <- climdat.hist[[1]]
#' map(mdr)
#' map(mdr, type = "stretch")
#' map(mdr, type = "sd", n = 2)
#'
#' @importFrom stats ecdf sd
#'
#'
#' @export

setGeneric("map", function(x, ...) {
  standardGeneric("map")
})

#' @rdname map
setMethod("map",
          signature(x = "RasterLayer"),
          function(x, type = "linear", n, legend.args, ...) {

            if(type == "stretch") y <- .stretch(x, type = "hist.equal")
            if(type == "linear") y <- x
            if(type == "sd") y <- .stretch(x, type = "sd", n = n)

            if(missing(legend.args)) legend.args <- list()
            legend.args$x <- x
            legend.args$legend.only <- TRUE

            plot(y, legend = F, ...)
            do.call(plot, args = legend.args)
          }
)


#' @keywords internal
.stretch <- function(x, type = "hist.equal", n) {
  if (type == "hist.equal") {
    ecdfun <- stats::ecdf(getValues(x))
    y <- calc(x, fun = function(x) round(ecdfun(x)*255, 0))
    return(y)
  } else if (type == "sd") {
    if (missing(n)) stop("number of standard deviations not specified")
    x.sd <- cellStats(x, sd)
    x.mean <- cellStats(x, mean)
    x.max <- x.mean + x.sd*n
    x.min <- x.mean - x.sd*n
    sdfun <- function(y) {
      y[y < 0] <- x.min
      y[y > x.max] <- x.max
      tt <- (y - x.min) / (x.max - x.min)
      return(tt*255)
    }
    return(calc(x, fun = function(x) sdfun(x)))
  }

}
