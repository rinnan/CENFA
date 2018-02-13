#' Predict methods
#'
#' Make a RasterLayer with predictions from a fitted model object.
#'
#' @param object model object
#' @param newdata optional new data
#' @param filename character optional filename
#' @param ... additional arguments
#'
#' @export

#' @rdname predict
setMethod("predict",
          signature(object = "cnfa"),
          function(object, newdata, filename = "", ...){

            x <- get(as.character(object@call$x))
            if (!is(x, "Raster")) x <- raster(x)
            nm <- names(x)
            if (!missing(newdata)){
              if (!all.equal(nm, names(newdata))) stop("layer names of newdata do not match layer names of model")
            }
            if (is.null(object@call$scale) | as.logical(as.character(object@call$scale))) {
              center <- cellStats(x, 'mean', na.rm = TRUE)
              sd <- cellStats(x, 'sd', na.rm = TRUE)
              if (missing(newdata)) x <- parScale(x, center = center, scale = sd)
            }

            if (!missing(newdata)) {
              if (is.null(object@call$scale) | as.logical(as.character(object@call$scale))) {
                x <- parScale(newdata, center = center, scale = sd)
              } else {
                x <- newdata
              }
            }

            m <- object@mf
            s <- object@sf
            filename <- trim(filename)
            if (!canProcessInMemory(x) && filename == '') {
              filename <- rasterTmpFile()
            }

            f1 <- function(x) (abs(x - m) %*%  s) / length(s)

            ras <- .calc(x, f1, forceapply = T, filename = filename, names = "Sensitivity", ...)

            return(ras)
          }
)

#' @rdname predict
setMethod("predict",
          signature(object = "enfa"),
          function(object, newdata, filename = "", ...){

            x <- get(as.character(object@call$x))
            if (!is(x, "Raster")) x <- raster(x)
            nm <- names(x)
            if (!missing(newdata)){
              if (!all.equal(nm, names(newdata))) stop("layer names of newdata do not match layer names of model")
            }
            U <- object@co
            if (is.null(object@call$scale) | as.logical(as.character(object@call$scale))) {
              center <- cellStats(x, 'mean', na.rm = TRUE)
              sd <- cellStats(x, 'sd', na.rm = TRUE)
              if (missing(newdata)) x <- parScale(x, center = center, scale = sd)
            }

            if (!missing(newdata)) {
              if (is.null(object@call$scale) | as.logical(as.character(object@call$scale))) {
                x <- parScale(newdata, center = center, scale = sd)
              } else {
                x <- newdata
              }
            }

            filename <- trim(filename)
            if (!canProcessInMemory(x) && filename == '') {
              filename <- rasterTmpFile()
            }

            ras <- .calc(x, function(x) {x %*% U}, forceapply = T, filename = filename, names = nm, ...)

            return(ras)
          }
)

#' @rdname predict
setMethod("predict",
          signature(object = "departure"),
          function(object, newdata, filename = "", ...){

            x <- get(as.character(object@call$x))
            if (is(x, "GLdeparture")) {
              x <- raster(x)
            } else if (is(x, "Raster")) {
              y <- get(as.character(object@call$y))
              if (is.null(object@call$center)) center <- TRUE else center <- as.logical(as.character(object@call$center))
              if (is.null(object@call$scale)) scale <- TRUE else scale <- as.logical(as.character(object@call$scale))
              gld <- GLdeparture(x = x, y = y, center = center, scale = scale)
              x <- raster(gld)
            }

            filename <- trim(filename)
            if (!canProcessInMemory(x) && filename == '') {
              filename <- rasterTmpFile()
            }

            d <- object@df
            f1 <- function(x) x %*% d

            ras <- .calc(x, fun = f1, forceapply = T, filename = filename, names = "Exposure", ...)

            return(ras)
          }
)

# mod <- cnfa(climdat.hist, s.dat = ABPR, field = "CODE", scale = T)
# as.logical(as.character(mod@call$scale))
#
# get(as.character(mod@call$x))
#
# glc <- GLcenfa(x = climdat.hist)
# mod2 <- cnfa(x = glc, s.dat = ABPR, field = "CODE")
#
# sensitivity_map(mod)
# test1 <- predict(mod2)
# test2 <- predict(mod2, newdata = climdat.fut)
# stretchPlot(test1, type = "sd", n = 2)
# stretchPlot(test2, type = "sd", n = 2)
#
# dep1 <- departure(x = climdat.hist, y = climdat.fut, s.dat = ABPR, field = "CODE")
# test3 <- predict(dep1)
# plot(test3)
