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
            if (is.null(object@call$scale) || as.logical(as.character(object@call$scale))) {
              center <- cellStats(x, 'mean', na.rm = TRUE)
              sd <- cellStats(x, 'sd', na.rm = TRUE)
              if (missing(newdata)) x <- parScale(x, center = center, scale = sd)
            }

            if (!missing(newdata)) {
              if (is.null(object@call$scale) || as.logical(as.character(object@call$scale))) {
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
            if (is.null(object@call$scale) || as.logical(as.character(object@call$scale))) {
              center <- cellStats(x, 'mean', na.rm = TRUE)
              sd <- cellStats(x, 'sd', na.rm = TRUE)
              if (missing(newdata)) x <- parScale(x, center = center, scale = sd)
            }

            if (!missing(newdata)) {
              if (is.null(object@call$scale) || as.logical(as.character(object@call$scale))) {
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
          function(object, filename = "", ...){

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
            f1 <- function(x) (x %*% d) / length(d)

            ras <- .calc(x, fun = f1, forceapply = T, filename = filename, names = "Exposure", ...)

            return(ras)
          }
)

#' @rdname predict
setMethod("predict",
          signature(object = "vulnerability"),
          function(object, newdata, filename = "", ...){

            x <- get(as.character(object@call$cnfa))
            y <- get(as.character(object@call$dep))
            if (is.null(object@call$w)) {
              w <- c(1, 1)
            } else {
              w <- as.numeric(object@call$w)
            }

            if (is.null(object@call$method)) {
              method <- "geometric"
            } else {
              method <- as.character(object@call$method)
            }

            ifelse(missing(newdata),
                   s.map <- predict(x),
                   s.map <- predict(x, newdata = newdata))
            e.map <- predict(y)

            if (method == "arithmetic") {
              ras <- (w[1] * s.map + w[2] * e.map) / sum(w)
            } else if (method == "geometric") {
              if(w[1] == w[2]) {
                ras <- sqrt(s.map * e.map)
              } else {
                w <- w / sum(w)
                ras <- ( (s.map^w[1]) * (e.map^w[2]) )^(1/sum(w))
              }
            }
            return(ras)
          }
)
