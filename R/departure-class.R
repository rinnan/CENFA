#' departure-class
#'
#' An object of class \code{departure} is created from performing ecological-niche
#' factor analysis on species presence data using the \code{enfa} function.
#'
#' @slot call Original function call
#' @slot df departure factor
#' @slot departure Magnitude of the departure factor
#' @slot ras Raster* object of transformed climate values
#' @slot weights Raster layer of weights used for departure calculation
#' @export

setClass("departure", slots = list(call = "call", df = "numeric", departure = "numeric", ras = "Raster", weights = "Raster"))

setMethod("show",
          signature = "departure",
          function(object){
            if (!inherits(object, "departure"))
              stop("Object of class 'departure' expected")
            cat("CLIMATIC DEPARTURE")
            cat("\n\nDeparture factor: \n")
            print(sort(round(object@df, 2), decreasing = T))
            cat("\nDeparture: ")
            cat(signif(object@departure, 4))
          }
)
