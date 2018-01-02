#' departure-class
#'
#' An object of class \code{departure} is created from performing ecological-niche factor analysis on species presence data using the \code{enfa} function.
#'
#' @slot call Original function call.
#' @slot df departure factor.
#' @slot departure Magnitude of the departure factor.
#' @slot present Number of cells in which species is present.
#' @slot ras Raster* object of transformed climate values, with number of layers equal to nf + 1.
#' @export

setClass("departure", slots = list(call = "call", df = "numeric", departure = "numeric", present = "numeric", ras = "Raster"))

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
            cat("\nNumber of cells present: ")
            cat(object@present)
          }
)
