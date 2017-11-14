#' departure-class
#'
#' An object of class \code{departure} is created from performing ecological-niche factor analysis on species presence data using the \code{enfa} function.
#'
#' @slot call marginality factor. Vector that describes the location of the species Hutchinsonian niche.
#' @slot df departure factor.
#' @slot departure Standardized magnitude of the marginality factor.
#' @slot departure_ras Raster* object of transformed climate values, with number of layers equal to nf + 1.
#' @slot present Number of cells in which species is present.
#' @export

# setClass("departure", slots = list(call = "call", departure = "numeric", departure_ras = "Raster", present = "numeric"))
#
# setMethod("show",
#           signature = "departure",
#           function(depart){
#             if (!inherits(depart, "departure"))
#               stop("Object of class 'departure' expected")
#             cat("CLIMATIC DEPARTURE")
#             cat("\ndeparture: ")
#             cat(signif(depart@departure, 4))
#             cat("\nnumber of cells present: ")
#             cat(depart@present)
#           }
# )
