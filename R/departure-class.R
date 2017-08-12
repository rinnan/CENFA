#' departure-class
#'
#' An object of class \code{departure} is created from performing ecological-niche factor analysis on species presence data using the \code{enfa} function.
#'
#' @slot call marginality factor. Vector that describes the location of the species Hutchinsonian niche.
#' @slot departure Standardized magnitude of the marginality factor.
#' @slot distances specialization factors. Matrix of specializations.
#' @slot departure_ras Raster* object of transformed climate values, with number of layers equal to nf + 1.
#' @slot present Number of cells in which species is present.
#' @export

setClass("departure", slots = list(call = "call", departure = "numeric", distances = "numeric", departure_ras = "Raster", present = "numeric"))

setMethod ("show", "departure", function(object){
  if (!inherits(object, "departure"))
    stop("Object of class 'departure' expected")
  cat("CLIMATIC DEPARTURE")
  cat("\ndeparture: ")
  cat(signif(object@departure, 4))
  cat("\nnumber of presence cells: ")
  cat(object@present)
}
)
