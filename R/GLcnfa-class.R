#' global.cnfa-class
#'
#' An object of class \code{global.cnfa} is created from performing climate-niche factor analysis on species presence data using the \code{cnfa} function.
#'
#' @slot  call original function call
#' @slot nf number of kept specialization factors
#' @slot m marginalities
#' @slot s specialization factors
#' @slot map Raster map of habitat suitability
#' @export

setClass("GLcnfa", slots = list(global_ras = "Raster", cov = "matrix", center = "numeric", sd = "numeric", ncells = "numeric"))

setMethod ("show", "GLcnfa", function(object){
  if (!inherits(object, "GLcnfa"))
    stop("Object of class 'GLcnfa' expected")
  cat("GLcnfa")
  cat("\n")
  cat("\nNumber of cells with climate data: ")
  cat(object@ncells)
  cat("\n")
  cat("\n")
  print(object@global_ras)
}
)
