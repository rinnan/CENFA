#' global.cnfamad-class
#'
#' An object of class \code{global.cnfamad} is created from performing climate-niche factor analysis on species presence data using the \code{cnfamad} function.
#'
#' @slot  call original function call
#' @slot nf number of kept specialization factors
#' @slot m marginalities
#' @slot s specialization factors
#' @slot map Raster map of habitat suitability
#' @export

setClass("GLcnfamad",slots=list(global_ras = "Raster", cov = "matrix", center = "numeric", mad = "numeric", ncells = "numeric", scale = "logical"))

setMethod ("show", "GLcnfamad", function(object){
  if (!inherits(object, "GLcnfamad"))
    stop("Object of class 'GLcnfamad' expected")
  cat("GLcnfa")
  cat("\n")
  cat("\nNumber of cells with climate data: ")
  cat(object@ncells)
  cat("\n")
  cat("\n")
  print(object@global_ras)
}
)
