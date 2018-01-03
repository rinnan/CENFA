#' @useDynLib CENFA
#' @importFrom Rcpp sourceCpp
#' @import raster
#' @import sp
#' @importFrom magrittr %>%
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

.onUnload <- function (libpath) {
  library.dynam.unload("CENFA", libpath)
}
