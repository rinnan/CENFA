#' CENFA: Climate- and Ecological-Niche Factor Analysis
#'
#' The CENFA package provides three categories of important functions:
#' foo, bar and baz.
#'
#' @section Foo functions:
#' The foo functions ...
#'
#' @docType package
#' @name CENFA-package
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
