#' @keywords internal

.onUnload <- function (libpath) {
  library.dynam.unload("CENFA", libpath)
}
