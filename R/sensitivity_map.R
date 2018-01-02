#' Create a sensitivity map
#'
#' Creates a sensitivity map of species habitat.
#'
#' @param cnfa Object of class \code{cnfa}.
#' @param parallel logical. If \code{TRUE} then multiple cores are utilized.
#' @param n numeric. Number of cores to use for calculation (optional).
#' @param filename character. Output filename (optional).
#' @param ... Additional arguments for file writing as for \code{\link[raster]{writeRaster}}.
#'
#' @details The values of the sensitivity raster are calculated by centering the
#' habitat's climate data around \strong{m} and projecting onto the sensitivity
#' factor \strong{s}, given by formula  |\bold{S} - \bold{m}|\bold{s}.
#'
#' @examples
#' mod1 <- cnfa(x = climdat.hist, s.dat = ABPR, field = "CODE")
#' sens.map <- sensitivity_map(mod1)
#'
#' @return RasterLayer of sensitivity values.
#'
#' @seealso \code{\link{cnfa}}
#'
#' @export

sensitivity_map <- function(cnfa, parallel = FALSE, n, filename = "", ...){

  ras <- cnfa@ras
  m <- cnfa@mf
  s <- cnfa@sf

  filename <- trim(filename)
  if (!canProcessInMemory(ras) && filename == '') {
    filename <- rasterTmpFile()
  }

  f1 <- function(x) abs(x - m) %*% s

  if(parallel) {
    beginCluster(n, exclude = "CENFA")
    sens.ras <- clusterR(ras, fun = .calc, args = list(fun = f1, forceapply = T, names = "Sensitivity"), filename = filename, ...)
    endCluster()
  } else {
    sens.ras <- .calc(ras, fun = f1, forceapply = T, filename = filename, names = "Sensitivity", ...)
  }

  return(sens.ras)
}
