#' Create a sensitivity map
#'
#' Creates a sensitivity map of species habitat.
#'
#' @param cnfa cnfa object.
#' @param cores numeric. Number of cores to use for calculation (optional).
#' @param filename character. Output filename (optional).
#' @param ... Additional arguments for file writing as for \code{\link[raster]{writeRaster}}.
#' @details The values of the sensitivity raster are calculated by centering the habitat's climate data around \strong{m} and projecting onto the sensitivity factor \strong{s}, given by formula  |\bold{S} - \bold{m}|\bold{s}.
#' @return RasterLayer.
#' @export

sensitivity_map <- function(cnfa, cores = 1, filename = "", ...){

  ras <- cnfa@ras
  m <- cnfa@mf
  s <- cnfa@sf

  filename <- trim(filename)
  if (!canProcessInMemory(ras) && filename == '') {
    filename <- rasterTmpFile()
  }

  f1 <- function(x) abs(x - m) %*% s

  if(cores == 1){
    sens.ras <- .calc(ras, fun = f1, forceapply = T, filename = filename, names = "Sensitivity", ...)
  }

  if(cores > 1){
    beginCluster(cores, exclude = "CENFA")
    sens.ras <- clusterR(ras, fun = .calc, args = list(fun = f1, forceapply = T, names = "Sensitivity"), filename = filename, ...)
    endCluster()
  }

  return(sens.ras)
}
