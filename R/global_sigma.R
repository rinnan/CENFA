#' Climate-Niche Factor Analysis
#'
#' Performs climate-niche factor analysis using climate raster data and species presence data in shapefile.
#'
#' @param climdat Raster* object, typically a brick or stack of climate raster layers
#' @export

global_sigma<-function(climdat){
  sigs<-cellStats(climdat,'sd')
  sqrt(sum(sigs^2))
}
