#' difRasterBrick-class
#'
#' An object of class \code{difRasterBrick} is created by the \code{\link{difRaster}}
#' function, and inherits from the \code{RasterBrick} class.
#'
#' @slot title Character
#' @slot file Object of class ".RasterFile"
#' @slot data Object of class ".SingleLayerData" or ".MultipleLayerData"
#' @slot history To record processing history, not yet in use
#' @slot legend Object of class .RasterLegend, Default legend. Should store
#'   preferences for plotting. Not yet implemented except that it stores the color
#'   table of images, if available
#' @slot extent Object of Extent-class
#' @slot ncols Integer
#' @slot nrows Integer
#' @slot crs Object of class "CRS", i.e. the coordinate reference system. In
#'   Spatial* objects this slot is called 'proj4string'
#'
#' @export

setClass("difRasterBrick", contains = "RasterBrick")
