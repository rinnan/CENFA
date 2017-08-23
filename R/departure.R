#' Climatic departure
#'
#' Calculates the climatic departure of a species using historical and future climate raster data and species presence data.
#'
#' @param x.hist Raster* object, typically a brick or stack of historical climate raster layers
#' @param x.fut  Raster* object, future climate values with the same layers as x.hist
#' @param s.dat SpatialPolygons* object detailing species presence or abundance
#' @param field field of \code{speciesdat} that specifies presence or abundance. This is equivalent to the \code{field} argument of \code{raster::rasterize}.
#' @param nf integer. Specifies the number of specialization axes to keep after transformation.
#' @param scale logical. If \code{TRUE} then the values of \code{x.hist} and \code{x.fut} will
#' be centered and scaled by the means and sds of the historical climate data. Depending on the resolution of the data and the
#' extent of the study area, this can be quite time consuming. If running this
#' function for multiple species, it is recommended that the climate data be scaled beforehand.
#' @param sp.prj character. Spatial projection of species data.
#' @return Returns an S4 object of class \code{cnfa} with the following slots:
#' @return call original function call
#' @return departure climatic departure D of species
#' @return distances vector of distances d_ij
#' @return departure_ras raster of distances d_ij
#' @return present number of cells in which species is present
#'
#' @export
#'
#'

setGeneric("departure", function(x.hist, x.fut, s.dat, ...) {
  standardGeneric("departure")
})

setClass("departure", slots = list(call = "call", departure = "numeric", distances = "numeric", departure_ras = "Raster", present = "numeric"))

#' @rdname departure
setMethod("departure",
          signature(x.hist = "RasterBrick", x.fut = "RasterBrick", s.dat = "SpatialPolygonsDataFrame"),
          function(x.hist, x.fut, s.dat, field,
                   scale = FALSE, depart.ras = TRUE){
            call <- match.call()
            if(!identicalCRS(x.hist, s.dat)) {stop("historical climate and species projections do not match")}
            if(!identicalCRS(x.hist, x.fut))     {stop("historical and future climate projections do not match")}
            if(!identicalCRS(x.fut, s.dat)) {stop("future climate and species projections do not match")}
            if(length(raster::intersect(extent(x.hist), extent(s.dat))) == 0) {stop("climate and species data to not overlap")}
            if(scale == TRUE) {
              x.hist <- raster::scale(x.hist)
              means <- cellStats(x.hist, mean)
              sds <- cellStats(x.fut, sd)
              x.fut <- scale(x.fut, center = means, scale = sds)
            }
            #gpres <- which(!is.na(values(x.hist[[1]])))
            #dat<-values(climdat)[gpres,]
            speciesdat.ras <- rasterize(s.dat, x.hist, field = field)
            pres <- which(!is.na(values(speciesdat.ras)))
            Ns <- length(pres)
            #prb<-values(speciesdat.ras)[pres]
            pres.dat <- values(x.hist)[pres,]
            d_ij <- values(x.fut)[pres,] - values(x.hist)[pres,]
            d <- sqrt(rowSums(d_ij^2))
            D <- 1/(1.96*Ns) * sum(d)

            if(depart.ras == T){
              ras <- speciesdat.ras
              values(ras)[pres] <- d
            }
            else ras <- NA



            depart <- methods::new("departure", call = call, departure = D, distances = d, departure_ras = ras, present = Ns)
            return(depart)
          }
)

#' @rdname departure
setMethod("departure",
          signature(x.hist = "Raster", x.fut = "Raster", s.dat = "enfa"),
          function(x.hist, x.fut, s.dat, depart.ras = TRUE, scale = FALSE){

            call <- match.call()
            sp.ras <- raster(s.dat)
            if(!identicalCRS(x.hist, sp.ras)) {stop("historical climate and species projections do not match")}
            if(!identicalCRS(x.hist, x.fut))     {stop("historical and future climate projections do not match")}
            if(!identicalCRS(x.fut, sp.ras)) {stop("future climate and species projections do not match")}
            if(length(raster::intersect(extent(x.hist), extent(sp.ras)))==0) {stop("climate and species data to not overlap")}
            if(scale == TRUE) {
              center <- cellStats(x.hist, mean)
              sds <- cellStats(x.hist, sd)
              x.hist <- raster::scale(x.hist, center = center, scale = sds)
              x.fut <- raster::scale(x.fut, center = center, scale = sds)
            }
            x.hist <- crop(x.hist, sp.ras)
            x.fut <- crop(x.fut, sp.ras)

            pres <- which(!is.na(values(sp.ras[[1]])))
            small <- canProcessInMemory(x.hist, 8)
            if(small){
              z_ij <- values(x.hist)[pres,] %*% as.matrix(s.dat@co)
              f_ij <- values(x.fut)[pres,] %*% as.matrix(s.dat@co)
              d_ij <- (f_ij - z_ij)^2
              d <- sqrt(rowSums(d_ij))
              D <- 1/(1.96) * mean(d, na.rm = T)
              sp.ras[pres] <- d
            } else {
              x.hist <- calc(x.hist, fun = function(x) {x %*% as.matrix(s.dat@co)})
              x.fut <- calc(x.fut, fun = function(x) {x %*% as.matrix(s.dat@co)})
              d_ij <- (x.fut - x.hist)^2
              d <- calc(d_ij, fun = function(x) {sqrt(sum(x))})
              values(d)[!pres] <- NA
              D <- 1/(1.96) * cellStats(d, mean)
              sp.ras <- d
            }

            # if(depart.ras == T){
            #   ras <- sp.ras[[1]]
            #   values(ras)[pres] <- d
            # }
            # else ras <- NA

            depart <- methods::new("departure", call = call, departure = D, distances = d, departure_ras = sp.ras, present = length(pres))
            return(depart)
          }
)

#' @rdname departure
setMethod("departure",
          signature(x.hist = "Raster", x.fut = "Raster", s.dat = "cnfa"),
          function(x.hist, x.fut, s.dat, depart.ras = TRUE, scale = FALSE){

            call <- match.call()
            sp.ras <- raster(s.dat)
            if(!identicalCRS(x.hist, sp.ras)) {stop("historical climate and species projections do not match")}
            if(!identicalCRS(x.hist, x.fut))     {stop("historical and future climate projections do not match")}
            if(!identicalCRS(x.fut, sp.ras)) {stop("future climate and species projections do not match")}
            if(length(raster::intersect(extent(x.hist), extent(sp.ras)))==0) {stop("climate and species data to not overlap")}
            if(scale == TRUE) {
              center <- cellStats(x.hist, mean)
              sds <- cellStats(x.hist, sd)
              x.hist <- raster::scale(x.hist, center = center, scale = sds)
              x.fut <- raster::scale(x.fut, center = center, scale = sds)
            }
            x.hist <- crop(x.hist, sp.ras)
            x.fut <- crop(x.fut, sp.ras)

            pres <- which(!is.na(values(sp.ras[[1]])))
            small <- canProcessInMemory(x.hist, 5)
            if(small){
              z_ij <- values(x.hist)[pres,] %*% as.matrix(s.dat@co)
              f_ij <- values(x.fut)[pres,] %*% as.matrix(s.dat@co)
              d_ij <- (f_ij - z_ij)^2
              d <- sqrt(rowSums(d_ij))
              D <- 1/(1.96) * mean(d, na.rm = T)
              sp.ras[pres] <- d
            } else {
              x.hist <- calc(x.hist, fun = function(x) {x %*% as.matrix(s.dat@co)})
              x.fut <- calc(x.fut, fun = function(x) {x %*% as.matrix(s.dat@co)})
              d_ij <- (x.fut - x.hist)^2
              d <- calc(d_ij, fun = function(x) {sqrt(sum(x))})
              values(d)[!pres] <- NA
              D <- 1/(1.96) * cellStats(d, mean)
              distances <- values(d)[pres]
            }

            # if(depart.ras == T){
            #   ras <- sp.ras[[1]]
            #   values(ras)[pres] <- d
            # }
            # else ras <- NA

            depart <- methods::new("departure", call = call, departure = D, distances = d, departure_ras = d, present = length(pres))
            return(depart)
          }
)
