#' Climatic departure
#'
#' Calculates the climatic departure of a species using historical and future climate raster data and species presence data.
#'
#' @aliases
#' @param hist.dat Raster* object, typically a brick or stack of historical climate raster layers
#' @param fut.dat  Raster* object, future climate values with the same layers as hist.dat
#' @param species.dat SpatialPolygons* object detailing species presence or abundance
#' @param field field of \code{speciesdat} that specifies presence or abundance. This is equivalent to the \code{field} argument of \code{raster::rasterize}.
#' @param nf integer. Specifies the number of specialization axes to keep after transformation.
#' @param scale logical. If \code{TRUE} then the values of \code{hist.dat} and \code{fut.dat} will
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

setGeneric("departure", function(hist.dat, fut.dat, species.dat, ...) {
  standardGeneric("departure")
})

setClass("departure", slots = list(call = "call", departure = "numeric", distances = "numeric", departure_ras = "Raster", present = "numeric"))

#' @rdname departure
setMethod("departure",
          signature(hist.dat = "RasterBrick", fut.dat = "RasterBrick", species.dat = "SpatialPolygonsDataFrame"),
          function(hist.dat, fut.dat, species.dat, field,
                   scale = FALSE, depart.ras = TRUE){
            call <- match.call()
            if(!identicalCRS(hist.dat, species.dat)) {stop("historical climate and species projections do not match")}
            if(!identicalCRS(hist.dat, fut.dat))     {stop("historical and future climate projections do not match")}
            if(!identicalCRS(fut.dat, species.dat)) {stop("future climate and species projections do not match")}
            if(length(raster::intersect(extent(hist.dat), extent(species.dat)))==0) {stop("climate and species data to not overlap")}
            if(scale == TRUE) {
              hist.dat <- raster::scale(hist.dat)
              means <- cellStats(hist.dat, mean)
              sds <- cellStats(fut.dat, sd)
              fut.dat <- scale(fut.dat, center = means, scale = sds)
            }
            #gpres <- which(!is.na(values(hist.dat[[1]])))
            #dat<-values(climdat)[gpres,]
            speciesdat.ras <- rasterize(species.dat, hist.dat, field = field)
            pres <- which(!is.na(values(speciesdat.ras)))
            Ns <- length(pres)
            #prb<-values(speciesdat.ras)[pres]
            pres.dat <- values(hist.dat)[pres,]
            d_ij <- values(fut.dat)[pres,] - values(hist.dat)[pres,]
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
          signature(hist.dat = "Raster", fut.dat = "Raster", species.dat = "enfa"),
          function(hist.dat, fut.dat, species.dat, depart.ras = TRUE, scale = FALSE){

            call <- match.call()
            sp.ras <- raster(species.dat)
            if(!identicalCRS(hist.dat, sp.ras)) {stop("historical climate and species projections do not match")}
            if(!identicalCRS(hist.dat, fut.dat))     {stop("historical and future climate projections do not match")}
            if(!identicalCRS(fut.dat, sp.ras)) {stop("future climate and species projections do not match")}
            if(length(raster::intersect(extent(hist.dat), extent(sp.ras)))==0) {stop("climate and species data to not overlap")}
            if(scale == TRUE) {
              center <- cellStats(hist.dat, mean)
              sds <- cellStats(hist.dat, sd)
              hist.dat <- raster::scale(hist.dat, center = center, scale = sds)
              fut.dat <- raster::scale(fut.dat, center = center, scale = sds)
            }
            hist.dat <- crop(hist.dat, sp.ras)
            fut.dat <- crop(fut.dat, sp.ras)

            pres <- which(!is.na(values(sp.ras[[1]])))
            small <- canProcessInMemory(hist.dat, 8)
            if(small){
              z_ij <- values(hist.dat)[pres,] %*% as.matrix(species.dat@co)
              f_ij <- values(fut.dat)[pres,] %*% as.matrix(species.dat@co)
              d_ij <- (f_ij - z_ij)^2
              d <- sqrt(rowSums(d_ij))
              D <- 1/(1.96) * mean(d, na.rm = T)
              sp.ras[pres] <- d
            } else {
              hist.dat <- calc(hist.dat, fun = function(x) {x %*% as.matrix(species.dat@co)})
              fut.dat <- calc(fut.dat, fun = function(x) {x %*% as.matrix(species.dat@co)})
              d_ij <- (fut.dat - hist.dat)^2
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
          signature(hist.dat = "Raster", fut.dat = "Raster", species.dat = "cnfa"),
          function(hist.dat, fut.dat, species.dat, depart.ras = TRUE, scale = FALSE){

            call <- match.call()
            sp.ras <- raster(species.dat)
            if(!identicalCRS(hist.dat, sp.ras)) {stop("historical climate and species projections do not match")}
            if(!identicalCRS(hist.dat, fut.dat))     {stop("historical and future climate projections do not match")}
            if(!identicalCRS(fut.dat, sp.ras)) {stop("future climate and species projections do not match")}
            if(length(raster::intersect(extent(hist.dat), extent(sp.ras)))==0) {stop("climate and species data to not overlap")}
            if(scale == TRUE) {
              center <- cellStats(hist.dat, mean)
              sds <- cellStats(hist.dat, sd)
              hist.dat <- raster::scale(hist.dat, center = center, scale = sds)
              fut.dat <- raster::scale(fut.dat, center = center, scale = sds)
            }
            hist.dat <- crop(hist.dat, sp.ras)
            fut.dat <- crop(fut.dat, sp.ras)

            pres <- which(!is.na(values(sp.ras[[1]])))
            small <- canProcessInMemory(hist.dat, 5)
            if(small){
              z_ij <- values(hist.dat)[pres,] %*% as.matrix(species.dat@co)
              f_ij <- values(fut.dat)[pres,] %*% as.matrix(species.dat@co)
              d_ij <- (f_ij - z_ij)^2
              d <- sqrt(rowSums(d_ij))
              D <- 1/(1.96) * mean(d, na.rm = T)
              sp.ras[pres] <- d
            } else {
              hist.dat <- calc(hist.dat, fun = function(x) {x %*% as.matrix(species.dat@co)})
              fut.dat <- calc(fut.dat, fun = function(x) {x %*% as.matrix(species.dat@co)})
              d_ij <- (fut.dat - hist.dat)^2
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
