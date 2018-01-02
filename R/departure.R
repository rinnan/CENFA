#' Climatic departure
#'
#' Calculates the climatic departure of a species using historical and future
#' climate raster data and species presence data.
#'
#' @param x.hist Raster* object, typically a brick or stack of historical climate
#'   raster layers
#' @param x.fut  Raster* object, future climate values with the same layers as x.hist
#' @param s.dat SpatialPolygons*, sf, or cnfa object detailing species presence
#' @param field field of \code{speciesdat} that specifies presence. This is
#'   equivalent to the \code{field} argument of \code{raster::rasterize}.
#' @param fun function or character. Determines what values to assign to cells
#'   with multiple spatial features, similar to the \code{fun} argument in
#'   \code{\link[raster]{rasterize}}.
#' @param scale logical. If \code{TRUE} then the values of \code{x.hist} and
#'   \code{x.fut} will be centered and scaled by the means and sds of the historical
#'   climate data. Depending on the resolution of the data and the extent of the
#'   study area, this can be quite time consuming. If running this function for
#'   multiple species, it is recommended that the climate data be scaled beforehand.
#' @param parallel logical. If \code{TRUE} then multiple cores are utilized.
#' @param n numeric. Optional number of CPU cores to utilize for parallel processing.
#' @param ... Additional arguments for \code{\link[raster]{clusterR}}.
#'
#' @examples
#' dep <- departure(x.hist = climdat.hist, x.fut = climdat.fut, s.dat = ABPR)
#'
#' @return Returns an S4 object of class \code{departure} with the following slots:
#' \describe{
#'   \item{call}{Original function call.}
#'   \item{df}{Departure factor. Vector of length p that describes the amount of
#'    departure between future and historical conditions for each climate variable.}
#'   \item{departure}{Magnitude of the departure factor.}
#'   \item{present}{Number of raster cells in which species is present.}
#'   \item{ras}{RasterBrick of climate departures, with p layers.}
#' }
#'
#' @export
#' @importFrom stats sd
#' @importFrom parallel detectCores
#'

setGeneric("departure", function(x.hist, x.fut, s.dat, ...) {
  standardGeneric("departure")
})

#' @rdname departure
setMethod("departure",
          signature(x.hist = "Raster", x.fut = "missing", s.dat = "cnfa"),
          function(x.hist, s.dat, field, fun = "last", ...){

            call <- sys.calls()[[1]]
            sp.ras <- raster(s.dat)

            warning("were rasters scaled before calculating differences?")
            if(!identicalCRS(x.hist, sp.ras)) {stop("climate and species projections do not match")}
            if(length(intersect(extent(x.hist), extent(sp.ras))) == 0) {stop("climate and species data to not overlap")}

            nm <- names(x.hist)
            x.dif <- crop(x.hist, extent(sp.ras))
            x.dif <- mask(x.dif, sp.ras[[1]])
            names(x.dif) <- nm
            d <- cellStats(x.dif, mean)
            D <- norm(d, "2")
            pres <- which(!is.na(values(max(x.dif))))

            depart <- methods::new("departure", call = call, df = d, departure = D,  present = s.dat@present, ras = x.dif)
            return(depart)
          }
)

#' @rdname departure
setMethod("departure",
          signature(x.hist = "Raster", x.fut = "Raster", s.dat = "cnfa"),
          function(x.hist, x.fut, s.dat, field, fun = "last", scale = TRUE, parallel = FALSE, n, ...) {

            call <- sys.calls()[[1]]
            sp.ras <- raster(s.dat)

            if(names(x.hist) != names(x.fut)) stop("historical and future rasters do not match")
            if(!identicalCRS(x.hist, sp.ras)) stop("climate and species projections do not match")
            if(length(intersect(extent(x.hist), extent(sp.ras))) == 0) stop("climate and species data to not overlap")

            nm <- names(x.hist)
            x.hist <- crop(x.hist, extent(sp.ras))
            x.fut <- crop(x.fut, extent(sp.ras))

            if(parallel == T) {
              if(scale == T) {
                if (missing(n)) {
                  n <- parallel::detectCores()
                  message(n, ' cores detected, using ', n-1)
                  n <- n-1
                }
                means <- cellStats(x.hist, mean)
                sds <- cellStats(x.hist, sd)
                beginCluster(n = n)
                x.hist <- clusterR(x.hist, fun = raster::scale, export = list("means", "sds"), args = list(center = means, scale = sds))
                x.fut <- clusterR(x.fut, raster::scale, export = list("means", "sds"), args = list(center = means, scale = sds))
                endCluster()
              }
            } else {
              if(scale == T) {
                means <- cellStats(x.hist, mean)
                sds <- cellStats(x.hist, sd)
                x.hist <- raster::scale(x.hist)
                x.fut <- scale(x.fut, center = means, scale = sds)
              }
            }

            x.dif <- abs(x.fut - x.hist)
            x.dif <- mask(x.dif, sp.ras[[1]])
            names(x.dif) <- nm
            d <- cellStats(x.dif, mean)
            D <- norm(d, "2")
            pres <- which(!is.na(values(max(x.dif))))

            # if(canProcessInMemory(x.dif)) {
            #   pres <- which(!is.na(values(sp.ras[[1]])))
            #   d_ij <- values(x.dif)[pres,]
            #   d <- colMeans(d_ij)
            #   D <- norm(d, "2")
            #   distances <- apply(d_ij, 1, norm, "2")
            #   ras <- raster(sp.ras[[1]])
            #   ras[pres] <- distances
            #   names(ras) <- paste0("D", 1:nlayers(ras))
            # } else {
            #   x.mask.h <- mask(x.dif, sp.ras[[1]])
            #   d <- cellStats(x.mask.h, mean)
            #   D <- norm(d, "2")
            #   beginCluster(cores, exclude = "CENFA")
            #   f1 <- function(x) norm(x, "2")
            #   ras <- clusterR(x.mask.h, fun = f1, args = list(fun = f1), ...)
            #   names(ras) <- paste0("D.", 1:nlayers(sp.ras))
            #   endCluster()
            # }

            depart <- methods::new("departure", call = call, df = d, departure = D, present = s.dat@present, ras = x.dif)
            return(depart)
          }
)

#' @rdname departure
setMethod("departure",
          signature(x.hist = "Raster", x.fut = "Raster", s.dat = "SpatialPolygons"),
          function(x.hist, x.fut, s.dat, field, fun = "last", scale = TRUE, parallel = FALSE, n, ...){

            call <- sys.calls()[[1]]

            if (! inherits(s.dat, 'SpatialPolygons')) stop('"s.dat" should be a "SpatialPolygons*" object')
            if(!identicalCRS(x.hist, s.dat)) {stop("historical climate and species projections do not match")}
            if(!identicalCRS(x.hist, x.fut))     {stop("historical and future climate projections do not match")}
            if(!identicalCRS(x.fut, s.dat)) {stop("future climate and species projections do not match")}
            if(length(intersect(extent(x.hist), extent(s.dat))) == 0) {stop("climate and species data to not overlap")}

            nm <- names(x.hist)
            x.hist <- crop(x.hist, extent(s.dat))
            x.fut <- crop(x.fut, extent(s.dat))

            if(parallel == T) {
              if(scale == T) {
                if (missing(n)) {
                  n <- parallel::detectCores()
                  message(n, ' cores detected, using ', n-1)
                  n <- n-1
                }
                means <- cellStats(x.hist, mean)
                sds <- cellStats(x.hist, sd)
                beginCluster(n = n)
                x.hist <- clusterR(x.hist, fun = raster::scale, export = list("means", "sds"), args = list(center = means, scale = sds))
                x.fut <- clusterR(x.fut, raster::scale, export = list("means", "sds"), args = list(center = means, scale = sds))
                endCluster()
              }
            } else {
              if(scale == T) {
                means <- cellStats(x.hist, mean)
                sds <- cellStats(x.hist, sd)
                x.hist <- raster::scale(x.hist)
                x.fut <- scale(x.fut, center = means, scale = sds)
              }
            }

            x.dif <- abs(x.fut - x.hist)
            s.dat.ras <- rasterize(s.dat, raster(x.dif), field = field, fun = fun, ...)
            x.dif <- mask(x.dif, s.dat.ras)
            names(x.dif) <- nm
            d <- cellStats(x.dif, mean)
            D <- norm(d, "2")
            pres <- which(!is.na(values(max(x.dif))))

            depart <- methods::new("departure", call = call, df = d, departure = D, present = length(pres), ras = x.dif)
            return(depart)
          }
)


# departure2 <- function(x.hist, s.dat, scale = FALSE, cores = 1, ...){
#
#   call <- sys.calls()[[1]]
#   sp.ras <- raster(s.dat)
#   if(!identicalCRS(x.hist, sp.ras)) {stop("climate and species projections do not match")}
#   if(length(raster::intersect(extent(x.hist), extent(sp.ras)))==0) {stop("climate and species data to not overlap")}
#   if(scale) {
#     center <- cellStats(x.hist, mean)
#     sds <- cellStats(x.hist, sd)
#     x.hist <- raster::scale(x.hist, center = center, scale = sds)
#   }
#   x.hist <- crop(x.hist, extent(sp.ras))
#   names(x.hist) <- paste0("D", 1:nlayers(x.hist))
#   #Rs.inv <- solve(covmat(sp.ras, ...), tol = 1e-20)
#
#   small <- canProcessInMemory(x.hist, 3)
#   if(small){
#     pres <- which(!is.na(values(sp.ras[[1]])))
#     #d_ij <- values(x.hist)[pres,] %*% as.matrix(s.dat@co)
#     d_ij <- values(x.hist)[pres,]
#     d <- colMeans(d_ij)
#     D <- as.numeric(t(d) %*% d)
#     distances <- apply(d_ij, 1, norm, "2")
#     #D <- mean(distances, na.rm = T)
#     ras <- raster(sp.ras[[1]])
#     ras[pres] <- distances
#     names(ras) <- paste0("D", 1:nlayers(ras))
#   } else {
#     x.mask.h <- mask(x.hist, sp.ras[[1]])
#     d <- cellStats(x.mask.h, mean)
#     D <- as.numeric(t(d) %*% d)
#     beginCluster(cores, exclude = "CENFA")
#     #cl <- makeCluster(getOption("cl.cores", cores))
#     ##clusterExport(cl, c("x.mask.h", "calc", "Rs.inv", "s.dat@co", "s"), envir = environment())
#     #registerDoSNOW(cl)
#
#     f1 <- function(x) norm(x, "2")
#     ras <- clusterR(x.mask.h, fun = f1, args = list(fun = f1), ...)
#     names(ras) <- paste0("D", 1:nlayers(ras))
#     #f2 <- function(x) sqrt(t(x) %*% Rs.inv %*% x)
#     #ras <- clusterR(x.mask.h2, fun = calc, args = list(fun = f2), export = "Rs.inv", m = cores)
#     endCluster()
#   }
#   #x.mask.h <- calc(x.mask.h, function(x) x %*% as.matrix(s.dat@co))
#   #ras <- calc(x.mask.h, function(x) sqrt(t(x) %*% Rs.inv %*% x))
#   #D <- cellStats(ras, mean)
#
#   depart <- methods::new("departure", call = call, df = d, departure = D, departure_ras = ras, present = s.dat@present)
#   return(depart)
# }
#
