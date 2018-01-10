#' Ecological-Niche Factor Analysis
#'
#' Performs ecological-niche factor analysis using environmental raster data and
#' species presence data.
#'
#' @aliases print.enfa, show.enfa
#'
#' @param x Raster* object, typically a brick or stack of ecological raster layers
#' @param s.dat matrix, SpatialPolygons*, SpatialPoints*, or sf object indicating
#'   species presence
#' @param field field of \code{s.dat} that specifies presence or abundance. This
#'   is equivalent to the \code{field} argument in the \code{raster} package.
#' @param fun function or character. Determines what values to assign to cells
#'   with multiple spatial features, similar to the \code{fun} argument in
#'   \code{\link[raster]{rasterize}}.
#' @param scale logical. If \code{TRUE} then the values of the Raster* object will
#'   be centered and scaled. Depending on the resolution of the climate data and
#'   the extent of the study area, this can be quite time consuming. If running
#'   this function for multiple species, it is recommended that the climate data
#'   be scaled beforehand using the \code{\link{GLcenfa}} function.
#' @param filename character. Optional filename to save the Raster* output to
#'   file. If this is not provided, a temporary file will be created for large \code{x}.
#' @param parallel logical. If \code{TRUE} then multiple cores are utilized for the
#'   calculation of the covariance matrices.
#' @param prj character. proj4string denoting projection of species observation data.
#' @param ... Additonal arguments for the \code{\link{covmat}} function, such as the
#'   number of cores \code{n}.
#'
#' @details Hirzel et al. (2002) defined the overall marginality \eqn{M} as a
#'   standardized distance between the centroid of the species' niche and the
#'   global niche. Basille and Calenge (2008), however, defined \eqn{M} as
#'   \deqn{M = m^Tm.} The default \code{mar.type} reflects Basille and Calenge's
#'   definition so that results will by default agree with those calculated using
#'   Basille and Calenge's \code{adehabitatHS} package.
#'
#' @examples
#' mod1 <- enfa(x = climdat.hist, s.dat = ABPR, field = "CODE")
#'
#' # using GLcenfa as an initial step
#' # for multi-species comparison
#'
#' glc <- GLcenfa(x = climdat.hist)
#' mod2 <- enfa(x = glc, s.dat = ABPR, field = "CODE")
#' all.equal(m.factor(mod1), m.factor(mod2))
#'
#' @return Returns an S4 object of class \code{enfa} with the following components:
#' \describe{
#'   \item{call}{Original function call.}
#'   \item{mf}{Marginality factor. Vector that describes the location of the
#'    species Hutchinsonian niche relative to the global niche.}
#'   \item{marginality}{Magnitude of the marginality factor.}
#'   \item{sf}{Specialization factor. Vector of eigenvalues of specialization.}
#'   \item{specialization}{Overall specialization, equal to the square root of
#'    the sum of eigenvalues of specialization.}
#'   \item{s.prop}{Vector representing the amount of specialization found in each
#'    ENFA factor.}
#'   \item{co}{A matrix describing the amount of marginality and specialization
#'    on each ENFA factor. (The marginality column is normalized.)}
#'   \item{ras}{RasterBrick of transformed climate values, with p layers.}
#'   \item{weights}{Raster layer of weights used for ENFA calculation.}
#' }
#'
#' @seealso \code{\link{GLcenfa}}, \code{\link{cnfa}}
#'
#' @export

setGeneric("enfa", function(x, s.dat, ...) {
  standardGeneric("enfa")
})

#' @rdname enfa
setMethod("enfa",
          signature(x = "GLcenfa", s.dat = "Spatial"),
          function(x, s.dat, field, fun = "last", filename = "", parallel = F, ...){

            call <- sys.calls()[[1]]

            if (! inherits(s.dat, c('SpatialPolygons', 'SpatialPoints'))) stop('"s.dat" should be a "SpatialPolygons*" or "SpatialPoints*" object')
            if (!identicalCRS(raster(x), s.dat)) stop("climate and species projections do not match")
            ras <- raster(x)
            ext <- extent(ras)
            ext.s <- extent(s.dat)
            if (is.null(intersect(ext, ext.s))) stop("climate and species data do not overlap")
            if (union(ext, ext.s) != ext) stop("extent of species data not contained within extent of climate data")

            x.crop <- crop(ras, ext.s)
            s.dat.ras <- rasterize(s.dat, x.crop, field = field, fun = fun)

            filename <- trim(filename)
            if (!canProcessInMemory(x.crop) && filename == '') {
              filename <- rasterTmpFile()
            }

            if (canProcessInMemory(x.crop)){
              pres <- which(!is.na(values(s.dat.ras)) & !is.na(values(max(x.crop))))
              S <- values(x.crop)[pres,]
              nS <- nrow(S)
              Rg <- x@cov
              p <- values(s.dat.ras)[pres]
              p <- p/(sum(p))
              mar <- apply(S, 2, function(x) sum(x * p))
              Sm <- sweep(S, 2, mar)
              DpSm <- apply(Sm, 2, function(x) x * p)
              Rs <- crossprod(Sm, DpSm)
            } else {
              x.mask <- mask(x.crop, s.dat.ras)
              Rg <- x@cov
              p.sum <- cellStats(s.dat.ras, sum)
              p <- s.dat.ras / (p.sum - 1)
              DpS <- x.mask * p
              mar <- cellStats(DpS, sum)
              Sm <- calc(x.mask, fun = function(x) x - mar, forceapply = T)
              Rs <- covmat(x = Sm, w = p, parallel = parallel, ...)
            }

            cZ <- nlayers(ras)
            m <- sqrt(c(t(mar) %*% mar))
            if(max(Im(eigen(Rs)$values)) > 1e-05) stop("complex eigenvalues. Try removing correlated variables.")
            eigRs <- lapply(eigen(Rs), Re)
            keep <- (eigRs$values > 1e-09)
            Rs12 <- eigRs$vectors[, keep] %*% diag(eigRs$values[keep]^(-0.5)) %*% t(eigRs$vectors[, keep])
            W <- Rs12 %*% Rg %*% Rs12
            z <- Rs12 %*% mar
            y <- z/sqrt(sum(z^2))
            H <- (diag(cZ) - y %*% t(y)) %*% W %*% (diag(cZ) - y %*% t(y))
            sf <- eigen(H)$values[-cZ]
            spec <- sqrt(sum(sf))/length(sf)
            s.p <- abs(sum(diag(W)) - sum(diag(H)))
            s <- c(s.p, sf)
            s.p <- abs(s)/sum(abs(s))
            v <- eigen(H)$vectors
            co <- matrix(nrow = cZ, ncol = cZ)
            co[, 1] <- mar/sqrt(t(mar) %*% mar)
            u <- as.matrix((Rs12 %*% v)[, 1:(cZ-1)])
            norw <- sqrt(diag(t(as.matrix(u)) %*% as.matrix(u)))
            co[, -1] <- sweep(as.matrix(u), 2, norw, "/")
            nm <- c("Marg", paste0("Spec", (1:(cZ-1))))
            if(canProcessInMemory(x.crop)){
              s.ras <- brick(x.crop)
              values(s.ras)[pres, ] <- S %*% co
              names(s.ras) <- nm
            } else{
              cat("\nCreating factor rasters...")
              s.ras <- .calc(x.mask, function(x) {x %*% co}, forceapply = T, filename = filename, names = nm, ...)
            }
            colnames(co) <- names(s.p) <- nm
            rownames(co) <- names(x)
            names(sf) <- nm[-1]

            enfa <- methods::new("enfa", call = call, mf = mar, marginality = m, sf = sf,
                                 specialization = spec, p.spec = s.p, co = co, cov = Rs, ras = s.ras, weights = s.dat.ras)
            return(enfa)
          }
)

#' @rdname enfa
setMethod("enfa",
          signature(x = "Raster", s.dat = "Spatial"),
          function(x, s.dat, field, fun = "last", scale = TRUE, filename = "", parallel = F, ...){

            call <- sys.calls()[[1]]

            if (! inherits(x, 'Raster')) stop('"x" should be a "Raster*" object')
            if (! inherits(s.dat, c('SpatialPolygons', 'SpatialPoints'))) stop('"s.dat" should be a "SpatialPolygons*" or "SpatialPoints*" object')
            if(!identicalCRS(x, s.dat)) stop("projections do not match")
            if(is.null(intersect(extent(x), extent(s.dat)))) stop("climate and species data do not overlap")
            if(union(extent(x), extent(s.dat)) != extent(x)) stop("extent of species data not contained within extent of climate data")

            s.dat.ras <- rasterize(s.dat, raster(x), field = field, fun = fun)

            if(scale == TRUE) x <- scale(x)

            filename <- trim(filename)
            if (!canProcessInMemory(x) && filename == '') {
              filename <- rasterTmpFile()
            }

            if(canProcessInMemory(x)){
              gpres <- which(!is.na(values(max(x))))
              pres <- which(!is.na(values(s.dat.ras)) & !is.na(values(max(x))))
              Z <- values(x)[gpres, ]
              S <- values(x)[pres, ]
              nZ <- nrow(Z)
              nS <- nrow(S)
              Rg <- crossprod(Z, Z)/nZ
              p <- values(s.dat.ras)[pres]
              p <- p/(sum(p))
              mar <- apply(S, 2, function(x) sum(x * p))
              Sm <- sweep(S, 2, mar)
              DpSm <- apply(Sm, 2, function(x) x * p)
              Rs <- crossprod(Sm, DpSm)
            } else {
              center <- cellStats(x, mean)
              x.mask <- mask(x, s.dat.ras)
              p.sum <- cellStats(s.dat.ras, sum)
              p <- s.dat.ras / p.sum
              DpS <- x.mask * p
              mar <- cellStats(DpS, sum)
              cat("\nCalculating study area covariance matrix...\n")
              Rg <- covmat(x, sample = F, parallel = parallel, ...)
              cat("\nCalculating species covariance matrix...\n")
              Sm <- calc(x.mask, fun = function(x) x - mar, forceapply = T)
              Rs <- covmat(x = Sm, w = p, ...)
            }

            cZ <- nlayers(x)
            m <- sqrt(c(t(mar) %*% mar))
            if(max(Im(eigen(Rs)$values)) > 1e-05) stop("complex eigenvalues. Try removing correlated variables.")
            eigRs <- lapply(eigen(Rs), Re)
            keep <- (eigRs$values > 1e-09)
            Rs12 <- eigRs$vectors[, keep] %*% diag(eigRs$values[keep]^(-0.5)) %*% t(eigRs$vectors[, keep])
            W <- Rs12 %*% Rg %*% Rs12
            z <- Rs12 %*% mar
            y <- z/sqrt(sum(z^2))
            H <- (diag(cZ) - y %*% t(y)) %*% W %*% (diag(cZ) - y %*% t(y))
            sf <- eigen(H)$values[-cZ]
            spec <- sqrt(sum(sf))/length(sf)
            s.p <- abs(sum(diag(W)) - sum(diag(H)))
            s <- c(s.p, sf)
            s.p <- abs(s)/sum(abs(s))
            v <- eigen(H)$vectors
            co <- matrix(nrow = cZ, ncol = cZ)
            co[, 1] <- mar/sqrt(t(mar) %*% mar)
            u <- as.matrix((Rs12 %*% v)[, 1:(cZ-1)])
            norw <- sqrt(diag(t(as.matrix(u)) %*% as.matrix(u)))
            co[, -1] <- sweep(as.matrix(u), 2, norw, "/")
            nm <- c("Marg", paste0("Spec", (1:(cZ-1))))
            if(canProcessInMemory(x)){
              s.ras <- brick(x)
              values(s.ras)[pres, ] <- S %*% co
              names(s.ras) <- nm
            } else{
              cat("\nCreating factor rasters...")
              s.ras <- .calc(x.mask, function(x) {x %*% co}, forceapply = T, filename = filename, names = nm, ...)
            }
            colnames(co) <- names(s.p) <- nm
            rownames(co) <- names(x)
            names(sf) <- nm[-1]

            enfa <- methods::new("enfa", call = call, mf = mar, marginality = m, sf = sf,
                                 specialization = spec, p.spec = s.p, co = co, cov = Rs, ras = s.ras, weights = s.dat.ras)
            return(enfa)
          }
)

#' @rdname enfa
setMethod("enfa",
          signature(x = "Raster", s.dat = "sf"),
          function(x, s.dat, field, fun = "last", scale = TRUE, filename = "", parallel = F, ...){
            if (!requireNamespace("sf")) {
              warning('cannot do this because sf is not available')
            }

            call <- sys.calls()[[1]]

            if (! inherits(x, 'Raster')) stop('"x" should be a "Raster*" object')
            if(!identicalCRS(x, s.dat)) stop("projections do not match")
            if(is.null(intersect(extent(x), extent(s.dat)))) stop("climate and species data do not overlap")
            if(union(extent(x), extent(s.dat)) != extent(x)) stop("extent of species data not contained within extent of climate data")

            s.dat <- as(s.dat, "Spatial")
            if (! inherits(s.dat, c('SpatialPolygons', 'SpatialPoints'))) stop('geometry of "s.dat" should be of class "sfc_POLYGON", "sfc_MULTIPOLYGON", "sfc_POINT", or "sfc_MULTIPOINT"')

            s.dat.ras <- rasterize(s.dat, raster(x), field = field, fun = fun)

            if(scale == TRUE) x <- scale(x)

            filename <- trim(filename)
            if (!canProcessInMemory(x) && filename == '') {
              filename <- rasterTmpFile()
            }

            if(canProcessInMemory(x)){
              gpres <- which(!is.na(values(max(x))))
              pres <- which(!is.na(values(s.dat.ras)) & !is.na(values(max(x))))
              Z <- values(x)[gpres, ]
              S <- values(x)[pres, ]
              nZ <- nrow(Z)
              nS <- nrow(S)
              Rg <- crossprod(Z, Z)/nZ
              p <- values(s.dat.ras)[pres]
              p <- p/(sum(p))
              mar <- apply(S, 2, function(x) sum(x * p))
              Sm <- sweep(S, 2, mar)
              DpSm <- apply(Sm, 2, function(x) x * p)
              Rs <- crossprod(Sm, DpSm)
            } else {
              center <- cellStats(x, mean)
              x.mask <- mask(x, s.dat.ras)
              p.sum <- cellStats(s.dat.ras, sum)
              p <- s.dat.ras / p.sum
              DpS <- x.mask * p
              mar <- cellStats(DpS, sum)
              cat("\nCalculating study area covariance matrix...\n")
              Rg <- covmat(x, sample = F, ...)
              cat("\nCalculating species covariance matrix...\n")
              Sm <- calc(x.mask, fun = function(x) x - mar, forceapply = T)
              Rs <- covmat(x = Sm, w = p, parallel = parallel, ...)
            }

            cZ <- nlayers(x)
            m <- sqrt(c(t(mar) %*% mar))
            if(max(Im(eigen(Rs)$values)) > 1e-05) stop("complex eigenvalues. Try removing correlated variables.")
            eigRs <- lapply(eigen(Rs), Re)
            keep <- (eigRs$values > 1e-09)
            Rs12 <- eigRs$vectors[, keep] %*% diag(eigRs$values[keep]^(-0.5)) %*% t(eigRs$vectors[, keep])
            W <- Rs12 %*% Rg %*% Rs12
            z <- Rs12 %*% mar
            y <- z/sqrt(sum(z^2))
            H <- (diag(cZ) - y %*% t(y)) %*% W %*% (diag(cZ) - y %*% t(y))
            sf <- eigen(H)$values[-cZ]
            spec <- sqrt(sum(sf))/length(sf)
            s.p <- abs(sum(diag(W)) - sum(diag(H)))
            s <- c(s.p, sf)
            s.p <- abs(s)/sum(abs(s))
            v <- eigen(H)$vectors
            co <- matrix(nrow = cZ, ncol = cZ)
            co[, 1] <- mar/sqrt(t(mar) %*% mar)
            u <- as.matrix((Rs12 %*% v)[, 1:(cZ-1)])
            norw <- sqrt(diag(t(as.matrix(u)) %*% as.matrix(u)))
            co[, -1] <- sweep(as.matrix(u), 2, norw, "/")
            nm <- c("Marg", paste0("Spec", (1:(cZ-1))))
            if(canProcessInMemory(x)){
              s.ras <- brick(x)
              values(s.ras)[pres, ] <- S %*% co
              names(s.ras) <- nm
            } else{
              cat("\nCreating factor rasters...")
              s.ras <- .calc(x.mask, function(x) {x %*% co}, forceapply = T, filename = filename, names = nm, ...)
            }
            colnames(co) <- names(s.p) <- nm
            rownames(co) <- names(x)
            names(sf) <- nm[-1]

            enfa <- methods::new("enfa", call = call, mf = mar, marginality = m, sf = sf,
                                 specialization = spec, p.spec = s.p, co = co, cov = Rs, ras = s.ras, weights = s.dat.ras)
            return(enfa)
          }
)

#' @rdname enfa
setMethod("enfa",
          signature(x = "Raster", s.dat = "matrix"),
          function(x, s.dat, prj, field, fun = "last", scale = TRUE, filename = "", parallel = F, ...){

            call <- match.call()

            if(!inherits(x, 'Raster')) stop('"x" should be a "Raster*" object')
            if(missing(prj)) {
              prj <- crs(x)
              warning("prj argument not supplied: inheriting projection from Raster* object.")
            }
            s.dat <- SpatialPointsDataFrame(coords = s.dat, data = 1, proj4string = crs(prj))
            if(!identicalCRS(x, s.dat)) stop("projections do not match")

            if(is.null(intersect(extent(x), extent(s.dat)))) stop("climate and species data do not overlap")

            s.dat.ras <- rasterize(s.dat, raster(x), field = field, fun = fun, ...)

            if(scale == TRUE) x <- scale(x)

            filename <- trim(filename)
            if (!canProcessInMemory(x) && filename == '') {
              filename <- rasterTmpFile()
            }

            if(canProcessInMemory(x)){
              gpres <- which(!is.na(values(max(x))))
              pres <- which(!is.na(values(s.dat.ras)) & !is.na(values(max(x))))
              Z <- values(x)[gpres, ]
              S <- values(x)[pres, ]
              nZ <- nrow(Z)
              nS <- nrow(S)
              Rg <- crossprod(Z, Z)/nZ
              p <- values(s.dat.ras)[pres]
              p <- p/(sum(p))
              mar <- apply(S, 2, function(x) sum(x * p))
              Sm <- sweep(S, 2, mar)
              DpSm <- apply(Sm, 2, function(x) x * p)
              Rs <- crossprod(Sm, DpSm)
            } else {
              center <- cellStats(x, mean)
              x.mask <- mask(x, s.dat.ras)
              #pres <- which(!is.na(values(s.dat.ras)) & !is.na(values(max(x.mask))))
              p.sum <- cellStats(s.dat.ras, sum)
              p <- s.dat.ras / p.sum
              DpS <- x.mask * p
              mar <- cellStats(DpS, sum)
              cat("\nCalculating study area covariance matrix...\n")
              Rg <- covmat(x, sample = F, ...)
              cat("\nCalculating species covariance matrix...\n")
              Sm <- calc(x.mask, fun = function(x) x - mar, forceapply = T)
              Rs <- covmat(x = Sm, w = p, parallel = parallel, ...)
              #Rs <- layerStats(Sm, 'weighted.cov', w = p, na.rm = T, asSample = F)[[1]]
              #Rs <- covmat(x.mask, center = T)
              #mar <- cellStats(x.mask, sum)/length(pres)
            }

            cZ <- nlayers(x)
            m <- sqrt(c(t(mar) %*% mar))
            if(max(Im(eigen(Rs)$values)) > 1e-05) stop("complex eigenvalues. Try removing correlated variables.")
            eigRs <- lapply(eigen(Rs), Re)
            keep <- (eigRs$values > 1e-09)
            Rs12 <- eigRs$vectors[, keep] %*% diag(eigRs$values[keep]^(-0.5)) %*% t(eigRs$vectors[, keep])
            W <- Rs12 %*% Rg %*% Rs12
            z <- Rs12 %*% mar
            y <- z/sqrt(sum(z^2))
            H <- (diag(cZ) - y %*% t(y)) %*% W %*% (diag(cZ) - y %*% t(y))
            sf <- eigen(H)$values[-cZ]
            spec <- sqrt(sum(sf))/length(sf)
            s.p <- abs(sum(diag(W)) - sum(diag(H)))
            s <- c(s.p, sf)
            s.p <- abs(s)/sum(abs(s))
            v <- eigen(H)$vectors
            co <- matrix(nrow = cZ, ncol = cZ)
            co[, 1] <- mar/sqrt(t(mar) %*% mar)
            u <- as.matrix((Rs12 %*% v)[, 1:(cZ-1)])
            norw <- sqrt(diag(t(as.matrix(u)) %*% as.matrix(u)))
            co[, -1] <- sweep(as.matrix(u), 2, norw, "/")
            nm <- c("Marg", paste0("Spec", (1:(cZ-1))))
            if(canProcessInMemory(x)){
              s.ras <- brick(x)
              values(s.ras)[pres, ] <- S %*% co
              names(s.ras) <- nm
            } else{
              cat("\nCreating factor rasters...")
              s.ras <- .calc(x.mask, function(x) {x %*% co}, forceapply = T, filename = filename, names = nm, ...)
            }
            colnames(co) <- names(s.p) <- nm
            rownames(co) <- names(x)
            names(sf) <- nm[-1]

            enfa <- methods::new("enfa", call = call, mf = mar, marginality = m, sf = sf,
                                 specialization = spec, p.spec = s.p, co = co, cov = Rs, ras = s.ras, weights = s.dat.ras)
            return(enfa)
          }
)
