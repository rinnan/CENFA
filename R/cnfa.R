#' Climate-Niche Factor Analysis
#'
#' Performs climate-niche factor analysis using climate raster data and species
#' presence data.
#'
#' @aliases print.cnfa, show.cnfa
#'
#' @param x Raster* object, typically a brick or stack with p climate
#'   raster layers, or a \code{GLcenfa} object
#' @param s.dat SpatialPolygons*, SpatialPoints*, or sf object indicating
#'   species presence
#' @param field field of \code{s.dat} that specifies presence or abundance. This
#'   is equivalent to the \code{field} argument in \code{\link[raster]{rasterize}}
#' @param fun function or character. Determines what values to assign to cells
#'   with multiple spatial features, similar to the \code{fun} argument in
#'   \code{\link[raster]{rasterize}}.  Options are 'first', 'last' (default),
#'   and 'count' (see Details)
#' @param scale logical. If \code{TRUE} then the values of \code{x} will get
#'   centered and scaled. Depending on the resolution of the climate data and
#'   the extent of the study area, this can be quite time consuming. If running
#'   this function for multiple species, it is recommended that the data be
#'   scaled beforehand using the \code{\link{GLcenfa}} function
#' @param filename character. Optional filename to save the Raster* output to
#'   file. If this is not provided, a temporary file will be created for large \code{x}
#' @param quiet logical. If \code{TRUE}, messages and progress bar will be
#'   suppressed
#' @param parallel logical. If \code{TRUE} then multiple cores are utilized for the
#'   calculation of the covariance matrices
#' @param n numeric. Number of CPU cores to utilize for parallel processing
#' @param ... Additional arguments for \code{\link[raster]{writeRaster}}
#'
#' @examples
#' mod1 <- cnfa(x = climdat.hist, s.dat = ABPR, field = "CODE")
#'
#' # using GLcenfa as an initial step
#' # for multi-species comparison
#'
#' glc <- GLcenfa(x = climdat.hist)
#' mod2 <- cnfa(x = glc, s.dat = ABPR, field = "CODE")
#'
#'# same results either way
#' all.equal(m.factor(mod1), m.factor(mod2))
#' all.equal(s.factor(mod1), s.factor(mod2))
#'
#' @return Returns an S4 object of class \code{cnfa} with the following components:
#' \describe{
#'   \item{call}{Original function call}
#'   \item{mf}{Marginality factor. Vector of length p that describes the location
#'   of the species Hutchinsonian niche relative to the global niche}
#'   \item{marginality}{Magnitude of the marginality factor, scaled by the
#'   global covariance matrix}
#'   \item{sf}{Sensitivity factor. Vector of length p that describes the amount of
#'    sensitivity for each climate variable}
#'   \item{sensitivity}{Magnitude of the sensitivity factor, scaled by the
#'   global covariance matrix}
#'   \item{s.prop}{Vector of length p representing the amount of specialization
#'   found in each CNFA factor}
#'   \item{co}{A p x p matrix describing the amount of marginality and specialization
#'    on each CNFA factor. (The marginality column is normalized)}
#'   \item{ras}{RasterBrick of transformed climate values, with p layers}
#'   \item{weights}{Raster layer of weights used for CNFA calculation}
#' }
#'
#' @details
#' The \code{cnfa} function is not to be confused with the
#' \code{\link{enfa}} function. \code{enfa} performs ENFA as described by Hirzel
#' et al. (2002) and Basille et al. (2008), and is offered as an alternative to
#' the \code{enfa} function in the \code{adehabitatHS} package. There are
#' several key differences between ENFA and CNFA.
#'
#' Whereas ENFA returns a \strong{specialization factor} that describes
#' the specialization in each \strong{ENFA factor}, CNFA returns a
#' \strong{sensitivity factor} \code{sf} that describes the sensitivity in each
#' \strong{environmental variable}. This makes the sensitivity factor more
#' directly comparable to the marginality factor \code{mf}, because their
#' dimensions are identical. Sensitivity is calculated by a weighted sum
#' of the amount of specialization found in each CNFA factor, \emph{including}
#' the marginality factor. As such, the sensitivity factor offers a more complete
#' measure of specialization than ENFA's specialization factor, which does
#' not calculate the amount of specialization found in the marginality factor.
#' As such, CNFA's overall sensitivity (found in the slot \code{sensitivity})
#' is likely more meaningful than ENFA's overall specialization (found in the
#' slot \code{specialization}).
#'
#' The default \code{fun = 'last'} gives equal weight to each occupied cell.
#' If multiple species observations occur in the same cell, the cell will only
#' be counted once. \code{fun = 'count'} will weight the cells by the number
#' of observations.
#'
#' If there is too much correlation between the layers of \code{x}, the global
#' covariance matrix will be singular, and the overall marginality and overall
#' sensitivity will not be meaningful. In this case, a warning is issued,
#' and \code{marginality} and \code{sensitivity} are both returned as \code{NA}.
#'
#' @references
#' Basille, Mathieu, et al. Assessing habitat selection using multivariate
#' statistics: Some refinements of the ecological-niche factor analysis. Ecological
#' Modelling 211.1 (2008): 233-240.
#'
#' Hirzel, Alexandre H., et al. Ecological-niche factor analysis: how to compute
#' habitat-suitability maps without absence data?. Ecology 83.7 (2002): 2027-2036.
#'
#' @seealso \code{\link{GLcenfa}}, \code{\link{enfa}}
#'
#' @export
#'
#' @importFrom stats cov
#' @importFrom methods as

setGeneric("cnfa", function(x, s.dat, ...){
  standardGeneric("cnfa")})

#' @rdname cnfa
setMethod("cnfa",
          signature(x = "GLcenfa", s.dat = "Spatial"),
          function(x, s.dat, field, fun = "last", filename = "", quiet = TRUE, parallel = FALSE, n = 1, ...){

            call <- sys.calls()[[1]]

            if (! inherits(s.dat, c('SpatialPolygons', 'SpatialPoints'))) stop('"s.dat" should be a "SpatialPolygons*" or "SpatialPoints*" object')
            if (!identicalCRS(raster(x), s.dat)) stop("climate and species projections do not match")
            ras <- raster(x)
            ext <- extent(ras)
            ext.s <- extent(s.dat)
            if (is.null(intersect(ext, ext.s))) stop("climate and species data do not overlap")
            if (raster::union(ext, ext.s) != ext) stop("extent of species data not contained within extent of climate data")
            if(parallel){
              if(n < 1 | !is.numeric(n)) {
                n <- parallel::detectCores() - 1
                message('incorrect number of cores specified, using ', n)
              } else if(n > parallel::detectCores()) {
                n <- parallel::detectCores() - 1
                message('too many cores specified, using ', n)
              }
            }

            x.crop <- crop(ras, ext.s)
            s.dat.ras <- rasterize(s.dat, x.crop, field = field, fun = fun)

            filename <- trim(filename)
            if (!canProcessInMemory(x.crop) && filename == '') {
              filename <- rasterTmpFile()
            }

            if (canProcessInMemory(x.crop) & !parallel) {
              pres <- which(!is.na(values(s.dat.ras)) & !is.na(values(max(x.crop))))
              S <- values(x.crop)[pres, ]
              nS <- nrow(S)
              Rg <- x@cov
              p <- values(s.dat.ras)[pres]
              p.sum <- sum(p)
              mar <- apply(S, 2, function(x) sum(x * p)) / p.sum
              if (!quiet) cat("Calculating species covariance matrix...\n")
              Sm <- sweep(S, 2, mar)
              DpSm <- apply(Sm, 2, function(x) x * p)
              Rs <- crossprod(Sm, DpSm) * 1 / (p.sum - 1)
            } else {
              x.mask <- mask(x.crop, s.dat.ras)
              Rg <- x@cov
              p.sum <- cellStats(s.dat.ras, sum)
              DpS <- x.mask * s.dat.ras
              mar <- cellStats(DpS, sum) / p.sum
              if (!quiet) cat("Calculating species covariance matrix...\n")
              Sm <- parScale(x.mask, center = mar, scale = F, parallel = parallel, n = n, quiet = T)
              Rs <- parCov(x = Sm, w = s.dat.ras, parallel = parallel, n = n, quiet = quiet)
            }

            cZ <- nlayers(ras)
            m <- tryCatch(sqrt(as.numeric(t(mar) %*% solve(Rg, tol = 1e-10) %*% mar)),
                          error = function(e){
                            message("Warning: global covariance matrix not invertible. Overall marginality and sensitivity will not be computed.")
                            return(as.numeric(NA))})
            if (max(Im(eigen(Rs)$values)) > 1e-05) stop("complex eigenvalues. Try removing correlated variables.")
            eigRs <- lapply(eigen(Rs), Re)
            Rs12 <- eigRs$vectors %*% diag(eigRs$values^(-0.5)) %*% t(eigRs$vectors)
            W <- Rs12 %*% Rg %*% Rs12
            z <- Rs12 %*% mar
            y <- z/sqrt(sum(z^2))
            H <- (diag(cZ) - y %*% t(y)) %*% W %*% (diag(cZ) - y %*% t(y))
            s <- Re(eigen(H)$values)[-cZ]
            s.p <- (t(mar) %*% Rg %*% mar) / (t(mar) %*% Rs %*% mar)
            s <- c(s.p, s)
            s.p <- abs(s)/sum(abs(s))
            v <- Re(eigen(H)$vectors)
            co <- matrix(nrow = cZ, ncol = cZ)
            u <- as.matrix((Rs12 %*% v)[, 1:(cZ-1)])
            norw <- sqrt(diag(t(u) %*% u))
            co[, -1] <- sweep(u, 2, norw, "/")
            if(is.na(m)) {
              co[, 1] <- mar / norm(mar, "2")
              sf <- as.numeric(abs(co) %*% s.p)
              sens <- as.numeric(NA)
            } else {
              co[, 1] <- mar / m
              sf <- as.numeric(abs(co) %*% s.p)
              sens <- sqrt(as.numeric(sf %*% solve(Rg) %*% sf))
            }
            nm <- c("Marg", paste0("Spec", (1:(cZ-1))))
            if (canProcessInMemory(x.crop) & !parallel){
              s.ras <- brick(x.crop)
              if (!quiet) cat("Creating factor rasters...")
              values(s.ras)[pres, ] <- S %*% co
              names(s.ras) <- nm
            } else {
              if (!quiet) cat("Creating factor rasters...")
              s.ras <- .calc(x.mask, function(x) {x %*% co}, forceapply = T, filename = filename, names = nm, ...)
            }
            colnames(co) <- names(s.p) <- nm
            rownames(co) <- names(sf) <- names(ras)

            cnfa <- methods::new("cnfa", call = call, mf = mar, marginality = m, sf = sf,
                                 sensitivity = sens, p.spec = s.p, co = co, cov = Rs, ras = s.ras, weights = s.dat.ras)
            return(cnfa)
          }
)

#' @rdname cnfa
setMethod("cnfa",
          signature(x = "Raster", s.dat = "Spatial"),
          function(x, s.dat, field, fun = "last", scale = TRUE, filename = "", quiet = TRUE, parallel = FALSE, n = 1, ...){

            call <- sys.calls()[[1]]

            if (! inherits(x, 'Raster')) stop('"x" should be a "Raster*" object')
            if (! inherits(s.dat, c('SpatialPolygons', 'SpatialPoints'))) stop('"s.dat" should be a "SpatialPolygons*" or "SpatialPoints*" object')
            if(!identicalCRS(x, s.dat)) stop("projections do not match")
            if(is.null(intersect(extent(x), extent(s.dat)))) stop("climate and species data do not overlap")
            if(raster::union(extent(x), extent(s.dat)) != extent(x)) stop("extent of species data not contained within extent of climate data")

            if(parallel){
              if(n < 1 | !is.numeric(n)) {
                n <- parallel::detectCores() - 1
                message('incorrect number of cores specified, using ', n)
              } else if(n > parallel::detectCores()) {
                n <- parallel::detectCores() - 1
                message('too many cores specified, using ', n)
              }
            }

            s.dat.ras <- rasterize(s.dat, raster(x), field = field, fun = fun)

            if(scale) {
              if (!quiet) cat("Scaling raster data...\n")
              x <- parScale(x, center = T, scale = T, parallel = parallel, n = n, quiet = quiet)
            }

            filename <- trim(filename)
            if (!canProcessInMemory(x) && filename == '') {
              filename <- rasterTmpFile()
            }

            if(canProcessInMemory(x) & !parallel){
              gpres <- which(!is.na(values(max(x))))
              pres <- which(!is.na(values(s.dat.ras)) & !is.na(values(max(x))))
              Z <- values(x)[gpres, ]
              S <- values(x)[pres, ]
              nZ <- nrow(Z)
              nS <- nrow(S)
              if (!quiet) cat("Calculating global covariance matrix...\n")
              Rg <- crossprod(Z, Z)/(nZ - 1)
              p <- values(s.dat.ras)[pres]
              p.sum <- sum(p)
              mar <- apply(S, 2, function(x) sum(x * p)) / p.sum
              if (!quiet) cat("Calculating species covariance matrix...\n")
              Sm <- sweep(S, 2, mar)
              DpSm <- apply(Sm, 2, function(x) x * p)
              Rs <- crossprod(Sm, DpSm) * 1 / (p.sum - 1)
            } else {
              center <- cellStats(x, mean)
              x.mask <- mask(x, s.dat.ras)
              p.sum <- cellStats(s.dat.ras, sum)
              DpS <- x.mask * s.dat.ras
              mar <- cellStats(DpS, sum) / p.sum
              if (!quiet) cat("Calculating global covariance matrix...\n")
              Rg <- parCov(x, parallel = parallel, n = n, quiet = quiet)
              if (!quiet) cat("Calculating species covariance matrix...\n")
              Sm <- parScale(x.mask, center = mar, scale = F, parallel = parallel, n = n, quiet = T)
              Rs <- parCov(x = Sm, w = s.dat.ras, parallel = parallel, n = n, quiet = quiet)
            }

            cZ <- nlayers(x)
            m <- tryCatch(sqrt(as.numeric(t(mar) %*% solve(Rg, tol = 1e-10) %*% mar)),
                          error = function(e){
                            message("Warning: global covariance matrix not invertible. Overall marginality and sensitivity will not be computed.")
                            return(as.numeric(NA))})
            if(max(Im(eigen(Rs)$values)) > 1e-05) stop("complex eigenvalues. Try removing correlated variables.")
            eigRs <- lapply(eigen(Rs), Re)
            Rs12 <- eigRs$vectors %*% diag(eigRs$values^(-0.5)) %*% t(eigRs$vectors)
            W <- Rs12 %*% Rg %*% Rs12
            z <- Rs12 %*% mar
            y <- z/sqrt(sum(z^2))
            H <- (diag(cZ) - y %*% t(y)) %*% W %*% (diag(cZ) - y %*% t(y))
            s <- Re(eigen(H)$values)[-cZ]
            s.p <- (t(mar) %*% Rg %*% mar) / (t(mar) %*% Rs %*% mar)
            s <- c(s.p, s)
            s.p <- abs(s)/sum(abs(s))
            v <- Re(eigen(H)$vectors)
            co <- matrix(nrow = cZ, ncol = cZ)
            u <- as.matrix((Rs12 %*% v)[, 1:(cZ-1)])
            norw <- sqrt(diag(t(u) %*% u))
            co[, -1] <- sweep(u, 2, norw, "/")
            if(is.na(m)) {
              co[, 1] <- mar / norm(mar, "2")
              sf <- as.numeric(abs(co) %*% s.p)
              sens <- as.numeric(NA)
            } else {
              co[, 1] <- mar / m
              sf <- as.numeric(abs(co) %*% s.p)
              sens <- sqrt(as.numeric(sf %*% solve(Rg) %*% sf))
            }
            nm <- c("Marg", paste0("Spec", (1:(cZ-1))))
            if (canProcessInMemory(x) & !parallel) {
              s.ras <- brick(x)
              if (!quiet) cat("Creating factor rasters...")
              values(s.ras)[pres, ] <- S %*% co
              names(s.ras) <- nm
            } else {
              if (!quiet) cat("Creating factor rasters...")
              s.ras <- .calc(x.mask, function(x) {x %*% co}, forceapply = T, filename = filename, names = nm, ...)
            }
            colnames(co) <- names(s.p) <- nm
            rownames(co) <- names(sf) <- names(x)

            cnfa <- methods::new("cnfa", call = call, mf = mar, marginality = m, sf = sf,
                                 sensitivity = sens, p.spec = s.p, co = co, cov = Rs, ras = s.ras, weights = s.dat.ras)
            return(cnfa)
          }
)

#' @rdname cnfa
setMethod("cnfa",
          signature(x = "Raster", s.dat = "sf"),
          function(x, s.dat, field, fun = "last", scale = TRUE, filename = "", quiet = TRUE, parallel = FALSE, n = 1, ...){
            if (!requireNamespace("sf")) {
              warning('cannot do this because sf is not available')
            }

            call <- sys.calls()[[1]]

            if (! inherits(x, 'Raster')) stop('"x" should be a "Raster*" object')
            if(!identicalCRS(x, s.dat)) stop("projections do not match")
            if(is.null(intersect(extent(x), extent(s.dat)))) stop("climate and species data do not overlap")
            if(raster::union(extent(x), extent(s.dat)) != extent(x)) stop("extent of species data not contained within extent of climate data")

            s.dat <- as(s.dat, "Spatial")
            if (! inherits(s.dat, c('SpatialPolygons', 'SpatialPoints'))) stop('geometry of "s.dat" should be of class "sfc_POLYGON", "sfc_MULTIPOLYGON", "sfc_POINT", or "sfc_MULTIPOINT"')
            if(parallel){
              if(n < 1 | !is.numeric(n)) {
                n <- parallel::detectCores() - 1
                message('incorrect number of cores specified, using ', n)
              } else if(n > parallel::detectCores()) {
                n <- parallel::detectCores() - 1
                message('too many cores specified, using ', n)
              }
            }

            s.dat.ras <- rasterize(s.dat, raster(x), field = field, fun = fun)

            if(scale) {
              if (!quiet) cat("Scaling raster data...\n")
              x <- parScale(x, center = T, scale = T, parallel = parallel, n = n, quiet = quiet)
            }

            filename <- trim(filename)
            if (!canProcessInMemory(x) && filename == '') {
              filename <- rasterTmpFile()
            }

            if(canProcessInMemory(x) & !parallel){
              gpres <- which(!is.na(values(max(x))))
              pres <- which(!is.na(values(s.dat.ras)) & !is.na(values(max(x))))
              Z <- values(x)[gpres, ]
              S <- values(x)[pres, ]
              nZ <- nrow(Z)
              nS <- nrow(S)
              if (!quiet) cat("Calculating global covariance matrix...\n")
              Rg <- crossprod(Z, Z)/(nZ - 1)
              p <- values(s.dat.ras)[pres]
              p.sum <- sum(p)
              mar <- apply(S, 2, function(x) sum(x * p)) / p.sum
              if (!quiet) cat("Calculating species covariance matrix...\n")
              Sm <- sweep(S, 2, mar)
              DpSm <- apply(Sm, 2, function(x) x * p)
              Rs <- crossprod(Sm, DpSm) / (p.sum - 1)
            } else {
              center <- cellStats(x, mean)
              x.mask <- mask(x, s.dat.ras)
              p.sum <- cellStats(s.dat.ras, sum)
              DpS <- x.mask * s.dat.ras
              mar <- cellStats(DpS, sum) / p.sum
              if (!quiet) cat("Calculating global covariance matrix...\n")
              Rg <- parCov(x, parallel = parallel, n = n, quiet = quiet)
              if (!quiet) cat("Calculating species covariance matrix...\n")
              Sm <- parScale(x.mask, center = mar, scale = F, parallel = parallel, n = n, quiet = T)
              Rs <- parCov(x = Sm, w = s.dat.ras, parallel = parallel, n = n, quiet = quiet)
            }

            cZ <- nlayers(x)
            m <- tryCatch(sqrt(as.numeric(t(mar) %*% solve(Rg, tol = 1e-10) %*% mar)),
                          error = function(e){
                            message("Warning: global covariance matrix not invertible. Overall marginality and sensitivity will not be computed.")
                            return(as.numeric(NA))})
            if(max(Im(eigen(Rs)$values)) > 1e-05) stop("complex eigenvalues. Try removing correlated variables.")
            eigRs <- lapply(eigen(Rs), Re)
            Rs12 <- eigRs$vectors %*% diag(eigRs$values^(-0.5)) %*% t(eigRs$vectors)
            W <- Rs12 %*% Rg %*% Rs12
            z <- Rs12 %*% mar
            y <- z/sqrt(sum(z^2))
            H <- (diag(cZ) - y %*% t(y)) %*% W %*% (diag(cZ) - y %*% t(y))
            s <- Re(eigen(H)$values)[-cZ]
            s.p <- (t(mar) %*% Rg %*% mar) / (t(mar) %*% Rs %*% mar)
            s <- c(s.p, s)
            s.p <- abs(s)/sum(abs(s))
            v <- Re(eigen(H)$vectors)
            co <- matrix(nrow = cZ, ncol = cZ)
            u <- as.matrix((Rs12 %*% v)[, 1:(cZ-1)])
            norw <- sqrt(diag(t(u) %*% u))
            co[, -1] <- sweep(u, 2, norw, "/")
            sf <- as.numeric( abs(co) %*% s.p )
            if(is.na(m)) {
              co[, 1] <- mar / norm(mar, "2")
              sens <- as.numeric(NA)
            } else {
              co[, 1] <- mar / m
              sens <- sqrt(as.numeric(sf %*% solve(Rg) %*% sf))
            }
            nm <- c("Marg", paste0("Spec", (1:(cZ-1))))
            if (canProcessInMemory(x) & !parallel){
              s.ras <- brick(x)
              if (!quiet) cat("Creating factor rasters...")
              values(s.ras)[pres, ] <- S %*% co
              names(s.ras) <- nm
            } else {
              if (!quiet) cat("\nCreating factor rasters...")
              s.ras <- .calc(x.mask, function(x) {x %*% co}, forceapply = T, filename = filename, names = nm, ...)
            }
            colnames(co) <- names(s.p) <- nm
            rownames(co) <- names(sf) <- names(x)

            cnfa <- methods::new("cnfa", call = call, mf = mar, marginality = m, sf = sf,
                                 sensitivity = sens, p.spec = s.p, co = co, cov = Rs, ras = s.ras, weights = s.dat.ras)
            return(cnfa)
          }
)


# maybe for a later version
#
# @rdname cnfa
# setMethod("cnfa",
#           signature(x = "Raster", s.dat = "matrix"),
#           function(x, s.dat, prj, field, fun = "last", scale = TRUE, filename = "", parallel = F, ...){
#
#             call <- match.call()
#
#             if(!inherits(x, 'Raster')) stop('"x" should be a "Raster*" object')
#             if(missing(prj)) {
#               prj <- crs(x)
#               warning("prj argument not supplied: inheriting projection from Raster* object.")
#             }
#             s.dat <- SpatialPointsDataFrame(coords = s.dat, data = 1, proj4string = crs(prj))
#             if(!identicalCRS(x, s.dat)) stop("projections do not match")
#
#             if(is.null(intersect(extent(x), extent(s.dat)))) stop("climate and species data do not overlap")
#
#             s.dat.ras <- rasterize(s.dat, raster(x), field = field, fun = fun, ...)
#
#             if(scale == TRUE) x <- scale(x)
#
#             filename <- trim(filename)
#             if (!canProcessInMemory(x) && filename == '') {
#               filename <- rasterTmpFile()
#             }
#
#             if(canProcessInMemory(x)){
#               gpres <- which(!is.na(values(max(x))))
#               pres <- which(!is.na(values(s.dat.ras)) & !is.na(values(max(x))))
#               Z <- values(x)[gpres, ]
#               S <- values(x)[pres, ]
#               nZ <- nrow(Z)
#               nS <- nrow(S)
#               Rg <- crossprod(Z, Z)/nZ
#               p <- values(s.dat.ras)[pres]
#               p <- p/(sum(p))
#               mar <- apply(S, 2, function(x) sum(x * p))
#               Sm <- sweep(S, 2, mar)
#               DpSm <- apply(Sm, 2, function(x) x * p)
#               Rs <- crossprod(Sm, DpSm)
#             } else {
#               center <- cellStats(x, mean)
#               x.mask <- mask(x, s.dat.ras)
#               #pres <- which(!is.na(values(s.dat.ras)) & !is.na(values(max(x.mask))))
#               p.sum <- cellStats(s.dat.ras, sum)
#               p <- s.dat.ras / p.sum
#               DpS <- x.mask * p
#               mar <- cellStats(DpS, sum)
#               cat("\nCalculating study area covariance matrix...\n")
#               Rg <- parCov(x, sample = F, ...)
#               cat("\nCalculating species covariance matrix...\n")
#               Sm <- calc(x.mask, fun = function(x) x - mar, forceapply = T)
#               Rs <- parCov(x = Sm, w = s.dat.ras, parallel = parallel, ...)
#               #Rs <- layerStats(Sm, 'weighted.cov', w = p, na.rm = T, asSample = F)[[1]]
#               #Rs <- parCov(x.mask, center = T)
#               #mar <- cellStats(x.mask, sum)/length(pres)
#             }
#
#             cZ <- nlayers(x)
#             m <- sqrt(c(t(mar) %*% mar))
#             if(max(Im(eigen(Rs)$values)) > 1e-05) stop("complex eigenvalues. Try removing correlated variables.")
#             eigRs <- lapply(eigen(Rs), Re)
#             keep <- (eigRs$values > 1e-09)
#             Rs12 <- eigRs$vectors[, keep] %*% diag(eigRs$values[keep]^(-0.5)) %*% t(eigRs$vectors[, keep])
#             W <- Rs12 %*% Rg %*% Rs12
#             z <- Rs12 %*% mar
#             y <- z/sqrt(sum(z^2))
#             H <- (diag(cZ) - y %*% t(y)) %*% W %*% (diag(cZ) - y %*% t(y))
#             s <- Re(eigen(H)$values)[-cZ]
#             s.p <- abs(sum(diag(W)) - sum(diag(H)))
#             s <- c(s.p, s)
#             s.p <- abs(s)/sum(abs(s))
#             v <- Re(eigen(H)$vectors)
#             co <- matrix(nrow = cZ, ncol = cZ)
#             co[, 1] <- mar/sqrt(t(mar) %*% mar)
#             u <- as.matrix((Rs12 %*% v)[, 1:(cZ-1)])
#             norw <- sqrt(diag(t(u) %*% u))
#             co[, -1] <- sweep(u, 2, norw, "/")
#             sf <- abs(co) %*% s.p
#             sf <- as.numeric(sf)
#             names(sf) <- names(x)
#             sens <- norm(sf, "2")
#             nm <- c("Marg", paste0("Spec", (1:(cZ-1))))
#             if(canProcessInMemory(x)){
#               s.ras <- brick(x)
#               values(s.ras)[pres, ] <- S %*% co
#               #setValues(s.ras, S %*% co, index = pres)
#               names(s.ras) <- nm
#             } else{
#               cat("\nCreating factor rasters...")
#               s.ras <- .calc(x.mask, function(x) {x %*% co}, forceapply = T, filename = filename, names = nm, ...)
#             }
#             colnames(co) <- names(s.p) <- nm
#             rownames(co) <- names(x)
#
#             cnfa <- methods::new("cnfa", call = call, mf = mar, marginality = m, sf = sf,
#                                  sensitivity = sens, p.spec = s.p, co = co, cov = Rs, ras = s.ras, weights = s.dat.ras)
#             return(cnfa)
#           }
# )
