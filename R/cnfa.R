#' Climate Niche Factor Analysis
#'
#' Performs climate niche factor analysis using climate raster data and species presence data.
#'
#' @aliases print.cnfa, show.cnfa
#' @param x Raster* or GLcenfa object, typically a brick or stack with p climate raster layers
#' @param s.dat Spatial* object detailing species presence or abundance
#' @param field field of \code{s.dat} that specifies presence or abundance. This is equivalent to the \code{field} argument in \code{raster::rasterize}.
#' @param fun function or character. Determines what values to assign to cells with multiple spatial features, similar to the \code{fun} argument in \code{raster::rasterize}.
#' @param scale logical. If \code{TRUE} then the values of \code{x} will
#' get centered and scaled. Depending on the resolution of the climate data and the
#' extent of the study area, this can be quite time consuming. If running this
#' function for multiple species, it is recommended that the data be scaled beforehand.
#' @details Hirzel et al. (2002) defined the overall marginality M as a standardized distance between the centroid of the species' niche and the global niche, given by \deqn{M = \frac{\sqrt{\sum\limits_{j=1}^P m_j^2}}{1.96}.} Basille and Calenge (2008), however, defined M as \deqn{M = \mathbf{m}^T\mathbf{m}.} The default \code{mar.type} reflects Basille and Calenge's definition so that results will by default agree with those calculated using Basille and Calenge's \code{adehabitatHS} package.
#' @return Returns an S4 object of class \code{cnfa} with the following components:
#' @return call Original function call.
#' @return mf marginality factor. Vector that describes the location of the species Hutchinsonian niche relative to the global niche.
#' @return marginality Magnitude of the marginality factor.
#' @return sf sensitivity factor.
#' @return sensitivity Magnitude of the sensitivity factor.
#' @return s.prop Vector representing the amount of specialization found in each ENFA factor.
#' @return co A matrix describing the amount of marginality and specialization on each factor. (The marginality column is normalized.)
#' @return present numeric. Number of raster cells in which species is present.
#' @return ras Raster* object of transformed climate values.
#' @export
#'
#' @importFrom raster brick canProcessInMemory cellStats crop extent mask nlayers raster rasterize rasterTmpFile trim values
#' @importFrom sp CRS identicalCRS SpatialPointsDataFrame SpatialPolygonsDataFrame
#' @importFrom stats cov

setGeneric("cnfa", function(x, s.dat, field, filename = "", ...) {
  standardGeneric("cnfa")
})

#' @rdname cnfa
setMethod("cnfa",
          signature(x = "GLcenfa", s.dat = "SpatialPolygonsDataFrame"),
          function(x, s.dat, field, filename = "", ...){

            call <- match.call()

            if(!identicalCRS(raster(x), s.dat)) stop("climate and species projections do not match")
            ras <- raster(x)
            ext <- extent(ras)
            if(is.null(raster::intersect(ext, extent(s.dat)))) stop("climate and species data do not overlap")

            if(raster::union(ext, extent(s.dat)) != ext) stop("extent of species data not contained within extent of climate data")
            x.crop <- crop(ras, extent(s.dat))
            s.dat.ras <- rasterize(s.dat, x.crop, field = field, ...)

            filename <- trim(filename)
            if (!canProcessInMemory(x.crop) && filename == '') {
              filename <- rasterTmpFile()
            }

            if(canProcessInMemory(x.crop)){
              pres <- which(!is.na(values(s.dat.ras)) & !is.na(values(max(x.crop))))
              S <- values(x.crop)[pres,]
              nS <- nrow(S)
              Rg <- x@cov
              Rs <- cov(S)
              mar <- colSums(S)/nS
            } else {
              x.mask <- mask(x.crop, s.dat.ras)
              pres <- which(!is.na(values(s.dat.ras)) & !is.na(values(max(x.mask))))
              Rg <- x@cov
              cat("\nCalculating species covariance matrix...\n")
              Rs <- covmat(x.mask, center = T)
              mar <- cellStats(x.mask, sum)/length(pres)
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
            s <- Re(eigen(H)$values)[-cZ]
            s.p <- abs(sum(diag(W)) - sum(diag(H)))
            s <- c(s.p, s)
            s.p <- abs(s)/sum(abs(s))
            v <- Re(eigen(H)$vectors)
            co <- matrix(nrow = cZ, ncol = cZ)
            co[, 1] <- mar/sqrt(t(mar) %*% mar)
            u <- as.matrix((Rs12 %*% v)[, 1:(cZ-1)])
            norw <- sqrt(diag(t(u) %*% u))
            co[, -1] <- sweep(u, 2, norw, "/")
            sf <- abs(co) %*% s.p
            sf <- as.numeric(sf)
            names(sf) <- names(ras)
            sens <- norm(sf, "2")
            nm <- c("Marg", paste0("Spec", (1:(cZ-1))))
             if(canProcessInMemory(x.crop)){
               s.ras <- brick(x.crop)
               values(s.ras)[pres, ] <- S %*% co
               names(s.ras) <- nm
             } else{
               cat("\nCreating factor rasters...")
               s.ras <- .calc(x.mask, function(x) {x %*% co}, forceapply = T, filename = filename, names = nm, ...)
             }
            colnames(co) <- c("Marg", paste0("Spec", (1:(cZ-1))))
            rownames(co) <- names(ras)

            cnfa <- methods::new("cnfa", call = call, mf = mar, marginality = m, sf = sf, sensitivity = sens, p.spec = s.p, co = co, cov = Rs, present = length(pres), ras = s.ras)
            return(cnfa)
          }
)

#' @rdname cnfa
setMethod("cnfa",
          signature(x = "GLcenfa", s.dat = "SpatialPoints"),
          function(x, s.dat, field, fun = "count", scale = FALSE, filename = "", ...){

            call <- match.call()

            if(!identicalCRS(raster(x), s.dat)) stop("climate and species projections do not match")
            ras <- raster(x)
            ext <- extent(ras)
            if(is.null(raster::intersect(ext, extent(s.dat)))) stop("climate and species data do not overlap")

            if(raster::union(ext, extent(s.dat)) != ext) stop("extent of species data not contained within extent of climate data")
            x.crop <- crop(ras, extent(s.dat))
            s.dat.ras <- rasterize(s.dat, x.crop, field = field, fun = fun)

            filename <- trim(filename)
            if (!canProcessInMemory(x) && filename == '') {
              filename <- rasterTmpFile()
            }

            if(canProcessInMemory(x)){
              pres <- which(!is.na(values(s.dat.ras)) & !is.na(values(max(x))))
              S <- values(x)[pres, ]
              nS <- nrow(S)
              Rg <- x@cov
              Rs <- cov(S)
              mar <- colSums(S)/nS
            } else {
              x.mask <- mask(x.crop, s.dat.ras)
              pres <- which(!is.na(values(s.dat.ras)) & !is.na(values(max(x.mask))))
              Rg <- x@cov
              cat("\nCalculating species covariance matrix...\n")
              Rs <- covmat(x.mask, center = T)
              mar <- cellStats(x.mask, sum)/length(pres)
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
            s <- Re(eigen(H)$values)[-cZ]
            s.p <- abs(sum(diag(W)) - sum(diag(H)))
            s <- c(s.p, s)
            s.p <- abs(s)/sum(abs(s))
            v <- Re(eigen(H)$vectors)
            co <- matrix(nrow = cZ, ncol = cZ)
            co[, 1] <- mar/sqrt(t(mar) %*% mar)
            u <- as.matrix((Rs12 %*% v)[, 1:(cZ-1)])
            norw <- sqrt(diag(t(u) %*% u))
            co[, -1] <- sweep(u, 2, norw, "/")
            sf <- abs(co) %*% s.p
            sf <- as.numeric(sf)
            names(sf) <- names(ras)
            sens <- norm(sf, "2")
            nm <- c("Marg", paste0("Spec", (1:(cZ-1))))
            if(canProcessInMemory(x)){
              s.ras <- brick(x)
              values(s.ras)[pres, ] <- S %*% co
              names(s.ras) <- nm
            } else{
              cat("\nCreating factor rasters...")
              s.ras <- .calc(x.mask, function(x) {x %*% co}, forceapply = T, filename = filename, names = nm, ...)
            }
            colnames(co) <- c("Marg", paste0("Spec", (1:(cZ-1))))
            rownames(co) <- names(ras)

            cnfa <- methods::new("cnfa", call = call, mf = mar, marginality = m, sf = sf, sensitivity = sens, p.spec = s.p, co = co, cov = Rs, present = length(pres), ras = s.ras)
            return(cnfa)
          }
)

#' @rdname cnfa
setMethod("cnfa",
          signature(x = "RasterBrick", s.dat = "SpatialPolygonsDataFrame"),
          function(x, s.dat, field, scale = FALSE, filename = "", ...){

            call <- match.call()

            if(!identicalCRS(x, s.dat)) stop("projections do not match")

            if(is.null(raster::intersect(extent(x), extent(s.dat)))) stop("climate and species data do not overlap")

            if(raster::union(extent(x), extent(s.dat)) != extent(x)) stop("extent of species data not contained within extent of climate data")

            if(scale) x <- raster::scale(x)

            s.dat.ras <- rasterize(s.dat, raster(x), field = field, ...)

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
              Rs <- cov(S)
              mar <- colSums(S)/nS
            } else {
              center <- cellStats(x, mean)
              x.mask <- mask(x, s.dat.ras)
              pres <- which(!is.na(values(s.dat.ras)) & !is.na(values(max(x.mask))))
              cat("\nCalculating study area covariance matrix...\n")
              Rg <- covmat(x, sample = F)
              cat("\nCalculating species covariance matrix...\n")
              Rs <- covmat(x.mask, center = T)
              mar <- cellStats(x.mask, sum)/length(pres)
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
            s <- Re(eigen(H)$values)[-cZ]
            s.p <- abs(sum(diag(W)) - sum(diag(H)))
            s <- c(s.p, s)
            s.p <- abs(s)/sum(abs(s))
            v <- Re(eigen(H)$vectors)
            co <- matrix(nrow = cZ, ncol = cZ)
            co[, 1] <- mar/sqrt(t(mar) %*% mar)
            u <- as.matrix((Rs12 %*% v)[, 1:(cZ-1)])
            norw <- sqrt(diag(t(u) %*% u))
            co[, -1] <- sweep(u, 2, norw, "/")
            sf <- abs(co) %*% s.p
            sf <- as.numeric(sf)
            names(sf) <- names(x)
            sens <- norm(sf, "2")
            nm <- c("Marg", paste0("Spec", (1:(cZ-1))))
            if(canProcessInMemory(x)){
              s.ras <- brick(x)
              values(s.ras)[pres, ] <- S %*% co
              names(s.ras) <- nm
            } else{
              cat("\nCreating factor rasters...")
              s.ras <- .calc(x.mask, function(x) {x %*% co}, forceapply = T, filename = filename, names = nm, ...)
            }
            colnames(co) <- c("Marg", paste0("Spec", (1:(cZ-1))))
            rownames(co) <- names(x)

            cnfa <- methods::new("cnfa", call = call, mf = mar, marginality = m, sf = sf, sensitivity = sens, p.spec = s.p, co = co, cov = Rs, present = length(pres), ras = s.ras)
            return(cnfa)
          }
)

#' @rdname cnfa
setMethod("cnfa",
          signature(x = "RasterBrick", s.dat = "SpatialPoints"),
          function(x, s.dat, field, fun = "count", scale = FALSE, filename = "", ...){

            call <- match.call()

            if(!identicalCRS(x, s.dat)) stop("projections do not match")

            if(is.null(raster::intersect(extent(x), extent(s.dat)))) stop("climate and species data do not overlap")

            if(raster::union(extent(x), extent(s.dat)) != extent(x)) stop("extent of species data not contained within extent of climate data")

            if(scale) x <- raster::scale(x)

            s.dat.ras <- rasterize(s.dat, raster(x), field = field, fun = fun)

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
              Rs <- cov(S)
              mar <- colSums(S)/nS
            } else {
              center <- cellStats(x, mean)
              x.mask <- mask(x, s.dat.ras)
              pres <- which(!is.na(values(s.dat.ras)) & !is.na(values(max(x.mask))))
              cat("\nCalculating study area covariance matrix...\n")
              Rg <- covmat(x, sample = F)
              cat("\nCalculating species covariance matrix...\n")
              Rs <- covmat(x.mask, center = T)
              mar <- cellStats(x.mask, sum)/length(pres)
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
            s <- Re(eigen(H)$values)[-cZ]
            s.p <- abs(sum(diag(W)) - sum(diag(H)))
            s <- c(s.p, s)
            s.p <- abs(s)/sum(abs(s))
            v <- Re(eigen(H)$vectors)
            co <- matrix(nrow = cZ, ncol = cZ)
            co[, 1] <- mar/sqrt(t(mar) %*% mar)
            u <- as.matrix((Rs12 %*% v)[, 1:(cZ-1)])
            norw <- sqrt(diag(t(u) %*% u))
            co[, -1] <- sweep(u, 2, norw, "/")
            sf <- abs(co) %*% s.p
            sf <- as.numeric(sf)
            names(sf) <- names(x)
            sens <- norm(sf, "2")
            nm <- c("Marg", paste0("Spec", (1:(cZ-1))))
            if(canProcessInMemory(x)){
              s.ras <- brick(x)
              values(s.ras)[pres, ] <- S %*% co
              names(s.ras) <- nm
            } else{
              cat("\nCreating factor rasters...")
              s.ras <- .calc(x.mask, function(x) {x %*% co}, forceapply = T, filename = filename, names = nm, ...)
            }
            colnames(co) <- c("Marg", paste0("Spec", (1:(cZ-1))))
            rownames(co) <- names(x)

            cnfa <- methods::new("cnfa", call = call, mf = mar, marginality = m, sf = sf, sensitivity = sens, p.spec = s.p, co = co, cov = Rs, present = length(pres), ras = s.ras)
            return(cnfa)
          }
)

#' @rdname cnfa
setMethod("cnfa",
          signature(x = "RasterBrick", s.dat = "matrix"),
          function(x, s.dat, prj, scale = FALSE, filename = "", ...){

            call <- match.call()

            s.dat <- SpatialPointsDataFrame(coords = s.dat, data = 1, proj4string = CRS(sp.prj))
            if(!identicalCRS(x, s.dat)) stop("projections do not match")

            if(is.null(raster::intersect(extent(x), extent(s.dat)))) stop("climate and species data do not overlap")

            if(scale == TRUE) x <- raster::scale(x)

            gpres <- which(!is.na(values(x[[1]])))
            dat <- values(x)[gpres, ]

            s.dat.ras <- rasterize(s.dat, raster(x), field = field)
            pres <- which(!is.na(values(s.dat.ras)))
            #prb<-values(speciesdat.ras)[pres]
            pres.dat <- values(x)[pres, ]
            #pr <- prb/sum(prb)
            #row.w<-rep(1,nrow(dat))/nrow(dat)
            #col.w <- rep(1,ncol(dat))
            center <- colMeans(dat)
            Z <- sweep(dat, 2, center)
            S <- sweep(pres.dat, 2, center)
            nZ <- nrow(Z)
            nS <- nrow(S)
            mar <- colSums(S)/nS
            if(mar.type == "BC") m <- c(t(mar) %*% mar)
            if(mar.type == "H")  m <- norm(mar, "2")/1.96
            Rg <- crossprod(Z, Z/nZ)
            Rs <- crossprod(S, S/nS)
            eigRs <- eigen(Rs)
            keep <- (eigRs$values > 1e-09)
            Rs12 <- eigRs$vectors[, keep] %*% diag(eigRs$values[keep]^(-0.5)) %*% t(eigRs$vectors[, keep])
            W <- Rs12 %*% Rg %*% Rs12
            z <- Rs12 %*% mar
            y <- z/sqrt(sum(z^2))
            H <- (diag(ncol(Z)) - y %*% t(y)) %*% W %*% (diag(ncol(Z)) - y %*% t(y))
            s <- eigen(H)$values
            s.p <- abs(s)/sum(abs(s))
            s.p[1] <- sum(diag(W)) - sum(diag(H))
            spec <- sqrt(sum(s))/ncol(Z)
            v <- eigen(H)$vectors
            if (nf == "BS") nf <- brStick(v)
            if (nf <= 0 | nf > (ncol(Z) - 1)) nf <- 1
            co <- matrix(nrow = ncol(Z), ncol = nf + 1)
            u <- (Rs12 %*% v)[, 1:nf]
            norw <- sqrt(diag(t(as.matrix(u)) %*% as.matrix(u)))
            co[, 2:(nf + 1)] <- sweep(as.matrix(u), 2, norw, "/")
            #co[, 2:(nf + 1)] <- as.matrix(u)#sweep(as.matrix(u), 2, norw, "/")
            co[, 1] <- mar
            ras <- brick(raster(s.dat.ras), nl = nf + 1)
            ss <- crop(s.dat.ras, extent(s.dat))
            pres.c <- which(!is.na(values(ss)))
            values(ras)[pres.c, ] <- S %*% co
            co <- as.data.frame(co)
            names(co) <- c("Marg", paste0("Spec", (1:nf)))
            row.names(co) <- names(x)
            co <- as.data.frame(co[order(abs(co$Marg), decreasing = T), ])

            cnfa <- methods::new("cnfa", call = call, nf = nf, mf = mar, marginality = m, sf = s, specialization = spec, s.prop = s.p, co = co, ras = ras, s.cov = Rs, present = length(pres))
            return(cnfa)
          }
)
