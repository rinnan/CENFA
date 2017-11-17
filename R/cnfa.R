#' Climate Niche Factor Analysis
#'
#' Performs climate niche factor analysis using climate raster data and species presence data.
#'
#' @aliases print.cnfa, show.cnfa
#' @param x Raster* or GLcnfa object, typically a brick or stack with p climate raster layers
#' @param s.dat Spatial* object detailing species presence or abundance
#' @param field field of \code{s.dat} that specifies presence or abundance. This is equivalent to the \code{field} argument in \code{raster::rasterize}.
#' @param fun function or character. Determines what values to assign to cells with multiple spatial features, similar to the \code{fun} argument in \code{raster::rasterize}.
#' @param nf integer or character. Specifies the number of specialization axes to keep after transformation. If \code{nf = "BS"} then the broken-stick method will be used to determine the number of significant factors.
#' @param scale logical. If \code{TRUE} then the values of \code{x} will
#' get centered and scaled. Depending on the resolution of the climate data and the
#' extent of the study area, this can be quite time consuming. If running this
#' function for multiple species, it is recommended that the data be scaled beforehand.
#' @param sp.prj character. Spatial projection of species data.
#' @param mar.type character. Choices are "BC" (Basille) or "H" (Hirzel). See details.
#' @details Hirzel et al. (2002) defined the overall marginality M as a standardized distance between the centroid of the species' niche and the global niche, given by \deqn{M = \frac{\sqrt{\sum\limits_{j=1}^P m_j^2}}{1.96}.} Basille and Calenge (2008), however, defined M as \deqn{M = \mathbf{m}^T\mathbf{m}.} The default \code{mar.type} reflects Basille and Calenge's definition so that results will by default agree with those calculated using Basille and Calenge's \code{adehabitatHS} package.
#' @return Returns an S4 object of class \code{cnfa} with the following components:
#' @return call Original function call.
#' @return nf Number of specialization factors kept.
#' @return mf marginality factor. Vector that describes the location of the species Hutchinsonian niche relative to the global niche.
#' @return marginality Standardized magnitude of the marginality factor.
#' @return sf specialization factors. Vector of eigenvalues.
#' @return specialization The square root of the sum of \code{sf} divided by the length of \code{sf}.
#' @return s.prop Vector representing the amount of specialization found in each factor.
#' @return co A data frame of variable loadings with p rows and nf + 1 columns.
#' @return ras Raster* object of transformed climate values, with number of layers equal to nf + 1.
#' @return present numeric. Number of raster cells in which species is present.
#' @export
#'
#'

setGeneric("cnfa", function(x, s.dat, nf = "BS", mar.type = "BC", ...) {
  standardGeneric("cnfa")
})

#' @rdname cnfa
setMethod("cnfa",
          signature(x = "GLcnfa", s.dat = "SpatialPolygonsDataFrame"),
          function(x, s.dat, field, nf = "BS", mar.type = "BC", filename = "", ...){

            call <- match.call()

            if(!identicalCRS(raster(x), s.dat)) stop("climate and species projections do not match")
            ext <- extent(raster(x))
            if(is.null(raster::intersect(ext, extent(s.dat)))) stop("climate and species data do not overlap")

            if(raster::union(ext, extent(s.dat)) != ext) stop("extent of species data not contained within extent of climate data")
            x.crop <- crop(raster(x), extent(s.dat))
            s.dat.ras <- rasterize(s.dat, x.crop, field = field)

            filename <- trim(filename)
            if (!canProcessInMemory(x.crop) && filename == '') {
              filename <- rasterTmpFile()
            }

            if(canProcessInMemory(x.crop)){
              pres <- which(!is.na(values(s.dat.ras)) & !is.na(values(max(x.crop))))
              S <- values(x.crop)[pres,]
              nS <- nrow(S)
              Rg <- x@cov
              #Rs <- crossprod(S, S)/(nS - 1)
              Rs <- cov(S)
              mar <- colSums(S)/nS
            } else {
              x.mask <- mask(x.crop, s.dat.ras)
              pres <- which(!is.na(values(s.dat.ras)) & !is.na(values(max(x.mask))))
              Rg <- x@cov
              cat("\nCalculating covariance matrix...\n")
              Rs <- covmat(x.mask, ...)
              mar <- cellStats(x.mask, sum)/length(pres)
            }

            cZ <- nlayers(raster(x))
            if(mar.type == "BC") m <- c(t(mar) %*% mar)
            if(mar.type == "H")  m <- norm(mar, "2")/1.96
            if(max(Im(eigen(Rs)$values)) > 1e-05) stop("complex eigenvalues")
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
            s.p <- abs(s.p)/sum(abs(s.p))
            spec <- sqrt(sum(s))/cZ
            v <- Re(eigen(H)$vectors)
            if (nf == "BS") nf <- brStick(s[-1])
            if (nf <= 0 | nf > (cZ - 1)) nf <- 1
            co <- matrix(nrow = cZ, ncol = nf + 1)
            co[, 1] <- mar
            u <- as.matrix((Rs12 %*% v)[, 1:nf])
            norw <- sqrt(diag(t(u) %*% u))
            co[, -1] <- sweep(u, 2, norw, "/")
            if(canProcessInMemory(x.crop)){
              ras <- brick(x.crop, nl = nf + 1)
              values(ras)[pres, ] <- S %*% co
            } else{
              cat("\nCreating raster of transformed variables...")
              ras <- calc(x.mask, function(x) {x %*% co}, forceapply = T, filename = filename, ...)
            }
            co <- as.data.frame(co)
            names(co) <- c("Marg", paste0("Spec", (1:nf)))
            row.names(co) <- names(raster(x))

            # pres <- which(!is.na(values(s.dat.ras)) & !is.na(values(x.crop[[1]])))
            # S <- values(x.crop)[pres, ]
            # rZ <- x@ncells
            # cZ <- nlayers(x@global_ras)
            # rS <- nrow(S)
            # mar <- colSums(S)/rS
            # if(mar.type == "BC") m <- c(t(mar) %*% mar)
            # if(mar.type == "H")  m <- norm(mar, "2")/1.96
            # Rg <- x@cov
            # Rs <- crossprod(S,S/rS)
            # eigRs <- eigen(Rs)
            # keep <- (eigRs$values > 1e-09)
            # Rs12 <- eigRs$vectors[, keep] %*% diag(eigRs$values[keep]^(-0.5)) %*% t(eigRs$vectors[, keep])
            # W <- Rs12 %*% Rg %*% Rs12
            # z <- Rs12 %*% mar
            # y <- z/sqrt(sum(z^2))
            # H <- (diag(cZ) - y %*% t(y)) %*% W %*% (diag(cZ) - y %*% t(y))
            # s <- eigen(H)$values[-cZ]
            # s.p <- abs(s)/sum(abs(s))
            # s.p[1] <- sum(diag(W)) - sum(diag(H))
            # spec <- sqrt(sum(s))/cZ
            # v <- eigen(H)$vectors
            # if (nf<=0 | nf>(cZ-1)){nf <- 1}
            # co <- matrix(nrow = cZ, ncol = nf + 1)
            # u <- (Rs12 %*% v)[, 1:nf]
            # norw <- sqrt(diag(t(as.matrix(u)) %*% as.matrix(u)))
            # co[, 2:(nf + 1)] <- sweep(as.matrix(u), 2, norw, "/")
            # co[, 1] <- mar
            # ras <- brick(x.crop, nl = nf + 1)
            # values(ras)[pres,] <- S %*% co
            # #names(ras) <- c("Marg", paste0("Spec", (1:nf)))
            # co <- as.data.frame(co)
            # names(co) <- c("Marg", paste0("Spec", (1:nf)))
            # row.names(co) <- names(x@global_ras)

            cnfa <- methods::new("cnfa", call = call, nf = nf, mf = mar, marginality = m, sf = s, specialization = spec, s.prop = s.p, co = co, ras = ras, s.cov = Rs, present = length(pres))
            return(cnfa)
          }
)

#' @rdname cnfa
setMethod("cnfa",
          signature(x = "RasterBrick", s.dat = "SpatialPolygonsDataFrame"),
          function(x, s.dat, field, scale = FALSE, mar.type = "BC", ...){
            call <- match.call()

            if(!identicalCRS(x, s.dat)) stop("projections do not match")

            if(is.null(raster::intersect(extent(x), extent(s.dat)))) stop("climate and species data do not overlap")

            if(raster::union(extent(x), extent(s.dat)) != extent(x)) stop("extent of species data not contained within extent of climate data")

            if(scale == TRUE) x <- raster::scale(x)

            s.dat.ras <- rasterize(s.dat, raster(x), field = field)

            if(canProcessInMemory(x)){
              gpres <- which(!is.na(values(max(x))))
              pres <- which(!is.na(values(s.dat.ras)))
              Z <- values(x)[gpres,]
              S <- values(x)[pres,]
              nZ <- nrow(Z)
              nS <- nrow(S)
              center <- colMeans(Z)
              Z <- sweep(Z, 2, center)
              S <- sweep(S, 2, center)
              mar <- colSums(S)/nS
              Rg <- crossprod(Z, Z)/nZ
              Rs <- crossprod(S, S)/(nS - 1)
            } else {
              center <- cellStats(x, mean)
              x.mask <- mask(x, s.dat.ras)
              pres <- which(!is.na(values(s.dat.ras)) & !is.na(values(max(x.mask))))
              Z <- calc(x, fun = function(p) {p - center})
              S <- calc(x.mask, fun = function(p) {p - center})
              mar <- cellStats(x.mask, sum)/length(pres)
              Rg <- covmat(Z, sample = F, ...)
              Rs <- covmat(S, ...)
            }


            if(mar.type == "BC") m <- c(t(mar) %*% mar)
            if(mar.type == "H")  m <- norm(mar, "2")/1.96
            eigRs <- eigen(Rs)
            keep <- (eigRs$values > 1e-09)
            Rs12 <- eigRs$vectors[, keep] %*% diag(eigRs$values[keep]^(-0.5)) %*% t(eigRs$vectors[, keep])
            W <- Rs12 %*% Rg %*% Rs12
            z <- Rs12 %*% mar
            y <- z/sqrt(sum(z^2))
            H <- (diag(ncol(Z)) - y %*% t(y)) %*% W %*% (diag(ncol(Z)) - y %*% t(y))
            s <- eigen(H)$values[-ncol(Z)]
            spec <- sqrt(sum(s))/ncol(Z)
            s.p <- abs(sum(diag(W)) - sum(diag(H)))
            s.p <- c(s.p, s)
            s.p <- abs(s.p)/sum(abs(s.p))
            v <- eigen(H)$vectors
            if (nf == "BS") nf <- brStick(s[-1])
            if (nf <= 0 | nf > (ncol(Z) - 1)) nf <- 1
            co <- matrix(nrow = ncol(Z), ncol = nf + 1)
            co[, 1] <- mar
            u <- as.matrix((Rs12 %*% v)[, 1:nf])
            norw <- sqrt(diag(t(u) %*% u))
            co[, -1] <- sweep(u, 2, norw, "/")
            if(canProcessInMemory(x)){
              ras <- brick(x, nl = nf + 1)
              values(ras)[pres, ] <- S %*% co
            } else{
              ras <- calc(x.mask, function(x) {x %*% co}, forceapply = T)
            }
            co <- as.data.frame(co)
            names(co) <- c("Marg", paste0("Spec", (1:nf)))
            row.names(co) <- names(x)

            cnfa <- methods::new("cnfa", call = call, nf = nf, mf = mar, marginality = m, sf = s, specialization = spec, s.prop = s.p, co = co, ras = ras, s.cov = Rs, present = pres)
            return(cnfa)
          }
)

#' @rdname cnfa
setMethod("cnfa",
          signature(x = "Raster", s.dat = "SpatialPoints"),
          function(x, s.dat, field, fun = "count", scale = FALSE){
            call <- match.call()


            if(!identicalCRS(x, s.dat)) stop("projections do not match")

            if(is.null(raster::intersect(extent(x), extent(s.dat)))) stop("climate and species data do not overlap")

            if(scale == TRUE) x <- raster::scale(x)

            gpres <- which(!is.na(values(x[[1]])))
            dat <- values(x)[gpres, ]

            s.dat.ras <- rasterize(s.dat, raster(x), field = field, fun = fun)
            pres <- which(!is.na(values(s.dat.ras)))
            #prb<-values(speciesdat.ras)[pres]
            pres.dat <- values(x)[pres, ]

            #prb<-na.omit(values(speciesdat.ras))
            #prb <- pr
            #pr <- prb/sum(prb)
            #row.w<-rep(1,nrow(dat))/nrow(dat)
            #col.w <- rep(1,ncol(dat))
            #n <- nrow(dat)
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
            class(cnfa) <- "cnfa"
            return(invisible(cnfa))
          }
)

#' @rdname cnfa
setMethod("cnfa",
          signature(x = "RasterBrick", s.dat = "matrix"),
          function(x, s.dat, field, nf = 1, scale = FALSE, sp.prj){
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
