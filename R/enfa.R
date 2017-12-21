#' Ecological-Niche Factor Analysis
#'
#' Performs ecological-niche factor analysis using environmental raster data and species presence data in shapefile.
#'
#' @aliases print.enfa, show.enfa
#' @param x Raster* object, typically a brick or stack of ecological raster layers
#' @param s.dat SpatialPolygons* object detailing species presence or abundance
#' @param field field of \code{s.dat} that specifies presence or abundance. This is equivalent to the \code{field} argument in the \code{raster} package.
#' @param fun function or character. Determines what values to assign to cells with multiple spatial features, similar to the \code{fun} argument in \code{\link[raster]{rasterize}}.
#' @param nf integer. Specifies the number of specialization axes to keep after transformation.
#' @param scale logical. If \code{TRUE} then the values of the Raster* object will
#' be centered and scaled. Depending on the resolution of the climate data and the
#' extent of the study area, this can be quite time consuming. If running this
#' function for multiple species, it is recommended that the climate data be scaled beforehand.
#' @param ... Additonal arguments for rasterizing as for \code{\link[raster]{rasterize}}.
#' @details Hirzel et al. (2002) defined the overall marginality $M$ as a standardized distance between the centroid of the species' niche and the global niche. Basille and Calenge (2008), however, defined $M$ as $$M = m^Tm$$. The default \code{mar.type} reflects Basille and Calenge's definition so that results will by default agree with those calculated using Basille and Calenge's \code{adehabitatHS} package.
#' @return Returns an S4 object of class \code{enfa} with the following components:
#' @return mf marginality factor. Vector that describes the location of the species Hutchinsonian niche.
#' @return marginality. Standardized magnitude of the marginality factor.
#' @return s specialization factors.  of specializations.
#' @return speciality specialization. Vector
#' @return ras Raster* object of transformed climate values, with number of layers equal to nf + 1.
#' @export
#'
#'

setGeneric("enfa", function(x, s.dat, ...) {
  standardGeneric("enfa")
})

#' @rdname enfa
setMethod("enfa",
          signature(x = "RasterBrick", s.dat = "SpatialPolygonsDataFrame"),
          function(x, s.dat, field,
                   nf = 1, scale = FALSE, ...){
            call <- match.call()
            if(!identicalCRS(x, s.dat)) {stop("spatial projections of environmental and species data do not match")}
            if(scale == TRUE) {x <- raster::scale(x)}

            x.mask <- mask(x, s.dat)
            gpres <- which(!is.na(values(x[[1]])))
            s.dat.ras <- rasterize(s.dat, x, field = field)
            pres <- which(!is.na(values(s.dat.ras)))
            prb <- c(s.dat.ras[pres])
            pr <- prb/sum(prb)

            small <- canProcessInMemory(x)
            if(small){
              dat <- values(x)[gpres, ]
              pres.dat <- values(x.mask)[pres, ]
              center <- colMeans(dat)
              Z <- sweep(dat, 2, center)
              S <- sweep(pres.dat, 2, center)
              nZ <- nrow(Z)
              nS <- nrow(S)
              Rg <- crossprod(Z, Z/nZ)
              Rs <- crossprod(S,S/nS)

            } else {
              center <- cellStats(x, mean)
              Z <- calc(x, fun = function(x) {x - center})
              S <- calc(x.mask, fun = function(x) {x - center})
              nZ <- nlayers(Z)
              nS <- nlayers(S)
              Rg <- covmat(Z, sample = F, ...)
              Rs <- covmat(S, ...)
            }

            #row.w<-rep(1,nrow(dat))/nrow(dat)
            #col.w <- rep(1,ncol(dat))
            #center <- colMeans(dat)
            mar <- colMeans(S)/sum(prb)
            eigRs <- eigen(Rs)
            keep <- (eigRs$values > 1e-09)
            Rs12 <- eigRs$vectors[, keep] %*% diag(eigRs$values[keep]^(-0.5)) %*% t(eigRs$vectors[, keep])
            W <- Rs12 %*% Rg %*% Rs12
            z <- Rs12 %*% mar
            y <- z/sqrt(sum(z^2))
            H <- (diag(ncol(Z)) - y %*% t(y)) %*% W %*% (diag(ncol(Z)) - y %*% t(y))
            s <- eigen(H)$values
            spec <- sqrt(sum(s[1:nf]))/nf
            v <- eigen(H)$vectors
            if (nf <= 0 | nf > (ncol(Z) - 1)) {nf <- 1}
            co <- matrix(nrow = ncol(Z), ncol = nf + 1)
            u <- (Rs12 %*% v)[, 1:nf]
            norw <- sqrt(diag(t(as.matrix(u)) %*% as.matrix(u)))
            co[, 2:(nf + 1)] <- sweep(as.matrix(u), 2, norw, "/")
            co[, 1] <- mar
            if(return_values == TRUE){
              li <- S %*% co
              ras <- raster::subset(x.mask, 1:(nf+1))
              values(ras)[pres, ] <- li
              names(ras) <- c("Marg", paste0("Spec", (1:nf)))
            } else {ras <- NA}
            co <- as.data.frame(co)
            names(co) <- c("Marg", paste0("Spec", (1:nf)))
            row.names(co) <- dimnames(x)[[2]]
            enfa <- methods::new("enfa", call = call, mf = mar, marginality = m, sf = sf, sensitivity = sens, p.spec = s.p, co = co, cov = Rs, present = length(pres), ras = s.ras)
            return(enfa)
          }
)

#' @rdname enfa
setMethod("enfa",
          signature(x = "RasterBrick", s.dat = "SpatialPoints"),
          function(x, s.dat, field,
                   nf = 1, scale = FALSE, fun = "count", ...){
            #call <- match.call()
            if(!identicalCRS(x, s.dat)) {stop("spatial projections of environmental and species data do not match")}
            if(scale == TRUE) {x <- raster::scale(x)}

            x.mask <- mask(x, s.dat)
            gpres <- which(!is.na(values(x[[1]])))
            s.dat.ras <- rasterize(s.dat, x.mask, field = field, fun = fun)
            pres <- which(!is.na(values(s.dat.ras)))
            prb <- c(s.dat.ras[pres])
            pr <- prb/sum(prb)

            small <- canProcessInMemory(x)
            if(small){
              dat <- values(x)[gpres, ]
              pres.dat <- values(x.mask)[pres, ]
              center <- colMeans(dat)
              Z <- sweep(dat, 2, center)
              S <- sweep(pres.dat, 2, center)
              nZ <- nrow(Z)
              nS <- nrow(S)
              Rg <- crossprod(Z, Z/nZ)
              Rs <- crossprod(S,S/nS)

            } else {
              center <- cellStats(x, mean)
              Z <- calc(x, fun = function(x) {x - center})
              S <- calc(x.mask, fun = function(x) {x - center})
              nZ <- nlayers(Z)
              nS <- nlayers(S)
              Rg <- covmat(Z, sample = F, ...)
              Rs <- covmat(S, ...)
            }

            #row.w<-rep(1,nrow(dat))/nrow(dat)
            #col.w <- rep(1,ncol(dat))
            #center <- colMeans(dat)
            mar <- colMeans(S)/sum(prb)
            eigRs <- eigen(Rs)
            keep <- (eigRs$values > 1e-09)
            Rs12 <- eigRs$vectors[, keep] %*% diag(eigRs$values[keep]^(-0.5)) %*% t(eigRs$vectors[, keep])
            W <- Rs12 %*% Rg %*% Rs12
            z <- Rs12 %*% mar
            y <- z/sqrt(sum(z^2))
            H <- (diag(ncol(Z)) - y %*% t(y)) %*% W %*% (diag(ncol(Z)) - y %*% t(y))
            s <- eigen(H)$values
            spec <- sqrt(sum(s[1:nf]))/nf
            v <- eigen(H)$vectors
            if (nf <= 0 | nf > (ncol(Z) - 1)) {nf <- 1}
            co <- matrix(nrow = ncol(Z), ncol = nf + 1)
            u <- (Rs12 %*% v)[, 1:nf]
            norw <- sqrt(diag(t(as.matrix(u)) %*% as.matrix(u)))
            co[, 2:(nf + 1)] <- sweep(as.matrix(u), 2, norw, "/")
            co[, 1] <- mar
            if(return_values == TRUE){
              li <- S %*% co
              ras <- raster::subset(x.mask, 1:(nf+1))
              values(ras)[pres, ] <- li
              names(ras) <- c("Marg", paste0("Spec", (1:nf)))
            } else {ras <- NA}
            co <- as.data.frame(co)
            names(co) <- c("Marg", paste0("Spec", (1:nf)))
            row.names(co) <- dimnames(x)[[2]]
            enfa <- methods::new("enfa", call = call, mf = mar, marginality = m, sf = sf, sensitivity = sens, p.spec = s.p, co = co, cov = Rs, present = length(pres), ras = s.ras)
            return(enfa)
          }
)

#' #' @rdname enfa
#' setMethod("enfa",
#'           signature(climdat = "GLenfa", s.dat = "SpatialPolygonsDataFrame"),
#'           function(climdat, s.dat, field, nf = 1){
#'
#'             call <- match.call()
#'
#'             if(!identicalCRS(climdat@global_ras, s.dat)){
#'               stop("climate and species projections do not match")
#'             }
#'
#'
#'             if(length(raster::intersect(extent(climdat@global_ras), extent(s.dat)))==0){
#'               stop("climate and species data to not overlap")
#'             }
#'
#'             rr <- crop(climdat@global_ras,extent(s.dat))
#'             s.dat.ras <- rasterize(s.dat,rr,field=field)
#'             pres <- which(!is.na(values(s.dat.ras)) & !is.na(values(rr[[1]])))
#'             pres.dat <- values(rr)[pres,]
#'             S <- sweep(pres.dat, 2, climdat@center)
#'             rZ <- climdat@ncells
#'             cZ <- nlayers(climdat@global_ras)
#'             rS <- nrow(S)
#'             mar <- colSums(S)/rS
#'             m <- norm(mar,"2")/1.96
#'             Rg <- climdat@cov
#'             Rs <- crossprod(S,S/rS)
#'             eigRs <- eigen(Rs)
#'             keep <- (eigRs$values > 1e-09)
#'             Rs12 <- eigRs$vectors[, keep] %*% diag(eigRs$values[keep]^(-0.5)) %*% t(eigRs$vectors[, keep])
#'             W <- Rs12 %*% Rg %*% Rs12
#'             z <- Rs12 %*% mar
#'             y <- z/sqrt(sum(z^2))
#'             H <- (diag(cZ) - y %*% t(y)) %*% W %*% (diag(cZ) - y %*% t(y))
#'             s <- eigen(H)$values[-cZ]
#'             s.p <- abs(s)/sum(abs(s))
#'             s.p[1] <- sum(diag(W)) - sum(diag(H))
#'             spec <- sqrt(sum(s))/cZ
#'             v <- eigen(H)$vectors
#'             if (nf<=0 | nf>(cZ-1)){nf <- 1}
#'             co <- matrix(nrow = cZ, ncol = nf + 1)
#'             u <- (Rs12 %*% v)[, 1:nf]
#'             norw <- sqrt(diag(t(as.matrix(u)) %*% as.matrix(u)))
#'             co[, 2:(nf + 1)] <- sweep(as.matrix(u), 2, norw, "/")
#'             co[, 1] <- mar
#'             ras<-brick(rr,nl=nf+1)
#'             values(ras)[pres,]<- S %*% co
#'             names(ras) <- c("Marg", paste0("Spec", (1:nf)))
#'             co <- as.data.frame(co)
#'             names(co) <- c("Marg", paste0("Spec", (1:nf)))
#'             row.names(co) <- dimnames(pres.dat)[[2]]
#'             cnfa<-methods::new("cnfa", call = call, mf = mar, marginality = m, s = s, specialization = spec, spec.account = s.p, co = co, species_ras = ras, present = length(pres))
#'             return(cnfa)
#'           }
#' )
