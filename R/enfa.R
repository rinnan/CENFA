#' Ecological-Niche Factor Analysis
#'
#' Performs ecological-niche factor analysis using environmental raster data and species presence data in shapefile.
#'
#' @aliases print.enfa, show.enfa
#' @param ecodat Raster* object, typically a brick or stack of ecological raster layers
#' @param speciesdat SpatialPolygons* object detailing species presence or abundance
#' @param field field of \code{speciesdat} that specifies presence or abundance. This is equivalent to the \code{field} argument in the \code{raster} package.
#' @param nf integer. Specifies the number of specialization axes to keep after transformation.
#' @param scale logical. If \code{TRUE} then the values of the Raster* object will
#' be centered and scaled. Depending on the resolution of the climate data and the
#' extent of the study area, this can be quite time consuming. If running this
#' function for multiple species, it is recommended that the climate data be scaled beforehand.
#' @param mar.type character. Choices are "Hirzel" or "Basille". See details.
#' @details Hirzel et al. (2002) defined the overall marginality $M$ as a standardized distance between the centroid of the species' niche and the global niche, given by $$M = \frac{\sqrt{\sum\limits_{j=1}^P m_j^2}}{1.96}.$$ Basille and Calenge (2008), however, defined $M$ as $$M = \mathbf{m}^T\mathbf{m}$$. The default \code{mar.type} reflects Basille and Calenge's definition so that results will by default agree with those calculated using Basille and Calenge's \code{adehabitatHS} package.
#' @return Returns an S4 object of class \code{enfa} with the following components:
#' @param mf marginality factor. Vector that describes the location of the species Hutchinsonian niche.
#' @param marginality. Standardized magnitude of the marginality factor.
#' @param s specialization factors. Matrix of specializations.
#' @param speciality specialization. Mahalanobis distance of something.
#' @param ras Raster* object of transformed climate values, with number of layers equal to nf + 1.
#' @export
#'
#'

setGeneric("enfa", function(ecodat, speciesdat, ...) {
  standardGeneric("enfa")
})

#' @rdname enfa
setMethod("enfa",
          signature(ecodat = "Raster", speciesdat = "SpatialPolygonsDataFrame"),
          function(ecodat, speciesdat, field,
                   nf = 1, scale = FALSE, mar.type = "Basille", return_values = FALSE, ...){
            call <- match.call()
            if(!identicalCRS(ecodat, speciesdat)) {stop("spatial projections of environmental and species data do not match")}
            if(scale == TRUE) {ecodat <- raster::scale(ecodat)}

            ecodat.mask <- mask(ecodat, speciesdat)
            gpres <- which(!is.na(values(ecodat[[1]])))
            speciesdat.ras <- rasterize(speciesdat, ecodat, field = field)
            pres <- which(!is.na(values(speciesdat.ras)))
            prb <- c(speciesdat.ras[pres])
            pr <- prb/sum(prb)

            small <- canProcessInMemory(ecodat)
            if(small){
              dat <- values(ecodat)[gpres, ]
              pres.dat <- values(ecodat.mask)[pres, ]
              center <- colMeans(dat)
              Z <- sweep(dat, 2, center)
              S <- sweep(pres.dat, 2, center)
              nZ <- nrow(Z)
              nS <- nrow(S)
              Rg <- crossprod(Z, Z/nZ)
              Rs <- crossprod(S,S/nS)

            } else {
              center <- cellStats(ecodat, mean)
              Z <- calc(ecodat, fun = function(x) {x - center})
              S <- calc(ecodat.mask, fun = function(x) {x - center})
              nZ <- nlayers(Z)
              nS <- nlayers(S)
              Rg <- covmat(Z, sample = F, ...)
              Rs <- covmat(S, ...)
            }

            #row.w<-rep(1,nrow(dat))/nrow(dat)
            #col.w <- rep(1,ncol(dat))
            #center <- colMeans(dat)
            mar <- colMeans(S)/sum(prb)
            if(mar.type == "Basille") m <- c(t(mar) %*% mar)
            if(mar.type == "Hirzel")  m <- norm(mar, "2")/1.96
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
              ras <- raster::subset(ecodat.mask, 1:(nf+1))
              values(ras)[pres, ] <- li
              names(ras) <- c("Marg", paste0("Spec", (1:nf)))
            } else {ras <- NA}
            co <- as.data.frame(co)
            names(co) <- c("Marg", paste0("Spec", (1:nf)))
            row.names(co) <- dimnames(ecodat)[[2]]
            enfa <- methods::new("enfa", call = call, mf = mar, marginality = m, s = s, specialization = spec, spec.account = s.p, co = co, ras = ras, present = length(pres))
            return(enfa)
          }
)

#' @rdname enfa
setMethod("enfa",
          signature(ecodat = "Raster", speciesdat = "SpatialPoints"),
          function(ecodat, speciesdat, field,
                   nf = 1, scale = FALSE, fun = "count", mar.type = "Basille", return_values = FALSE, ...){
            #call <- match.call()
            if(!identicalCRS(ecodat, speciesdat)) {stop("spatial projections of environmental and species data do not match")}
            if(scale == TRUE) {ecodat <- raster::scale(ecodat)}

            ecodat.mask <- mask(ecodat, speciesdat)
            gpres <- which(!is.na(values(ecodat[[1]])))
            speciesdat.ras <- rasterize(speciesdat, ecodat.mask, field = field, fun = fun)
            pres <- which(!is.na(values(speciesdat.ras)))
            prb <- c(speciesdat.ras[pres])
            pr <- prb/sum(prb)

            small <- canProcessInMemory(ecodat)
            if(small){
              dat <- values(ecodat)[gpres, ]
              pres.dat <- values(ecodat.mask)[pres, ]
              center <- colMeans(dat)
              Z <- sweep(dat, 2, center)
              S <- sweep(pres.dat, 2, center)
              nZ <- nrow(Z)
              nS <- nrow(S)
              Rg <- crossprod(Z, Z/nZ)
              Rs <- crossprod(S,S/nS)

            } else {
              center <- cellStats(ecodat, mean)
              Z <- calc(ecodat, fun = function(x) {x - center})
              S <- calc(ecodat.mask, fun = function(x) {x - center})
              nZ <- nlayers(Z)
              nS <- nlayers(S)
              Rg <- covmat(Z, sample = F, ...)
              Rs <- covmat(S, ...)
            }

            #row.w<-rep(1,nrow(dat))/nrow(dat)
            #col.w <- rep(1,ncol(dat))
            #center <- colMeans(dat)
            mar <- colMeans(S)/sum(prb)
            ifelse(mar.type = "Basille",
                   m <- t(mar) %*% mar,
                   m <- norm(mar,"2")/1.96)
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
              ras <- raster::subset(ecodat.mask, 1:(nf+1))
              values(ras)[pres, ] <- li
              names(ras) <- c("Marg", paste0("Spec", (1:nf)))
            } else {ras <- NA}
            co <- as.data.frame(co)
            names(co) <- c("Marg", paste0("Spec", (1:nf)))
            row.names(co) <- dimnames(ecodat)[[2]]
            enfa <- methods::new("enfa", call = call, mf = mar, marginality = m, s = s, specialization = spec, co = co, ras = ras, present = length(pres))
            return(enfa)
          }
)

#' #' @rdname enfa
#' setMethod("enfa",
#'           signature(climdat = "GLenfa", speciesdat = "SpatialPolygonsDataFrame"),
#'           function(climdat, speciesdat, field, nf = 1){
#'
#'             call <- match.call()
#'
#'             if(!identicalCRS(climdat@global_ras, speciesdat)){
#'               stop("climate and species projections do not match")
#'             }
#'
#'
#'             if(length(raster::intersect(extent(climdat@global_ras), extent(speciesdat)))==0){
#'               stop("climate and species data to not overlap")
#'             }
#'
#'             rr <- crop(climdat@global_ras,extent(speciesdat))
#'             speciesdat.ras <- rasterize(speciesdat,rr,field=field)
#'             pres <- which(!is.na(values(speciesdat.ras)) & !is.na(values(rr[[1]])))
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
