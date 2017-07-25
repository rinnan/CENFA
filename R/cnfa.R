#' Climate-Niche Factor Analysis
#'
#' Performs climate-niche factor analysis using climate raster data and species presence data in shapefile.
#'
#' @aliases print.cnfa, show.cnfa
#' @param climdat Raster* object, typically a brick or stack of climate raster layers
#' @param speciesdat SpatialPolygons* object detailing species presence or abundance
#' @param field field of \code{speciesdat} that specifies presence or abundance. This is equivalent to the \code{field} argument in the \code{raster} package.
#' @param nf integer. Specifies the number of specialization axes to keep after transformation.
#' @param scale logical. If \code{TRUE} then the values of the Raster* object will
#' be centered and scaled. Depending on the resolution of the climate data and the
#' extent of the study area, this can be quite time consuming. If running this
#' function for multiple species, it is recommended that the climate data be scaled beforehand.
#' @param sp.prj character. Spatial projection of species data.
#' @return Returns an S4 object of class \code{cnfa} with the following components:
#' @return call original function call
#' @return tab a dataframe
#' @export
#'
#'

setGeneric("cnfa", function(climdat, speciesdat,...) {
  standardGeneric("cnfa")
})

#' @rdname cnfa
setMethod("cnfa",
          signature(climdat = "RasterBrick", speciesdat = "SpatialPolygonsDataFrame"),
          function(climdat, speciesdat, field,
                   nf = 1, scale = FALSE){
            call <- match.call()
            if(!identicalCRS(climdat,speciesdat)) {stop("climate and species projections do not match")}
            if(length(raster::intersect(extent(climdat),extent(speciesdat)))==0) {stop("climate and species data to not overlap")}
            if(scale==TRUE) {
              climdat<-raster::scale(climdat)
            }
            gpres<-which(!is.na(values(climdat[[1]])))
            dat<-values(climdat)[gpres,]
            speciesdat.ras<-rasterize(speciesdat,climdat,field=field)
            pres<-which(!is.na(values(speciesdat.ras)))
            #prb<-values(speciesdat.ras)[pres]
            pres.dat<-values(climdat)[pres,]
            #pr <- prb/sum(prb)
            #row.w<-rep(1,nrow(dat))/nrow(dat)
            #col.w <- rep(1,ncol(dat))
            center <- colMeans(dat)
            Z <- sweep(dat, 2, center)
            S <- sweep(pres.dat, 2, center)
            nZ <- nrow(Z)
            nS <- nrow(S)
            mar <- colSums(S)/nS
            m <- norm(mar,"2")/1.96
            Rg <- crossprod(Z, Z/nZ)
            Rs <- crossprod(S,S/nS)
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
            if (nf <= 0 | nf > (ncol(Z) - 1)) {nf <- 1}
            co <- matrix(nrow = ncol(Z), ncol = nf + 1)
            u <- (Rs12 %*% v)[, 1:nf]
            norw <- sqrt(diag(t(as.matrix(u)) %*% as.matrix(u)))
            co[, 2:(nf + 1)] <- sweep(as.matrix(u), 2, norw, "/")
            co[, 2:(nf + 1)] <- as.matrix(u)#sweep(as.matrix(u), 2, norw, "/")
            co[, 1] <- mar
            ras<-brick(climdat,nl=nf+1)
            ras<-crop(ras,extent(speciesdat))
            ss<-crop(speciesdat.ras,extent(speciesdat))
            pres.c<-which(!is.na(values(ss)))
            values(ras)[pres.c,]<- S %*% co
            #global_ras<-raster::subset(climdat,1:(nf+1))
            #values(global_ras)[gpres,]<- Z %*% co
            names(ras) <- c("Marg", paste0("Spec", (1:nf)))
            co <- as.data.frame(co)
            names(co) <- c("Marg", paste0("Spec", (1:nf)))
            row.names(co) <- dimnames(dat)[[2]]
            co<-as.data.frame(co[order(abs(co$Marg),decreasing = T),])
            cnfa<-methods::new("cnfa", call = call, mf = mar, marginality = m, s = s, specialization = spec, spec.account = s.p, co = co, species_ras = ras, present = length(pres))
            return(cnfa)
          }
)

#' @rdname cnfa
setMethod("cnfa",
          signature(climdat = "GLcnfa", speciesdat = "SpatialPolygonsDataFrame"),
          function(climdat, speciesdat, field, nf=1){

            call <- match.call()

            if(!identicalCRS(climdat@global_ras, speciesdat)){
              stop("climate and species projections do not match")
            }


            if(length(raster::intersect(extent(climdat@global_ras), extent(speciesdat)))==0){
              stop("climate and species data to not overlap")
            }

            rr <- crop(climdat@global_ras,extent(speciesdat))
            speciesdat.ras <- rasterize(speciesdat,rr,field=field)
            pres <- which(!is.na(values(speciesdat.ras)) & !is.na(values(rr[[1]])))
            pres.dat <- values(rr)[pres,]
            S <- sweep(pres.dat, 2, climdat@center)
            rZ <- climdat@ncells
            cZ <- nlayers(climdat@global_ras)
            rS <- nrow(S)
            mar <- colSums(S)/rS
            m <- norm(mar,"2")/1.96
            Rg <- climdat@cov
            Rs <- crossprod(S,S/rS)
            eigRs <- eigen(Rs)
            keep <- (eigRs$values > 1e-09)
            Rs12 <- eigRs$vectors[, keep] %*% diag(eigRs$values[keep]^(-0.5)) %*% t(eigRs$vectors[, keep])
            W <- Rs12 %*% Rg %*% Rs12
            z <- Rs12 %*% mar
            y <- z/sqrt(sum(z^2))
            H <- (diag(cZ) - y %*% t(y)) %*% W %*% (diag(cZ) - y %*% t(y))
            s <- eigen(H)$values[-cZ]
            s.p <- abs(s)/sum(abs(s))
            s.p[1] <- sum(diag(W)) - sum(diag(H))
            spec <- sqrt(sum(s))/cZ
            v <- eigen(H)$vectors
            if (nf<=0 | nf>(cZ-1)){nf <- 1}
            co <- matrix(nrow = cZ, ncol = nf + 1)
            u <- (Rs12 %*% v)[, 1:nf]
            norw <- sqrt(diag(t(as.matrix(u)) %*% as.matrix(u)))
            co[, 2:(nf + 1)] <- sweep(as.matrix(u), 2, norw, "/")
            co[, 1] <- mar
            ras<-brick(rr,nl=nf+1)
            values(ras)[pres,]<- S %*% co
            names(ras) <- c("Marg", paste0("Spec", (1:nf)))
            co <- as.data.frame(co)
            names(co) <- c("Marg", paste0("Spec", (1:nf)))
            row.names(co) <- dimnames(pres.dat)[[2]]
            cnfa<-methods::new("cnfa", call = call, mf = mar, marginality = m, s = s, specialization = spec, spec.account = s.p, co = co, species_ras = ras, present = length(pres))
            return(cnfa)
          }
)


#' @rdname cnfa
setMethod("cnfa",
          signature(climdat = "Raster", speciesdat = "SpatialPoints"),
          function(climdat, speciesdat, field, fun='count',
                   nf = 1,scale=FALSE){
            call <- match.call()
            if(!identicalCRS(climdat,speciesdat)) stop("climate and species projections do not match")
            if(scale==TRUE){climdat<-scale(climdat)}
            dat<-na.omit(values(climdat))
            #ext<-intersect(extent(climdat),extent(speciesdat))
            #climdat.crop<-crop(climdat,ext) %>% mask(.,speciesdat)
            climdat.crop<- mask(climdat,speciesdat)

            climdat.present<-na.omit(values(climdat.crop))
            #speciesdat.ras<-rasterize(speciesdat,climdat[[1]],field=field,background=0) %>%mask(.,climdat[[1]])
            speciesdat.ras<-rasterize(speciesdat,climdat.crop[[1]],field=field,fun=fun)# %>%mask(.,climdat[[1]])
            prb<-na.omit(values(speciesdat.ras))
            #prb <- pr
            pr <- prb/sum(prb)
            row.w<-rep(1,nrow(dat))/nrow(dat)
            col.w <- rep(1,ncol(dat))
            n <- nrow(dat)
            center<-colMeans(dat)
            Z <- sweep(dat, 2, center)
            A<-sweep(climdat.present,2,center)/sum(prb) #equiv to DpZ
            #Ze <- sweep(Z, 2, sqrt(col.w), "*")
            #DpZ <- apply(Z, 2, function(x) x * pr)
            #mar <- apply(Z, 2, function(x) sum(x * pr))
            mar<-colSums(A)
            #me <- mar
            #Se <- crossprod(Z, DpZ)
            Se<-crossprod(climdat.present,A)
            #Ge <- crossprod(Ze, apply(Ze, 2, function(x) x * row.w))
            Ge <- crossprod(Z, Z/n)
            eS <- eigen(Se)
            kee <- (eS$values > 1e-09)
            S12 <- eS$vectors[, kee] %*% diag(eS$values[kee]^(-0.5)) %*% t(eS$vectors[, kee])
            W <- S12 %*% Ge %*% S12
            x <- S12 %*% mar
            b <- x/sqrt(sum(x^2))
            H <- (diag(ncol(Z)) - b %*% t(b)) %*% W %*% (diag(ncol(Z)) - b %*% t(b))
            s <- eigen(H)$values[-ncol(Z)]
            # if (scannf) {
            #   barplot(s)
            #   cat("Select the number of specialization axes: ")
            #   nf <- as.integer(readLines(n = 1))
            # }
            if (nf <= 0 | nf > (ncol(Z) - 1))
              nf <- 1
            co <- matrix(nrow = ncol(Z), ncol = nf + 1)
            tt <- data.frame((S12 %*% eigen(H)$vectors)[, 1:nf])
            #ww <- apply(tt, 2, function(x) x/sqrt(col.w))
            norw <- sqrt(diag(t(as.matrix(tt)) %*% as.matrix(tt)))
            co[, 2:(nf + 1)] <- sweep(as.matrix(tt), 2, norw, "/")
            m <- mar
            co[, 1] <- m/sqrt(sum(m^2))
            m <- sum(m^2)
            li <- Z %*% co
            co <- as.data.frame(co)
            li <- as.data.frame(li)
            names(co) <- c("Mar", paste0("Spe", (1:nf)))
            row.names(co) <- dimnames(dat)[[2]]
            names(li) <- c("Mar", paste0("Spe", (1:nf)))
            cnfa <- list(call = call, pr = prb,
                         nf = nf, m = m, s = s, li = li,
                         co = co, mar = mar)
            class(cnfa) <- "cnfa"
            return(invisible(cnfa))
          }
)

#' @rdname cnfa
setMethod("cnfa",
          signature(climdat = "RasterBrick", speciesdat = "matrix"),
          function(climdat, speciesdat, field,
                   nf = 1,scale = FALSE,sp.prj){
            #call <- match.call()
            #             if(inherits(speciesdat, "data.frame")){
            #               speciesdat <- as.matrix(speciesdat)
            #             }
            speciesdat<-SpatialPointsDataFrame(coords=speciesdat, data=1,proj4string=CRS(sp.prj))
            if(!identicalCRS(climdat,speciesdat)) {stop("climate and species projections do not match")}
            if(length(intersect(extent(climdat),extent(speciesdat)))==0) {stop("climate and species data to not overlap")}
            if(scale==TRUE) {climdat<-raster::scale(climdat)}
            dat<-na.omit(values(climdat))
            climdat.crop<-mask(climdat,speciesdat)
            pres<-which(!is.na(values(climdat.crop[[1]])))
            pres.dat<-values(climdat.crop)[pres,]
            #pres.dat<-na.omit(values(climdat.crop))
            speciesdat.ras<-rasterize(speciesdat,climdat.crop,field=field)
            prb<-c(na.omit(values(speciesdat.ras)))
            pr <- prb/sum(prb)
            #row.w<-rep(1,nrow(dat))/nrow(dat)
            #col.w <- rep(1,ncol(dat))
            center <- colMeans(dat)
            Z <- sweep(dat, 2, center)
            S <- sweep(pres.dat, 2, center)
            nZ <- nrow(Z)
            nS <- nrow(S)
            mar <- colMeans(S)/sum(prb)
            m <- norm(mar,"2")/1.96
            Rg <- crossprod(Z, Z/nZ)
            Rs <- crossprod(S,S/nS)
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
            li <- S %*% co
            ras<-raster::subset(climdat.crop,1:(nf+1))
            values(ras)[pres,]<-li
            names(ras) <- c("Marg", paste0("Spec", (1:nf)))
            co <- as.data.frame(co)
            names(co) <- c("Marg", paste0("Spec", (1:nf)))
            row.names(co) <- dimnames(dat)[[2]]
            cnfa<-methods::new("cnfa", mf = mar, marginality = m, s = s, speciality = spec, co = co, ras = ras)
            return(cnfa)
          }
)
