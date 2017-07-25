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
          signature(hist.dat = "RasterBrick", fut.dat = "RasterBrick", species.dat = "cnfa"),
          function(hist.dat, fut.dat, species.dat, depart.ras = TRUE){

            call <- match.call()
            if(!identicalCRS(hist.dat, species.dat@species_ras)) {stop("historical climate and species projections do not match")}
            if(!identicalCRS(hist.dat, fut.dat))     {stop("historical and future climate projections do not match")}
            if(!identicalCRS(fut.dat, species.dat@species_ras)) {stop("future climate and species projections do not match")}
            if(length(raster::intersect(extent(hist.dat), extent(species.dat@species_ras)))==0) {stop("climate and species data to not overlap")}
            #             if(scale == TRUE) {
            #               hist.dat <- raster::scale(hist.dat)
            #               means <- cellStats(hist.dat, mean)
            #               sds <- cellStats(fut.dat, sd)
            #               fut.dat <- scale(fut.dat, center = means, scale = sds)
            #             }

            speciesdat.ras <- species.dat@species_ras
            pres <- which(!is.na(values(speciesdat.ras)))
            Ns <- length(pres)
            p.dat <- values(hist.dat)[pres,]
            f.dat <- values(fut.dat)[pres,]
            z_ij <- p.dat %*% as.matrix(species.dat@co[,-1])
            f_ij <- f.dat %*% as.matrix(species.dat@co[,-1])
            d_ij <- f_ij - z_ij
            d <- sqrt(rowSums(d_ij^2))
            D <- 1/(1.96*Ns) * sum(d, na.rm = T)

            if(depart.ras == T){
              ras <- speciesdat.ras
              values(ras)[pres] <- d
            }
            else ras <- NA

            depart <- methods::new("departure", call = call, departure = D, distances = d, departure_ras = ras, present = Ns)
            return(depart)
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
