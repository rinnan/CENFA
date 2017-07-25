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
#' @return Returns an S4 object of class \code{enfa} with the following components:
#' @param mf marginality factor. Vector that describes the location of the species Hutchinsonian niche.
#' @param marginality. Standardized magnitude of the marginality factor.
#' @param s specialization factors. Matrix of specializations.
#' @param speciality specialization. Mahalanobis distance of something.
#' @param ras Raster* object of transformed climate values, with number of layers equal to nf + 1.
#' @export
#'
#'

setGeneric("enfa", function(ecodat, speciesdat,...) {
  standardGeneric("enfa")
})

#' @rdname enfa
setMethod("enfa",
          signature(ecodat = "RasterBrick", speciesdat = "SpatialPolygonsDataFrame"),
          function(ecodat, speciesdat, field,
                   nf = 1,scale = FALSE,return_values = FALSE){
            #call <- match.call()
            if(!identicalCRS(ecodat,speciesdat)) {stop("spatial projections of environmental and species data do not match")}
            if(scale==TRUE) {ecodat<-raster::scale(ecodat)}
            dat<-na.omit(values(ecodat))
            ecodat.crop<-mask(ecodat,speciesdat)
            pres<-which(!is.na(values(ecodat.crop[[1]])))
            pres.dat<-values(ecodat.crop)[pres,]
            #pres.dat<-na.omit(values(ecodat.crop))
            speciesdat.ras<-rasterize(speciesdat,ecodat.crop,field=field)
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
            if(return_values == TRUE){
              li <- S %*% co
              ras<-raster::subset(ecodat.crop,1:(nf+1))
              values(ras)[pres,]<-li
              names(ras) <- c("Marg", paste0("Spec", (1:nf)))
            } else {ras<-NA}
            co <- as.data.frame(co)
            names(co) <- c("Marg", paste0("Spec", (1:nf)))
            row.names(co) <- dimnames(dat)[[2]]
            enfa<-methods::new("enfa", mf = mar, marginality = m, s = s, speciality = spec, co = co, ras = ras)
            return(enfa)
          }
)

#' @rdname enfa
setMethod("enfa",
          signature(ecodat = "Raster", speciesdat = "SpatialPoints"),
          function(ecodat, speciesdat, field, fun='count',
                   nf = 1,scale=FALSE){
            call <- match.call()
            if(!identicalCRS(ecodat,speciesdat)) stop("environmental and species projections do not match")
            if(scale==TRUE){ecodat<-scale(ecodat)}
            dat<-na.omit(values(ecodat))
            #ext<-intersect(extent(ecodat),extent(speciesdat))
            #ecodat.crop<-crop(ecodat,ext) %>% mask(.,speciesdat)
            ecodat.crop<- mask(ecodat,speciesdat)

            ecodat.present<-na.omit(values(ecodat.crop))
            #speciesdat.ras<-rasterize(speciesdat,ecodat[[1]],field=field,background=0) %>%mask(.,ecodat[[1]])
            speciesdat.ras<-rasterize(speciesdat,ecodat.crop[[1]],field=field,fun=fun)# %>%mask(.,ecodat[[1]])
            prb<-na.omit(values(speciesdat.ras))
            #prb <- pr
            pr <- prb/sum(prb)
            row.w<-rep(1,nrow(dat))/nrow(dat)
            col.w <- rep(1,ncol(dat))
            n <- nrow(dat)
            center<-colMeans(dat)
            Z <- sweep(dat, 2, center)
            A<-sweep(ecodat.present,2,center)/sum(prb) #equiv to DpZ
            #Ze <- sweep(Z, 2, sqrt(col.w), "*")
            #DpZ <- apply(Z, 2, function(x) x * pr)
            #mar <- apply(Z, 2, function(x) sum(x * pr))
            mar<-colSums(A)
            #me <- mar
            #Se <- crossprod(Z, DpZ)
            Se<-crossprod(ecodat.present,A)
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
            enfa <- list(call = call, pr = prb,
                         nf = nf, m = m, s = s, li = li,
                         co = co, mar = mar)
            class(enfa) <- "enfa"
            return(invisible(enfa))
          }
)
