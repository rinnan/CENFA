#' Climate-Niche Factor Analysis
#'
#' Performs climate-niche factor analysis using climate raster data and species presence data in shapefile.
#'
#' @aliases print.cnfamad, show.cnfamad
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

setGeneric("cnfamad", function(climdat, speciesdat,...) {
  standardGeneric("cnfamad")
})

#' @rdname cnfamad
setMethod("cnfamad",
          signature(climdat = "RasterBrick", speciesdat = "SpatialPolygonsDataFrame"),
          function(climdat, speciesdat, field,
                   nf = 1, scale = FALSE){
            call <- match.call()
            if(!identicalCRS(climdat,speciesdat)) {stop("climate and species projections do not match")}
            if(length(raster::intersect(extent(climdat),extent(speciesdat)))==0) {stop("climate and species data to not overlap")}
            if(scale==TRUE) {
              meds<-cellStats(climdat,'median')
              mads<-cellStats(climdat,'mad')
              climdat<-raster::scale(climdat,center=meds,scale=mads)
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
            center <- apply(dat,2,median,na.rm=T)
            Z <- sweep(dat, 2, center)
            S <- sweep(pres.dat, 2, center)
            nZ <- nrow(Z)
            nS <- nrow(S)
            mar <- apply(pres.dat,2,median,na.rm=T)
            m <- norm(mar,"2")/(1.96*1.4826)
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
            cnfa<-methods::new("cnfamad", call = call, mf = mar, marginality = m, s = s, specialization = spec, co = co, species_ras = ras, present = length(pres))
            return(cnfa)
          }
)

#' @rdname cnfamad
setMethod("cnfamad",
          signature(climdat = "GLcnfamad", speciesdat = "SpatialPolygonsDataFrame"),
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
            cnfa<-methods::new("cnfamad", call = call, mf = mar, marginality = m, s = s, specialization = spec, co = co, species_ras = ras, present = length(pres))
            return(cnfa)
          }
)

cnfa.mad<-function(climdat, speciesdat, field, nf = 1, scale = FALSE){
  call <- match.call()
  if(!identicalCRS(climdat,speciesdat)) {stop("climate and species projections do not match")}
  if(length(raster::intersect(extent(climdat),extent(speciesdat)))==0) {stop("climate and species data to not overlap")}
  if(scale==TRUE) {
    meds<-cellStats(climdat,'median')
    mads<-cellStats(climdat,'mad')
    climdat<-raster::scale(climdat,center=meds,scale=mads)
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
  center <- apply(dat,2,median,na.rm=T)
  Z <- sweep(dat, 2, center)
  S <- sweep(pres.dat, 2, center)
  nZ <- nrow(Z)
  nS <- nrow(S)
  mar <- apply(pres.dat,2,median,na.rm=T)
  m <- norm(mar,"2")/(1.96*1.4826)
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
  cnfa<-methods::new("cnfa.mad", call = call, mf = mar, marginality = m, s = s, specialization = spec, co = co, species_ras = ras, present = length(pres))
  return(cnfa)
}
