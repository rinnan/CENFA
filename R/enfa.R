#' Ecological-Niche Factor Analysis
#'
#' Performs ecological-niche factor analysis using environmental raster data and
#' species presence data in shapefile.
#'
#' @aliases print.enfa, show.enfa
#' @param x Raster* object, typically a brick or stack of ecological raster layers
#' @param s.dat matrix, SpatialPolygons*, or SpatialPoints* object detailing
#'   species presence or abundance
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
#' @param ... Additional arguments for rasterizing as for \code{\link[raster]{rasterize}}.
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
#' \dontrun{
#' mod2 <- enfa(x = glc, s.dat = ABPR, field = "CODE")
#' all.equal(mf(mod1), mf(mod2))}
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
#'   \item{present}{Number of raster cells in which species is present.}
#'   \item{ras}{Raster* object of transformed environmental values.}
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
          signature(x = "Raster", s.dat = "Spatial"),
          function(x, s.dat, field, fun = "last", scale = TRUE, filename = "", ...){

            call <- match.call()

            if (! inherits(x, 'Raster')) stop('"x" should be a "Raster*" object')
            if (! inherits(s.dat, c('SpatialPolygons', 'SpatialPoints'))) stop('"s.dat" should be a "SpatialPolygons*" or "SpatialPoints*" object')
            if(!identicalCRS(x, s.dat)) stop("projections do not match")
            if(is.null(intersect(extent(x), extent(s.dat)))) stop("climate and species data do not overlap")
            if(union(extent(x), extent(s.dat)) != extent(x)) stop("extent of species data not contained within extent of climate data")

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
            sf <- eigen(H)$values
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
              #setValues(s.ras, S %*% co, index = pres)
              names(s.ras) <- nm
            } else{
              cat("\nCreating factor rasters...")
              s.ras <- .calc(x.mask, function(x) {x %*% co}, forceapply = T, filename = filename, names = nm, ...)
            }
            colnames(co) <- c("Marg", paste0("Spec", (1:(cZ-1))))
            rownames(co) <- names(x)

            enfa <- methods::new("enfa", call = call, mf = mar, marginality = m, sf = sf, specialization = spec,
                                 p.spec = s.p, co = co, cov = Rs, present = length(pres), ras = s.ras)
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
