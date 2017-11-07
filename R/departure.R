#' Climatic departure
#'
#' Calculates the climatic departure of a species using historical and future climate raster data and species presence data.
#'
#' @param x.hist Raster* object, typically a brick or stack of historical climate raster layers
#' @param x.fut  Raster* object, future climate values with the same layers as x.hist
#' @param s.dat SpatialPolygons* object detailing species presence or abundance
#' @param field field of \code{speciesdat} that specifies presence or abundance. This is equivalent to the \code{field} argument of \code{raster::rasterize}.
#' @param nf integer. Specifies the number of specialization axes to keep after transformation.
#' @param scale logical. If \code{TRUE} then the values of \code{x.hist} and \code{x.fut} will
#' be centered and scaled by the means and sds of the historical climate data. Depending on the resolution of the data and the
#' extent of the study area, this can be quite time consuming. If running this
#' function for multiple species, it is recommended that the climate data be scaled beforehand.
#' @param sp.prj character. Spatial projection of species data.
#' @return Returns an S4 object of class \code{cnfa} with the following slots:
#' @return call original function call
#' @return departure climatic departure D of species
#' @return departure_ras raster of distances d_i
#' @return present number of cells in which species is present
#'
#' @export
#'
#'

setGeneric("departure", function(x.hist, x.fut, s.dat, ...) {
  standardGeneric("departure")
})

setClass("departure", slots = list(call = "call", departure = "numeric", departure_ras = "Raster", present = "numeric"))

#' @rdname departure
setMethod("departure",
          signature(x.hist = "RasterBrick", x.fut = "RasterBrick", s.dat = "SpatialPolygonsDataFrame"),
          function(x.hist, x.fut, s.dat, field,
                   scale = FALSE, depart.ras = TRUE){
            call <- match.call()
            if(!identicalCRS(x.hist, s.dat)) {stop("historical climate and species projections do not match")}
            if(!identicalCRS(x.hist, x.fut))     {stop("historical and future climate projections do not match")}
            if(!identicalCRS(x.fut, s.dat)) {stop("future climate and species projections do not match")}
            if(length(raster::intersect(extent(x.hist), extent(s.dat))) == 0) {stop("climate and species data to not overlap")}
            if(scale) {
              x.hist <- raster::scale(x.hist)
              means <- cellStats(x.hist, mean)
              sds <- cellStats(x.fut, sd)
              x.fut <- scale(x.fut, center = means, scale = sds)
            }
            #gpres <- which(!is.na(values(x.hist[[1]])))
            #dat<-values(climdat)[gpres,]
            speciesdat.ras <- rasterize(s.dat, x.hist, field = field)
            pres <- which(!is.na(values(speciesdat.ras)))
            Ns <- length(pres)
            #prb<-values(speciesdat.ras)[pres]
            pres.dat <- values(x.hist)[pres,]
            d_ij <- values(x.fut)[pres,] - values(x.hist)[pres,]
            d <- sqrt(rowSums(d_ij^2))
            D <- 1/(1.96*Ns) * sum(d)

            depart <- methods::new("departure", call = call, departure = D, departure_ras = ras, present = Ns)
            return(depart)
          }
)

#' @rdname departure
setMethod("departure",
          signature(x.hist = "Raster", x.fut = "Raster", s.dat = "enfa"),
          function(x.hist, x.fut, s.dat, scale = FALSE){

            call <- match.call()
            sp.ras <- raster(s.dat)
            if(!identicalCRS(x.hist, sp.ras)) {stop("historical climate and species projections do not match")}
            if(!identicalCRS(x.hist, x.fut))     {stop("historical and future climate projections do not match")}
            if(!identicalCRS(x.fut, sp.ras)) {stop("future climate and species projections do not match")}
            if(length(raster::intersect(extent(x.hist), extent(sp.ras)))==0) {stop("climate and species data to not overlap")}
            if(scale) {
              center <- cellStats(x.hist, mean)
              sds <- cellStats(x.hist, sd)
              x.hist <- raster::scale(x.hist, center = center, scale = sds)
              x.fut <- raster::scale(x.fut, center = center, scale = sds)
            }
            x.hist <- crop(x.hist, sp.ras)
            x.fut <- crop(x.fut, sp.ras)

            pres <- which(!is.na(values(sp.ras[[1]])))
            small <- canProcessInMemory(x.hist, 8)
            if(small){
              z_ij <- values(x.hist)[pres,] %*% as.matrix(s.dat@co)
              f_ij <- values(x.fut)[pres,] %*% as.matrix(s.dat@co)
              d_ij <- (f_ij - z_ij)^2
              d <- sqrt(rowSums(d_ij))
              D <- 1/(1.96) * mean(d, na.rm = T)
              sp.ras[pres] <- d
            } else {
              x.hist <- calc(x.hist, fun = function(x) {x %*% as.matrix(s.dat@co)})
              x.fut <- calc(x.fut, fun = function(x) {x %*% as.matrix(s.dat@co)})
              d_ij <- (x.fut - x.hist)^2
              d <- calc(d_ij, fun = function(x) {sqrt(sum(x))})
              values(d)[!pres] <- NA
              D <- 1/(1.96) * cellStats(d, mean)
              sp.ras <- d
            }

            depart <- methods::new("departure", call = call, departure = D, departure_ras = sp.ras, present = length(pres))
            return(depart)
          }
)

#' @rdname departure
setMethod("departure",
          signature(x.hist = "RasterBrick", x.fut = "RasterBrick", s.dat = "cnfa"),
          function(x.hist, x.fut, s.dat, scale = FALSE){

            call <- match.call()
            sp.ras <- raster(s.dat)
            if(!identicalCRS(x.hist, sp.ras)) {stop("historical climate and species projections do not match")}
            if(!identicalCRS(x.hist, x.fut))     {stop("historical and future climate projections do not match")}
            if(!identicalCRS(x.fut, sp.ras)) {stop("future climate and species projections do not match")}
            if(length(raster::intersect(extent(x.hist), extent(sp.ras)))==0) {stop("climate and species data to not overlap")}
            if(scale) {
              center <- cellStats(x.hist, mean)
              sds <- cellStats(x.hist, sd)
              x.hist <- raster::scale(x.hist, center = center, scale = sds)
              x.fut <- raster::scale(x.fut, center = center, scale = sds)
            }
            x.hist <- crop(x.hist, sp.ras)
            x.fut <- crop(x.fut, sp.ras)
            Rs.inv <- solve(s.dat@s.cov)

            small <- canProcessInMemory(x.hist, 5)
            if(small){
              pres <- which(!is.na(values(sp.ras[[1]])))
              z_ij <- values(x.hist)[pres,]
              f_ij <- values(x.fut)[pres,]
              d_ij <- (f_ij - z_ij)
              distances <- apply(d_ij, 1, function(x) sqrt(t(x) %*% Rs.inv %*% x))
              D <- mean(distances, na.rm = T)/1.96
              ras <- raster(sp.ras[[1]])
              ras[pres] <- distances
            } else {
              x.mask.h <- mask(x.hist, sp.ras[[1]])
              x.mask.f <- mask(x.fut, sp.ras[[1]])
              d_ij <- (x.mask.f - x.mask.h)
              ras <- calc(d_ij, function(x) sqrt(t(x) %*% Rs.inv %*% x))
              D <- cellStats(ras, mean)/1.96
              distances <- na.omit(values(ras))
            }

            depart <- methods::new("departure", call = call, departure = D, departure_ras = ras, present = length(pres))
            return(depart)
          }
)

#' @rdname departure
setMethod("departure",
          signature(x.hist = "RasterBrick", x.fut = "missing", s.dat = "cnfa"),
          function(x.hist, s.dat, scale = FALSE){

            call <- match.call()
            sp.ras <- raster(s.dat)
            if(!identicalCRS(x.hist, sp.ras)) {stop("climate and species projections do not match")}
            if(length(raster::intersect(extent(x.hist), extent(sp.ras)))==0) {stop("climate and species data to not overlap")}
            if(scale) {
              center <- cellStats(x.hist, mean)
              sds <- cellStats(x.hist, sd)
              x.hist <- raster::scale(x.hist, center = center, scale = sds)
            }
            x.hist <- crop(x.hist, extent(sp.ras))
            Rs.inv <- solve(s.dat@s.cov, tol = 1e-20)
            #Rs.inv <- solve(covmat(sp.ras), tol = 1e-20)

            small <- canProcessInMemory(x.hist, 5)
            if(small){
              pres <- which(!is.na(values(sp.ras[[1]])))
              #d_ij <- values(x.hist)[pres,]
              # ref <- values(sp.ras)[pres,]
              # ref[,1] <- ref[,1] - (s.dat@marginality*1.96)^2
              # ref <- mean(apply(ref, 1, norm, "2"))
              d_ij <- values(x.hist)[pres,] #%*% as.matrix(s.dat@co)
              distances <- apply(d_ij, 1, function(x) sqrt(t(x) %*% Rs.inv %*% x))
              D <- mean(distances, na.rm = T)
              ras <- raster(sp.ras[[1]])
              ras[pres] <- distances
            } else {
              x.mask.h <- mask(x.hist, sp.ras[[1]])
              ras <- calc(x.mask.h, function(x) sqrt(t(x) %*% Rs.inv %*% x))
              D <- cellStats(ras, mean)
            }

            depart <- methods::new("departure", call = call, departure = D, departure_ras = ras, present = s.dat@present)
            return(depart)
          }
)

setMethod("show",
          signature = "departure",
          function(depart){
  if (!inherits(depart, "departure"))
    stop("Object of class 'departure' expected")
  cat("CLIMATIC DEPARTURE")
  cat("\ndeparture: ")
  cat(signif(depart@departure, 4))
  cat("\nnumber of cells present: ")
  cat(depart@present)
}
)

departure2 <- function(x.hist, s.dat, scale = FALSE, ...){

  call <- match.call()
  sp.ras <- raster(s.dat)
  if(!identicalCRS(x.hist, sp.ras)) {stop("climate and species projections do not match")}
  if(length(raster::intersect(extent(x.hist), extent(sp.ras)))==0) {stop("climate and species data to not overlap")}
  if(scale) {
    center <- cellStats(x.hist, mean)
    sds <- cellStats(x.hist, sd)
    x.hist <- raster::scale(x.hist, center = center, scale = sds)
  }
  x.hist <- crop(x.hist, extent(sp.ras))
  Rs.inv <- solve(covmat(sp.ras, ...), tol = 1e-20)

  small <- canProcessInMemory(x.hist, 5)
  if(small){
    pres <- which(!is.na(values(sp.ras[[1]])))
    d_ij <- values(x.hist)[pres,] %*% as.matrix(s.dat@co)
    distances <- apply(d_ij, 1, function(x) sqrt(t(x) %*% Rs.inv %*% x))
    D <- mean(distances, na.rm = T)
    ras <- raster(sp.ras[[1]])
    ras[pres] <- distances
  } else {
    x.mask.h <- mask(x.hist, sp.ras[[1]])
    cl <- makeCluster(getOption("cl.cores", cores))
    #clusterExport(cl, c("x.mask.h", "calc", "Rs.inv", "s.dat@co", "s"), envir = environment())
    registerDoSNOW(cl)

    f1 <- function(x) x %*% as.matrix(s.dat@co)
    x.mask.h2 <- clusterR(x.mask.h, fun = calc, args = list(fun = f1), export = "s.dat@co", m = cores)
    f2 <- function(x) sqrt(t(x) %*% Rs.inv %*% x)
    ras <- clusterR(x.mask.h2, fun = calc, args = list(fun = f2), export = "Rs.inv", m = cores)
    stopCluster(cl)
  }
  #x.mask.h <- calc(x.mask.h, function(x) x %*% as.matrix(s.dat@co))
  #ras <- calc(x.mask.h, function(x) sqrt(t(x) %*% Rs.inv %*% x))
  D <- cellStats(ras, mean)

  depart <- methods::new("departure", call = call, departure = D, departure_ras = ras, present = s.dat@present)
  return(depart)
}

departure2 <- function(x.hist, s.dat, scale = FALSE, cores = 1, filename = "", ...){

  out <- raster(x.hist)

  cl <- getCluster(cores)
  on.exit( returnCluster() )

  nodes <- length(cl)

  bs <- blockSize(x, minblocks=nodes*4)
  pb <- pbCreate(bs$n)

  # the function to be used (simple example)
  f1 <- function(x) x %*% as.matrix(s.dat@co)

  # get all nodes going
  for (i in 1:nodes) {
    sendCall(cl[[i]], f1, i, tag=i)
  }

  filename <- trim(filename)
  if (!canProcessInMemory(out) & filename == "") {
    filename <- rasterTmpFile()
  }

  if (filename != "") {
    out <- writeStart(out, filename=filename, ... )
  } else {
    vv <- matrix(ncol=nrow(out), nrow=ncol(out))
  }
  for (i in 1:bs$n) {
    # receive results from a node
    d <- recvOneData(cl)

    # error?
    if (! d$value$success) {
      stop('cluster error')
    }

    # which block is this?
    b <- d$value$tag
    cat('received block: ',b,'\n'); flush.console();

    if (filename != "") {
      out <- writeValues(out, d$value$value, bs$row[b])
    } else {
      cols <- bs$row[b]:(bs$row[b] + bs$nrows[b]-1)
      vv[,cols] <- matrix(d$value$value, nrow=out@ncols)
    }

    # need to send more data?
    ni <- nodes + i
    if (ni <= bs$n) {
      sendCall(cl[[d$node]], clFun, ni, tag=ni)
    }
    pbStep(pb)
  }
  if (filename != "") {
    out <- writeStop(out)
  } else {
    out <- setValues(out, as.vector(vv))
  }
  pbClose(pb)

  return(out)
}
