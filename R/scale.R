#' @keywords internal

.scale <- function(x, filename = filename, progress = progress, median = F, overwrite = T){
  if(!median){
    ifelse(canProcessInMemory(x), {
      br <- brick(x)
      if(progress) cat("Scaling data...")
      v <- values(x)
      t <- scale(v, center = T, scale = T)
      br <- setValues(br, t)
      center <- attr(t, "scaled:center")
      sds <- attr(t, "scaled:scale")
      y <- list(br, center, sds)
    }, {

      br <- brick(x)
      if(progress) cat("Calculating layer means...")
      center <- cellStats(x, 'mean', na.rm = T)
      if(progress) cat("Calculating layer sds...")
      sds <- cellStats(x, 'sd', na.rm = T)
      if(progress) cat("Scaling data...")
      t <- raster::scale(x, center = center, scale = sds)
      if(progress) cat("Writing data to file...")
      br <- writeRaster(t, filename = filename, overwrite = overwrite)
      y <- list(br, center, sds)
    }
    )
  }

  # if(median){
  #   ifelse(canProcessInMemory(x), {
  #     br <- brick(x)
  #     if(progress) cat("Scaling data...")
  #     v <- values(x)
  #     center <- apply(v, 2, median, na.rm = T)
  #     mads<-apply(v,2,mad,na.rm=T)
  #     t <- (v-center)/mads
  #     br <- setValues(br, t)
  #     y<-list(br,center,sds)
  #   }, {
  #
  #     br <- brick(x)
  #     if(progress) cat("Calculating layer medians...")
  #     center <- cellStats(x, 'median', na.rm=TRUE)
  #     if(progress) cat("Calculating layer MADs...")
  #     mads <- cellStats(x, 'mad', na.rm=TRUE)
  #     if(progress) cat("Scaling data...")
  #     t <- raster::scale(x, center=center, scale=mads)
  #     if(progress) cat("Writing data to file...")
  #     br <- writeRaster(t, filename = filename,overwrite=T)
  #     y<-list(br,center,mads)
  #   }
  #   )
  # }

  return(y)
}

.covmat <- function(x, cores = 1, .scale = FALSE, progress = TRUE, ...){
  stopifnot(is.numeric(cores) & cores >= 0)

  small <- canProcessInMemory(x)
  if(small){
    dat <- na.omit(values(x))
    mat <- cov(dat, method = "pearson")
    return(mat)
  }
  nl <- nlayers(x)
  mat <- matrix(NA, nrow = nl, ncol = nl)
  colnames(mat) <- rownames(mat) <- names(x)

  ii <- rep(1, nl)
  for(i in 2:nl) ii <- c(ii, rep(i, each = nl-i+1))
  jj <- 1:nl
  for(i in 2:nl) jj <- c(jj, i:nl)
  s <- 1:length(ii)

  if(.scale){
    means <- cellStats(x, stat = 'mean', na.rm = T)
    sds <- cellStats(x, stat = 'sd', na.rm = T)
    x <- (x - means)/sds
  }

  if(cores == 1){
    if(progress){
      pboptions(type = "txt", char = "=", txt.width = NA)
      result <- pbsapply(s, function(p) do.call(.covij, list(x = x, i = ii[p], j = jj[p])))
    } else {result <- sapply(s, function(p) do.call(covij, list(x = x, i = ii[p], j = jj[p])))}
  }

  if(cores > 1){
    cl <- makeCluster(getOption("cl.cores", cores))
    clusterExport(cl, c(".covij", "raster", "cellStats", "x", "ii", "jj", "s", "nl", "canProcessInMemory", "values"), envir = environment())
    registerDoSNOW(cl)
    if(progress){
      pb <- txtProgressBar(min = 0, max = length(s), style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      result <- foreach(p = s, .combine = c, .options.snow = opts) %dopar% {
        do.call(.covij, list(x = x, i = ii[p], j = jj[p]))
      }
      close(pb)
    }   else if(!progress){
      result <- foreach(p = s, .combine = c, .options.snow = opts) %dopar% {
        do.call(.covij, list(x = x, i = ii[p], j = jj[p]))
      }
    }
    stopCluster(cl)
  }

  for(p in s){
    mat[ii[p], jj[p]] <- mat[jj[p], ii[p]] <- result[p]
  }

  closeAllConnections()
  return(mat)
}
