.scale<-function(x, filename=filename, progress=progress,median=F){
  if(!median){
  ifelse(canProcessInMemory(x), {
    br <- brick(x)
    if(progress) cat("Scaling data...")
    v <- values(x)
    t <- scale(v, center=TRUE, scale=TRUE)
    br <- setValues(br, t)
    center<-attr(t, "scaled:center")
    sds<-attr(t,"scaled:scale")
    y<-list(br,center,sds)
  }, {

    br <- brick(x)
    if(progress) cat("Calculating layer means...")
    center <- cellStats(x, 'mean', na.rm=TRUE)
    if(progress) cat("Calculating layer sds...")
    sds <- cellStats(x, 'sd', na.rm=TRUE)
    if(progress) cat("Scaling data...")
    t <- raster::scale(x, center=center, scale=sds)
    if(progress) cat("Writing data to file...")
    br <- writeRaster(t, filename = filename,overwrite=T)
    y<-list(br,center,sds)
  }
  )
  }

  if(median){
    ifelse(canProcessInMemory(x), {
      br <- brick(x)
      if(progress) cat("Scaling data...")
      v <- values(x)
      center<-apply(v,2,median,na.rm=T)
      mads<-apply(v,2,mad,na.rm=T)
      t <- (v-center)/mads
      br <- setValues(br, t)
      y<-list(br,center,sds)
    }, {

      br <- brick(x)
      if(progress) cat("Calculating layer medians...")
      center <- cellStats(x, 'median', na.rm=TRUE)
      if(progress) cat("Calculating layer MADs...")
      mads <- cellStats(x, 'mad', na.rm=TRUE)
      if(progress) cat("Scaling data...")
      t <- raster::scale(x, center=center, scale=mads)
      if(progress) cat("Writing data to file...")
      br <- writeRaster(t, filename = filename,overwrite=T)
      y<-list(br,center,mads)
    }
    )
  }

  return(y)
}
