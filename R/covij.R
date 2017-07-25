.covij<-function(x,i,j,n=ncell(x)){
  sm <- canProcessInMemory(raster(x,layer=i))
  if(sm){
    r <- values(raster(x, layer=i) * raster(x, layer=j))
    nn <- length(r[!is.na(r)])
    v <- sum(r,na.rm = T)/(nn-1)
  }
  if(!sm){
    r <- raster(x, layer=i) * raster(x, layer=j)
    v <- cellStats(r, stat='sum', na.rm=T) / (n - cellStats(r, stat='countNA') - 1)
    f <- filename(r)
    file.remove(c(f, extension(f, '.gri')))
  }
  return(v)
}
