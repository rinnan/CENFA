.covij <- function(x, i, j, n = ncell(x), sample = T){
  sm <- canProcessInMemory(raster(x, layer = i))
  if(sm){
    r <- values(raster(x, layer=i) * raster(x, layer=j))
    nn <- length(r[!is.na(r)])
    ifelse(sample == T,
           v <- sum(r, na.rm = T)/(nn - 1),
           v <- sum(r, na.rm = T)/nn)
  }
  if(!sm){
    r <- raster(x, layer=i) * raster(x, layer=j)
    ifelse(sample == T,
      v <- cellStats(r, stat='sum', na.rm=T) / (n - cellStats(r, stat='countNA') - 1),
      v <- cellStats(r, stat='sum', na.rm=T) / (n - cellStats(r, stat='countNA')))
    f <- filename(r)
    file.remove(c(f, extension(f, '.gri')))
  }
  return(v)
}
