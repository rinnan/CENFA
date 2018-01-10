#' @keywords internal

.covij <- function(x, y, w, sample = sample){
  sm <- canProcessInMemory(x)
  if(sm){
    x <- values(x)
    x <- x - mean(x, na.rm = T)
    y <- values(y)
    y <- y - mean(y, na.rm = T)
    if(!is.null(w)){
      w <- values(w)
      sumw <- sum(w, na.rm = T)
      #w <- w / sumw
      nn <- sumw
    r <- na.omit(w * x * y)
    } else {
      r <- na.omit(x * y)
      nn <- length(r)
    }
 #   ifelse(sample == T,
           v <- sum(r, na.rm = T)/(nn - sample)
#           v <- sum(r, na.rm = T)/nn)
  }
  if(!sm){
    x <- scale(x, scale = F)
    y <- scale(y, scale = F)
    if(!is.null(w)){
      sumw <- cellStats(w, sum, na.rm = T)
      #w <- w / sumw
      nn <- sumw
      #w <- w / cellStats(w, sum, na.rm = T)
      r <- w * x * y
    } else {
      r <- x * y
      nn <- length(r[!is.na(r)])
    }
    #ifelse(sample == T,
           v <- cellStats(r, stat='sum', na.rm = T) / (nn - sample)
           #v <- cellStats(r, stat='sum', na.rm=T) / nn)
    f <- filename(r)
    file.remove(c(f, extension(f, '.gri')))
  }
  return(v)
}
