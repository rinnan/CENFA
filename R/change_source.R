#change_source<-function(x,filename){
#  x@file@name <- filename
#}

change_source <- function(x, filename){
  x@ras@file@name <- filename
  x
}
