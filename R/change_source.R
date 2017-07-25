change_source<-function(x,filename){
  x@file@name <- filename
}

change_source<-function(x,filename){
  x@global_ras@file@name <- filename
}
