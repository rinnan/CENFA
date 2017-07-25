#' cnfamad-class
#'
#' An object of class \code{cnfamad} is created from performing climate-niche factor analysis on species presence data using the \code{cnfa} function.
#'
#' @slot  call original function call
#' @slot nf number of kept specialization factors
#' @slot m marginalities
#' @slot s specialization factors
#' @slot map Raster map of habitat suitability
#' @export

setClass("cnfamad",slots=list(call = "call", mf = "numeric", marginality = "numeric", s = "numeric",
                           specialization = "numeric", co = "data.frame",
                           species_ras = "Raster", present = "numeric"))

setMethod ("show", "cnfamad", function(object){
  if (!inherits(object, "cnfamad"))
    stop("Object of class 'cnfamad' expected")
  cat("CNFA")
  cat("\nMarginality: ")
  cat(signif(object@marginality, 4))
  cat("\nSpecialization: ")
  cat(signif(object@specialization, 4))
  cat("\nEigenvalues of specialization: ")
  cat("\n")
  v <- signif(abs(object@s)/sum(abs(object@s)),2)
  v <- paste0(100*v,"%")
  s<-data.frame(t(object@s))
  names(s)<-v
  row.names(s)<-""
  l0 <- length(object@s)
  print(format(s[,1:(min(5, l0))]))
  if (l0 > 5)
    cat(" ...")
  cat("\n")
  cat("\n")
  co<-as.data.frame(object@co[order(abs(object@co$Marg),decreasing = T),])
  print(co)
}
)
