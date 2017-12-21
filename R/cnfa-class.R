#' cnfa-class
#'
#' An object of class \code{cnfa} is created from performing climate-niche factor analysis on species presence data using the \code{cnfa} function.
#'
#' @slot call Original function call.
#' @slot mf numeric. Named vector representing the marginality factor, describing the location of the species niche relative to the global niche.
#' @slot marginality numeric. Magnitude of the marginality factor \code{mf}.
#' @slot sf numeric. Named vector representing the sensitivity factor.
#' @slot sensitivity numeric. The magnitude of the sensitivity factor \code{sf}.
#' @slot p.spec numeric. Named vector representing the proportion of specialization found on each factor.
#' @slot co p x p matrix of standardized variable loadings.
#' @slot cov p x p species covariance matrix.
#' @slot present numeric. Number of raster cells in which species is present.
#' @slot ras RasterBrick of transformed climate values, with p layers.
#' @export

setClass("cnfa", slots = list(call = "call", mf = "numeric", marginality = "numeric", sf = "numeric",
                           sensitivity = "numeric", p.spec = "numeric", co = "matrix", cov = "matrix", present = "numeric", ras = "Raster"))

setMethod ("show", "cnfa", function(object){
  if (!inherits(object, "cnfa"))
    stop("Object of class 'cnfa' expected")
  cat("CNFA")
  cat("\nOriginal function call: ")
  cat(object@call)
  cat("\nMarginality factor: ")
  cat(signif(object@mf, 3))
  cat("\nSensitivity factor: ")
  cat(signif(object@sf, 3))
  cat("\nEigenvalues of specialization: ")
  cat("\n")
  v <- paste0(100*object@p.spec, "%")
  #s <- data.frame(t(object@sf))
  #names(s) <- v
  #row.names(s) <- ""
  #l0 <- length(object@sf)
  #print(format(s[ , 1:(min(5, l0))]))
  #if (l0 > 5)
  #  cat(" ...")
  cat("\n")
  cat("\n")
  co <- as.data.frame(object@co[order(abs(object@co$Marg), decreasing = T), ])
  print(co)
}
)
