#' cnfa-class
#'
#' An object of class \code{cnfa} is created by performing climate-niche factor
#' analysis on species presence data using the \code{cnfa} function.
#'
#' @slot call Original function call
#' @slot mf numeric. Named vector representing the marginality factor, describing
#'   the location of the species niche relative to the global niche
#' @slot marginality numeric. Magnitude of the marginality factor \code{mf}, scaled
#'   by the global covariance matrix
#' @slot sf numeric. Named vector representing the sensitivity factor
#' @slot sensitivity numeric. The magnitude of the sensitivity factor \code{sf},
#'   scaled by the global covariance matrix
#' @slot p.spec numeric. Named vector representing the proportion of specialization
#'   found on each factor
#' @slot co p x p matrix of standardized variable loadings
#' @slot cov p x p species covariance matrix
#' @slot ras RasterBrick of transformed climate values, with p layers
#' @slot weights Raster layer of weights used for CNFA calculation
#'
#' @export

setClass("cnfa", slots = list(call = "call", mf = "numeric", marginality = "numeric", sf = "numeric",
                           sensitivity = "numeric", p.spec = "numeric", co = "matrix", cov = "matrix", ras = "Raster", weights = "Raster"))

setMethod ("show", "cnfa", function(object){
  if (!inherits(object, "cnfa"))
    stop("Object of class 'cnfa' expected")
  cat("CNFA\n")
  cat("\nOriginal function call: ")
  print(object@call)
  cat("\nMarginality factor: \n")
  j <- order(abs(object@mf), decreasing = T)
  print(round(object@mf[j], 2))
  cat("\nSensitivity factor: \n")
  print(sort(round(object@sf, 2), decreasing = T))
  cat("\nPercentage of specialization contained in factors: \n")
  print(round(100*object@p.spec, 2))
  cat("\nOverall marginality: ", round(object@marginality, 3), "\n")
  cat("\nOverall sensitivity: ", round(object@sensitivity, 3), "\n")
  cat("\nSignificant CNFA factors: \n")
  n <- brStick(object@p.spec[-1])
  co <- as.data.frame(object@co[order(abs(object@co[,1]), decreasing = T), ])
  print(round(co[, 1:(n+1)], 2))
}
)
