#' enfa-class
#'
#' An object of class \code{enfa} is created from performing ecological-niche
#' factor analysis on species presence data using the \code{enfa} function.
#'
#' @slot call Original function call
#' @slot mf numeric. Named vector representing the marginality factor, describing
#'   the location of the species niche relative to the global niche
#' @slot marginality numeric. Magnitude of the marginality factor \code{mf}
#' @slot sf numeric. Named vector representing the specialization factor,
#'   equivalent to the eigenvalues of specialization
#' @slot specialization numeric. The square of the sum of eigenvalues, divided
#'   by the length of \code{sf}
#' @slot p.spec numeric. Named vector representing the proportion of
#'   specialization found on each factor
#' @slot co p x p matrix of standardized variable loadings
#' @slot cov p x p species covariance matrix
#' @slot ras RasterBrick of transformed climate values, with p layers
#' @slot weights Raster layer of weights used for ENFA calculation
#' @export

setClass("enfa", slots = list(call = "call", mf = "numeric", marginality = "numeric",
                              sf = "numeric", specialization = "numeric", p.spec = "numeric",
                              co = "matrix", cov = "matrix", ras = "Raster", weights = "Raster"))

setMethod ("show", "enfa", function(object){
  if (!inherits(object, "enfa"))
    stop("Object of class 'enfa' expected")
  cat("ENFA\n")
  cat("\nOriginal function call: ")
  print(object@call)
  cat("\nMarginality factor: \n")
  j <- order(abs(object@mf), decreasing = T)
  print(round(object@mf[j], 2))
  cat("\nEigenvalues of specialization: \n")
  print(sort(round(object@sf, 2), decreasing = T))
  cat("\nPercentage of specialization contained in factors: \n")
  print(round(100*object@p.spec, 2))
  cat("\nOverall marginality: ", round(object@marginality, 3), "\n")
  cat("\nOverall specialization: ", round(object@specialization, 3), "\n")
  cat("\nSignificant ENFA factors: \n")
  n <- brStick(object@p.spec[-1])
  co <- as.data.frame(object@co[order(abs(object@co[,1]), decreasing = T), ])
  print(round(co[, 1:(n+1)], 2))
}
)
