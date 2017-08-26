#' cnfa-class
#'
#' An object of class \code{cnfa} is created from performing climate-niche factor analysis on species presence data using the \code{cnfa} function.
#'
#' @slot call Original function call.
#' @slot nf Number of kept specialization factors.
#' @slot mf marginality factor. Vector that describes the location of the species Hutchinsonian niche relative to the global niche.
#' @slot marginality Standardized magnitude of the marginality factor.
#' @slot sf specialization factors. Vector of eigenvalues.
#' @slot specialization The square root of the sum of \code{sf} divided by the length of \code{sf}.
#' @slot s.prop Vector representing the proportion of specialization found in each factor.
#' @slot co A data frame of variable loadings with p rows and nf + 1 columns.
#' @slot ras Raster* object of transformed climate values, with number of layers equal to nf + 1.
#' @slot present numeric. Number of raster cells in which species is present.
#' @export

setClass("cnfa", slots = list(call = "call", nf = "numeric", mf = "numeric", marginality = "numeric", sf = "numeric",
                           specialization = "numeric", s.prop = "numeric", co = "data.frame",
                           ras = "Raster", s.cov = "matrix", present = "numeric"))

setMethod ("show", "cnfa", function(object){
  if (!inherits(object, "cnfa"))
    stop("Object of class 'cnfa' expected")
  cat("CNFA")
  cat("\nOriginal function call: ")
  cat(object@call)
  cat("\nMarginality: ")
  cat(signif(object@marginality, 4))
  cat("\nSpecialization: ")
  cat(signif(object@specialization, 4))
  cat("\nEigenvalues of specialization: ")
  cat("\n")
  v <- paste0(100*object@s.prop, "%")
  s <- data.frame(t(object@sf))
  names(s) <- v
  row.names(s) <- ""
  l0 <- length(object@sf)
  print(format(s[ , 1:(min(5, l0))]))
  if (l0 > 5)
    cat(" ...")
  cat("\n")
  cat("\n")
  co <- as.data.frame(object@co[order(abs(object@co$Marg), decreasing = T), ])
  print(co)
}
)
