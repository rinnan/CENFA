#' enfa-class
#'
#' An object of class \code{enfa} is created from performing ecological-niche factor analysis on species presence data using the \code{enfa} function.
#'
#' @slot mf marginality factor. Vector that describes the location of the species Hutchinsonian niche.
#' @slot marginality. Standardized magnitude of the marginality factor.
#' @slot s specialization factors. Matrix of specializations.
#' @slot specialization specialization. Mahalanobis distance of something.
#' @slot ras Raster* object of transformed climate values, with number of layers equal to nf + 1.
#' @slot present Number of cells in which species is present.
#' @export

setClass("enfa", slots=list(call = "call", mf = "numeric", marginality = "numeric", s = "numeric",
                           specialization = "numeric", spec.account = "numeric", co = "data.frame", ras = "Raster", present = "numeric"))

setMethod ("show", "enfa", function(object){
  if (!inherits(object, "enfa"))
    stop("Object of class 'enfa' expected")
  cat("ENFA")
  cat("\nmarginality: ")
  cat(signif(object@marginality, 4))
  cat("\nspeciality: ")
  cat(signif(object@specialization, 4))
  cat("\neigen values of specialization: ")
  l0 <- length(object@s)
  cat(signif(object@s, 4)[1:(min(5, l0))])
  if (l0 > 5)
    cat(" ...")
  cat("\n")
  cat("\n")
  sumry <- array("", c(5, 4), list(1:5, c("vector", "length",
                                          "mode", "content")))
  sumry[1, ] <- c("@mf", length(object@mf),
                  mode(object@mf), "coordinates of the marginality vector")
  sumry[2, ] <- c("@s", length(object@s),
                  mode(object@s), "eigen values of specialization")
  class(sumry) <- "table"
  print(sumry)
  cat("\n")
  sumry <- array("", c(3, 4), list(1:3, c("data.frame", "nrow",
                                          "ncol", "content")))
  sumry[1, ] <- c("@co", nrow(object@co), ncol(object@co), "column coordinates")
  class(sumry) <- "table"
  print(sumry)
  if (length(names(object)) > 11) {
    cat("\nother elements: ")
    cat(names(object)[12:(length(object))], "\n")
  }
}
)
