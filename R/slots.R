#' Accessing CENFA slots
#'
#' Functions for extracting data from slots of objects of classes \code{cnfa} and \code{enfa}.
#'
#' @param x cnfa or enfa object
#'
#' @examples
#' mod1 <- cnfa(x = climdat.hist, s.dat = ABPR, field = "CODE")
#' m.factor(mod1)
#'
#' @name slot-access
NULL

#' @rdname slot-access
#' @export
m.factor <- function(x){
  return(x@mf)
}

#' @rdname slot-access
#' @export
s.factor <- function(x){
  return(x@sf)
}

#' @rdname slot-access
#' @export
marginality <- function(x){
  if (!inherits(x, c("cnfa", "enfa"))) stop("Object of class 'cnfa' or 'enfa' expected")
  return(x@marginality)
}

#' @rdname slot-access
#' @export
specialization <- function(x){
  if (!inherits(x, "enfa"))
    stop("Object of class 'enfa' expected")
  return(x@specialization)
}

#' @rdname slot-access
#' @export
sensitivity <- function(x){
  if (!inherits(x, "cnfa"))
    stop("Object of class 'cnfa' expected")
  return(x@sensitivity)
}
