#' Climate Niche Factor Analysis
#'
#' Functions for extracting data from slots of objects of classes \code{cnfa} and \code{enfa}.
#'
#' @param x cnfa or enfa object
#'
#' @examples
#' mod1 <- cnfa(x = climdat.hist, s.dat = ABPR, field = "CODE")
#' mf(mod1)
#'
#' @export

mf <- function(x){
  return(x@mf)
}

sf <- function(x){
  return(x@sf)
}

marginality <- function(x){
  if (!inherits(x, c("cnfa", "enfa"))) stop("Object of class 'cnfa' or 'enfa' expected")
  return(x@marginality)
}

specialization <- function(x){
  if (!inherits(x, "enfa"))
    stop("Object of class 'enfa' expected")
  return(x@specialization)
}

sensitivity <- function(x){
  if (!inherits(x, "cnfa"))
    stop("Object of class 'cnfa' expected")
  return(x@sensitivity)
}
