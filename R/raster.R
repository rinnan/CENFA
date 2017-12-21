#' CENFA raster extraction
#'
#' Extracts Raster* objects from objects of classes \code{\link{cnfa}}, \code{\link{enfa}}, and \code{\link{GLcenfa}}.
#'
#' @param x object of class \code{cnfa}, \code{enfa}, or \code{GLcenfa}
#' @return Raster* object.
#'
#' @export

#' @rdname raster
setMethod("raster",
          signature(x = "enfa"),
          function(x){
            temp <- x@ras
            return(temp)
          }
)

#' @rdname raster
setMethod("raster",
          signature(x = "cnfa"),
          function(x){
            temp <- x@ras
            return(temp)
          }
)

#' @rdname raster
setMethod("raster",
          signature(x = "GLcenfa"),
          function(x){
            temp <- x@global_ras
            return(temp)
          }
)
