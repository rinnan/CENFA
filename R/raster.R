#' CNFA raster extraction
#'
#' Extracts Raster* object from an object of class \code{cnfa}.
#'
#' @aliases
#' @param cnfa object of class \code{cnfa}
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
