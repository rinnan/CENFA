#' CNFA raster extraction
#'
#' Extracts Raster* object from an object of class \code{cnfa}.
#'
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
          signature(x = "GLcnfa"),
          function(x){
            temp <- x@global_ras
            return(temp)
          }
)
