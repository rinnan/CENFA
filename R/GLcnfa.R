#' Climate Niche Factor Analysis
#'
#' Performs climate niche factor analysis using climate raster data and species presence data in shapefile.
#'
#' @aliases print.cnfa, show.cnfa
#' @param climdat Raster* object, typically a brick or stack of climate raster layers
#' @param speciesdat SpatialPolygons* object detailing species presence or abundance
#' @param field field of \code{speciesdat} that specifies presence or abundance. This is equivalent to the \code{field} argument in the \code{raster} package.
#' @param nf integer. Specifies the number of specialization axes to keep after transformation.
#' @param scale logical. If \code{TRUE} then the values of the Raster* object will
#' be centered and scaled. Depending on the resolution of the climate data and the
#' extent of the study area, this can be quite time consuming. If running this
#' function for multiple species, it is recommended that the climate data be scaled beforehand.
#' @param sp.prj character. Spatial projection of species data.
#' @return Returns an S4 object of class \code{cnfa} with the following components:
#' @return call original function call
#' @return tab a dataframe
#' @export
#'
#'

setGeneric("GLcnfa", function(x, ...) {
  standardGeneric("GLcnfa")
})

#' @rdname GLcnfa
setMethod("GLcnfa",
          signature(x = "RasterBrick"),
          function(x, scale = FALSE, filename = '', progress = TRUE, cores = 1){

            out <- brick(x)

            small <- canProcessInMemory(x, 4)
            filename <- trim(filename)

            if (!small & filename == ''){
              filename <- rasterTmpFile()
            }

            filename <- raster:::.fullFilename(filename, expand = TRUE)

            if (!file.exists(dirname(filename))) {
              stop("Attempting to write a file to a path that does not exist:\n  ", dirname(filename))
            }

            if (scale){
              Z <- .scale(x, filename = filename, progress = progress)
              x <- Z[[1]]
              center <- Z[[2]]
              sds <- Z[[3]]
            }

            if (!scale){
              cat("Warning: new raster will not be written to file.")
              if(progress) cat("\nCalculating layer means...")
              center <- cellStats(x, 'mean')
              if(progress) cat("\nCalculating layer sds...")
              sds <- cellStats(x, 'sd')
            }

            gpres <- which(!is.na(values(x[[1]])))

            if (small){
              if(progress) cat("\nCalculating covariance matrix...")
              Z <- values(x)[gpres, ]
              cov <- crossprod(Z, Z/length(gpres))
            }

            if (!small){
              if(progress) cat("\nCalculating covariance matrix...\n")
              cov <- .covmat(x, cores = cores, .scale = F)
            }
            rownames(cov) <- colnames(cov) <- names(x)

            GLcnfa <- methods::new("GLcnfa", global_ras = x, cov = cov, center = center, sd = sds, ncells = length(gpres), scale = scale)
            return(GLcnfa)
          }
)
