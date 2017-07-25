#' Climate-Niche Factor Analysis
#'
#' Performs climate-niche factor analysis using climate raster data and species presence data in shapefile.
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

setGeneric("GLcnfamad", function(climdat,...) {
  standardGeneric("GLcnfamad")
})

#' @rdname GLcnfamad
setMethod("GLcnfamad",
          signature(climdat = "RasterBrick"),
          function(climdat, scale = FALSE, filename='',progress=TRUE,cores=1){

            out<-brick(climdat)
            #             Z <- getValues(climdat)
            #             rZ <- nrow(Z)
            #             cZ <- ncol(Z)
            #             center <- .colMeans(Z, rZ, cZ, na.rm=T)
            #             sds <- apply(Z, 2, sd, na.rm=T)
            small <- canProcessInMemory(climdat,3)
            filename <- trim(filename)



            if (!small & filename == ''){
              filename <- rasterTmpFile()
            }

            filename <- raster:::.fullFilename(filename, expand=TRUE)

            if (!file.exists(dirname(filename))) {
              stop("Attempting to write a file to a path that does not exist:\n  ", dirname(filename))
            }

            if (scale){
              Z <- .scale(climdat, filename = filename,progress=progress)
              climdat <- Z[[1]]
              center <- Z[[2]]
              mads <- Z[[3]]
            }

            if (!scale){
              #out <- climdat
              cat("Warning: new raster will not be written to file.")
              if(progress) cat("Calculating layer medians...")
              center <- cellStats(x, 'median', na.rm=TRUE)
              if(progress) cat("Calculating layer MADs...")
              mads <- cellStats(x, 'mad', na.rm=TRUE)
            }

            gpres <- which(!is.na(values(climdat[[1]])))

            if (small){
              if(progress) cat("\nCalculating covariance matrix...")
              Z <- values(climdat)[gpres,]
              cov <- crossprod(Z, Z/length(gpres))
            }

            if (!small){
              if(progress) cat("\nCalculating covariance matrix...\n")
              cov<-.covmat(climdat,cores=cores,.scale=F)
            }
            rownames(cov) <- colnames(cov) <- names(climdat)

            gl.cnfa <- methods::new("GLcnfamad", global_ras = climdat, cov = cov, center = center, mad = mads, ncells = length(gpres), scale = scale)
            return(gl.cnfa)
          }
)
