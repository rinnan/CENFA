#' GLcenfa
#'
#' This function is used to facilitate comparisons between species in the same study area. It speeds up the computation of multiple CNFAs or ENFAs by calculating the global covariance matrix as a first step, which can then be fed into the \code{\link{cnfa}} or \code{\link{enfa}} functions as their first argument. This saves the user from having to calculate the global covariance matrix for each species, which can take quite a bit of time.
#'
#' @aliases print.GLcenfa, show.GLcenfa
#' @param x Raster* object, typically a brick or stack of p environmental raster layers.
#' @param scale logical. If \code{TRUE} then the values of \code{x} will
#' be centered and scaled. Depending on the resolution of the data and the
#' extent of the study area, this can be quite time consuming.
#' @param progress logical. If \code{TRUE} then progress updates are printed.
#' @param cores numeric. Number of cores to utilize for speedier parallel computation of the covariance matrix.
#' @param filename character. Optional filename to save the RasterBrick output to file. If this is not provided, a temporary file will be created for large \code{x}.
#' @return Returns an S4 object of class \code{GLcenfa} with the following components:
#' @return global_ras Raster* brick \code{x} with p layers, possibly centered and scaled.
#' @return cov matrix. Global p x p covariance matrix.
#' @return center numeric. Layer means of \code{x}.
#' @return sd numeric. Layer standard deviations of \code{x}.
#' @return ncells numeric. Total number of raster cells with environmental data.
#' @seealso \code{\link{cnfa}}, \code{\link{enfa}}
#' @export
#'

setGeneric("GLcenfa", function(x, ...) {
  standardGeneric("GLcenfa")
})

#' @rdname GLcenfa
setMethod("GLcenfa",
          signature(x = "RasterBrick"),
          function(x, scale = FALSE, filename = '', progress = TRUE, cores = 1, ...){

            out <- brick(x)

            small <- canProcessInMemory(x, 4)
            filename <- trim(filename)

            if (!small & filename == ''){
              filename <- rasterTmpFile()
            }

            filename <- raster:::.fullFilename(filename, expand = T)

            if (!file.exists(dirname(filename))) {
              stop("Attempting to write a file to a path that does not exist:\n  ", dirname(filename))
            }

            if (scale){
              Z <- .scale(x, filename = filename, progress = progress, ...)
              names(Z[[1]]) <- names(x)
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

            gpres <- which(!is.na(values(max(x))))

            if (small){
              if(progress) cat("\nCalculating global covariance matrix...")
              Z <- values(x)[gpres, ]
              cov <- crossprod(Z, Z/length(gpres))
            }

            if (!small){
              if(progress) cat("\nCalculating global covariance matrix...\n")
              cov <- .covmat(x, cores = cores, .scale = F)
            }
            rownames(cov) <- colnames(cov) <- names(x)

            GLcenfa <- methods::new("GLcenfa", global_ras = x, cov = cov, center = center, sd = sds, ncells = length(gpres))
            return(GLcenfa)
          }
)

#' @rdname GLcenfa
setMethod("GLcenfa",
          signature(x = "RasterStack"),
          function(x, scale = FALSE, filename = '', progress = TRUE, cores = 1, ...){

            out <- brick(x)

            small <- canProcessInMemory(x, 4)
            filename <- trim(filename)

            if (!small & filename == ''){
              filename <- rasterTmpFile()
            }

            filename <- raster:::.fullFilename(filename, expand = T)

            if (!file.exists(dirname(filename))) {
              stop("Attempting to write a file to a path that does not exist:\n  ", dirname(filename))
            }

            if (scale){
              Z <- .scale(x, filename = filename, progress = progress, ...)
              names(Z[[1]]) <- names(x)
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

            gpres <- which(!is.na(values(max(x))))

            if (small){
              if(progress) cat("\nCalculating global covariance matrix...")
              Z <- values(x)[gpres, ]
              cov <- crossprod(Z, Z/length(gpres))
            }

            if (!small){
              if(progress) cat("\nCalculating global covariance matrix...\n")
              cov <- .covmat(x, cores = cores, .scale = F)
            }
            rownames(cov) <- colnames(cov) <- names(x)

            GLcenfa <- methods::new("GLcenfa", global_ras = x, cov = cov, center = center, sd = sds, ncells = length(gpres))
            return(GLcenfa)
          }
)
