
# News

## CENFA 1.1.2

Published on CRAN, 2023-07-18.

- Updated methods for detecting projection mismatches between climate
  and species data. (`sp::identicalCRS()` no longer works reliably with
  Raster\* objects.)

## CENFA 1.1.1

Published on CRAN, 2021-06-08.

## CENFA 1.1.0.9000

In current development.

- Fixed bug in `enfa()` function that returned error for large datasets.

- Fixed bug in `predict()` methods for `enfa` objects that returned
  raster brick instead of layer.

- Removed `closeAllConnections()` command that may cause unwanted
  behavior.

## CENFA 1.1.0

Published on CRAN, 2020-02-15.

## CENFA 1.0.0.9000

- Updated descriptions of CNFA sensitivity and ENFA specialization.

- Fixed small bug in `parCov()` function that caused NA observations to
  be discarded.

- Added more informative error messages for `cnfa()` function.

- Forward compatibility with upcoming R 4.0 release.

## CENFA 1.0.0

Published on CRAN, 2018-11-06.

## CENFA 0.1.0.9000

- Changed `quiet = TRUE` arguments to `progress = FALSE` for more
  intuitive interface.

- Different definition of CNFA sensitivity factor and overall
  sensitivity, to better agree with ENFA specialization factor.

- Different definition of vulnerability factor and overall
  vulnerability.

- Fixed bugs in `parScale` function involving parallel methods.

- Added predict methods for `enfa`, `cnfa`, `departure` and
  `vulnerability` objects.

- Fixed bug in `GLcenfa` that prevented writing of Raster\* objects.

- Imports `raster` package instead of depends.

- Fixes bug in `parCov` function involving covariance calculations
  between two Raster\* objects.

## CENFA 0.1.0

Published on CRAN, as of 2018-02-06!
