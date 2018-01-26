
-   [CENFA: Climate and Ecological Niche Factor Analysis](#cenfa-climate-and-ecological-niche-factor-analysis)
    -   [Installation](#installation)
    -   [Examples](#examples)
        -   [`enfa`](#enfa)
        -   [`cnfa`](#cnfa)
    -   [Useful raster functions](#useful-raster-functions)
        -   [`parScale`](#parscale)
        -   [`parCov`](#parcov)

CENFA: Climate and Ecological Niche Factor Analysis
===================================================

`CENFA` provides tools for performing ecological-niche factor analysis (ENFA) and climate-niche factor analysis (CNFA). This package was created with three goals in mind:

-   To update the ENFA method for use with large datasets and modern data formats.
-   To correct a critical error in the ENFA method itself, that has managed to persist since Hirzel et al. first introduced ENFA in 2002.
-   To expand the application of ENFA in the context of climate change in order to quantify different aspects of species vulnerability to climate change, and to facilitate quantitative comparisons of vulnerability between species.

`CENFA` takes advantage of the `raster` and `sp` packages, allowing the user to conduct analyses directly with raster, shapefile, and point data, and to handle large datasets efficiently via partial data loading and parallelization.

Installation
------------

You can install CENFA from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("rinnan/CENFA")
```

Examples
--------

### `enfa`

We will use some example datasets to perform a basic ENFA. The historical climate dataset `climdat.hist` is a RasterBrick of 10 climate variables, covering much of the western US coast. `QUGA` is a SpatialPolygonsDataFrame of the historical range map of Garry oak (*Quercus garryana*).

A plot of the data, using the one of the layers of `climdat.hist`:

![](README-QUGA-plot-1.png)

The `enfa` function takes three basic arguments: the dataset of ecological variables (`climdat.hist`), the map of species presence (`QUGA`), and the values of `QUGA` that specify presence (in this case, a column named "CODE"). Calling the `enfa` object by name provides a standard summary of the ENFA results.

``` r
mod.enfa <- enfa(x = climdat.hist, s.dat = QUGA, field = "CODE")
mod.enfa
#> ENFA
#> 
#> Original function call: enfa(x = climdat.hist, s.dat = QUGA, field = "CODE")
#> 
#> Marginality factor: 
#>   MDR   ISO    TS HMmax CMmin   PWM   PDM    PS   PWQ   PDQ 
#> -0.13  0.51 -0.67 -0.03  0.67  0.71 -0.13  0.70  0.72  0.13 
#> 
#> Eigenvalues of specialization: 
#>  Marg Spec1 Spec2 Spec3 Spec4 Spec5 Spec6 Spec7 Spec8 Spec9 
#>  5.12  8.92  4.82  3.01  2.55  2.01  1.28  0.77  0.68  0.36 
#> 
#> Percentage of specialization contained in ENFA factors: 
#>  Marg Spec1 Spec2 Spec3 Spec4 Spec5 Spec6 Spec7 Spec8 Spec9 
#> 17.36 30.24 16.32 10.20  8.64  6.80  4.33  2.61  2.30  1.21 
#> 
#> Overall marginality:  1.191 
#> 
#> Overall specialization:  0.543 
#> 
#> Significant ENFA factors: 
#>        Marg Spec1 Spec2 Spec3
#> PWQ    0.61  0.03 -0.19 -0.25
#> PWM    0.60 -0.08  0.26  0.27
#> PS     0.59  0.12 -0.11 -0.02
#> TS    -0.56 -0.34 -0.54  0.39
#> CMmin  0.56 -0.54  0.00  0.49
#> ISO    0.43  0.11 -0.46 -0.05
#> MDR   -0.11 -0.46  0.59  0.38
#> PDM   -0.11  0.10  0.09  0.03
#> PDQ    0.11  0.00 -0.11 -0.04
#> HMmax -0.02  0.58 -0.07 -0.57
```

For larger datasets, we can speed up the computation via parallelization. We provide two additional arguments, `parallel = TRUE`, and `n`, which specifies the number of cores to use. `n` has a default value of 1, so only setting `parallel = TRUE` will not parallelize the function by itself.

``` r
# does not enable parallelization
mod <- enfa(x = climdat.hist, s.dat = QUGA, field = "CODE", parallel = TRUE)

# enables parallelization across 4 cores
mod <- enfa(x = climdat.hist, s.dat = QUGA, field = "CODE", parallel = TRUE, n = 4)
```

The function will attempt to match the value provided to `n` with the number of cores detected on the local device via `parallel::detectCores()`; if the provided `n` is greater than the number of available cores `k`, a warning will be issued and `n` will be set to `k - 1`.

### `cnfa`

The `cnfa` function is very similar to `enfa`, but performs a slightly different analysis. Whereas ENFA returns a *specialization factor* (the eigenvalues of specialization) describing the amount of specialization found in each *ENFA factor*, CNFA returns a *sensitivity factor* that reflects the amount of specialization found in each *ecological variable*. This makes the sensitivity factor more directly comparable to the marginality factor, and more interpretable in the context of species' sensitivity to a given variable.

``` r
mod.cnfa <- cnfa(x = climdat.hist, s.dat = QUGA, field = "CODE")
mod.cnfa
#> CNFA
#> 
#> Original function call: cnfa(x = climdat.hist, s.dat = QUGA, field = "CODE")
#> 
#> Marginality factor: 
#>   MDR   ISO    TS HMmax CMmin   PWM   PDM    PS   PWQ   PDQ 
#> -0.13  0.51 -0.67 -0.03  0.67  0.71 -0.13  0.70  0.72  0.13 
#> 
#> Sensitivity factor: 
#>   MDR   ISO    TS HMmax CMmin   PWM   PDM    PS   PWQ   PDQ 
#>  0.38  0.26  0.37  0.31  0.36  0.30  0.12  0.17  0.27  0.08 
#> 
#> Percentage of specialization contained in CNFA factors: 
#>  Marg Spec1 Spec2 Spec3 Spec4 Spec5 Spec6 Spec7 Spec8 Spec9 
#> 17.36 30.24 16.32 10.20  8.64  6.80  4.33  2.61  2.30  1.21 
#> 
#> Overall marginality:  1.191 
#> 
#> Overall sensitivity:  9.213 
#> 
#> Significant CNFA factors: 
#>        Marg Spec1 Spec2 Spec3
#> PWQ    0.61 -0.03  0.19  0.25
#> PWM    0.60  0.08 -0.26 -0.27
#> PS     0.59 -0.12  0.11  0.02
#> TS    -0.56  0.34  0.54 -0.39
#> CMmin  0.56  0.54  0.00 -0.49
#> ISO    0.43 -0.11  0.46  0.05
#> MDR   -0.11  0.46 -0.59 -0.38
#> PDM   -0.11 -0.10 -0.09 -0.03
#> PDQ    0.11  0.00  0.11  0.04
#> HMmax -0.02 -0.58  0.07  0.57
```

Using the `sensitivity_map` function, we can create a habitat map that identifies where we expect the species to be most sensitivite to changes in climate.

``` r
s.map <- sensitivity_map(mod.cnfa)
```

![](README-sensitivity-map-1.png)

Useful raster functions
-----------------------

The `raster` package contains the `clusterR` function, which enables parallelization methods for certain raster operations. `clusterR` only works on functions that operate on a cell-by-cell basis, however, which limits its usefulness. The `CENFA` package contains a few functions that speed up some basic `raster` functions considerably by parallelizing on a layer-by-layer basis rather than a cell-by-cell basis.

### `parScale`

The `parScale` function is identical to `raster::scale`, but has a parallelization option that will scale each raster layer in parallel. The `center` and `scale` arguments can be logical (`TRUE` or `FALSE`) or numeric vectors.

``` r
clim.scaled <- parScale(x = climdat.hist, parallel = TRUE, n = 4)
```

### `parCov`

The `parCov` function returnds the covariance matrix of a Raster\* object `x`, computing the covariance between each layer of `x`. This is similar to `raster::layerStats(x, stat = 'cov')`, but much faster when parallelization is employed.

``` r
mat <- parCov(x = climdat.hist, parallel = TRUE, n = 4)
```

Additionally, `parCov` can accept two Raster\* objects as arguments, similar to `stats::cov(x, y)`. If two Raster\* objects are supplied, then the covariance is calculated between the layers of `x` and the layers of `y`.

``` r
mat <- parCov(x = climdat.hist, y = climdat.fut, parallel = TRUE, n = 4)
```
