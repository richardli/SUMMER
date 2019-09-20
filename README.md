# SUMMER
Spatio-temporal Under-five Mortality Models of Estimation in R

[![Build Status](https://travis-ci.org/bryandmartin/SUMMER.svg?branch=master)](https://travis-ci.org/bryandmartin/SUMMER) [![](https://cranlogs.r-pkg.org/badges/SUMMER)](https://cran.r-project.org/package=SUMMER) [![](https://cranlogs.r-pkg.org/badges/grand-total/SUMMER?color=orange)](https://cran.r-project.org/package=SUMMER)

## Major update (version 0.3.0)
Version 0.3.0 contains some major updates from the previous versions. Some of the substantial change to existing functions are listed here. For a complete log of changes, see the [News](https://github.com/bryandmartin/SUMMER/NEWS.md) section.

#### Function name changes
The following functions have been renamed. Most of the function arguments remain the same:

+ ``countrySummary`` is now ``getDirect``
+ ``cuontrySummary_mult`` is now ``getDirectList``
+ ``fitspace`` is now ``fitGeneric``
+ ``projINLA`` is now ``getSmooth``

For now, the old function names can still be called and they are internally linked to the new functions. But we encourage users to switch to the new functions, as they will be more actively maintained. In the next major release, the old functions will be removed.

#### New functions
The following new functions are added:

+ ``getDiag``: produce diagnostic plots for the fitted model.
+ ``getAdjusted``: produce adjusted estimates for a fitted model

#### New methods
``fitINLA`` function now implements new smoothing methods based on binomial models at cluster level. 



## Installation - CRAN

The package is now available on CRAN. The easiest way to download is to install directly using the code below.

``` r 
install.packages("SUMMER")
```

## Installation - Development Version

To download the development version of the SUMMER package, use the code below.

``` r
# install.packages("devtools")
devtools::install_github("bryandmartin/SUMMER")
```
 
## Usage
To see example usage of all main functions, build the package vignette. Note that the vignette will take a few minutes to compile.

``` r
# install.packages("devtools")
devtools::install_github("bryandmartin/SUMMER", build_vignettes = TRUE)
# Use this to view the vignette in the SUMMER HTML help
help(package = "SUMMER", help_type = "html")
# Use this to view the vignette as an isolated HTML file
utils::browseVignettes(package = "SUMMER")
```

## Bug Reports / Change Requests
If you encounter a bug or would like make a change request, please file it as an issue [here](https://github.com/bryandmartin/SUMMER/issues).
