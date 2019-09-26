# SUMMER
Spatio-temporal Under-five Mortality Models of Estimation in R

[![Build Status](https://travis-ci.org/bryandmartin/SUMMER.svg?branch=master)](https://travis-ci.org/bryandmartin/SUMMER) [![](https://cranlogs.r-pkg.org/badges/SUMMER)](https://cran.r-project.org/package=SUMMER) [![](https://cranlogs.r-pkg.org/badges/grand-total/SUMMER?color=orange)](https://cran.r-project.org/package=SUMMER)

## Major update (version 0.3.0)
Version 0.3.0 contains some major updates from the previous versions. Some of the substantial change to existing functions are listed here. For a complete log of changes, see the [News](https://github.com/bryandmartin/SUMMER/blob/master/NEWS.md) section.

#### Function name changes
The following functions have been renamed. Most of the function arguments remain the same:

+ ``countrySummary`` is now ``getDirect``
+ ``cuontrySummary_mult`` is now ``getDirectList``
+ ``fitspace`` is now ``fitGeneric``
+ ``projINLA`` is now ``getSmooth``

#### New functions
The following new functions are added:

+ ``getDiag``: produce diagnostic plots for the fitted model.
+ ``getAdjusted``: produce adjusted estimates for a fitted model
+ ``getAmat``: automatic extract spatial adjacency matrix from the polygon file.
+ ``hatchPlot``: plot variables on a map with hatching indicating the width of the credible interval.

#### New methods
``fitINLA2``: implements new smoothing methods based on binomial models at cluster level. 



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




### Random notes for package authors:
+ In order to use `devtools::check` to check the package with static vignettes, use the hidden option `devtools::check(clean_doc = FALSE)` to avoid deleting the inst/doc folder.
+ Also for large PDF vignette, run `tools::compactPDF(gs_quality = "ebook", paths = "inst/doc/")` to please CRAN.

