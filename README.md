# SUMMER
Spatio-temporal Under-five Mortality Models of Estimation in R

[![Build Status](https://travis-ci.org/bryandmartin/SUMMER.svg?branch=master)](https://travis-ci.org/bryandmartin/SUMMER)

## Installation - CRAN

The package is now available on CRAN. The easiest way to downlaod is to install directly using the code below.

``` r 
install.packages("SUMMER")
```

## Installation - Development Version

To download the development version of the SUMMER package, use the code below.

``` r
# install.packages("devtools")
devtools::install_github("martinbryan/SUMMER")
```

**Development Version Updates**

We are currently incorporating a new yearly model into the development version of SUMMER, improving the implementation of our the five yearly model, and adding new options to the functions. The development package may be unstable until we finalize this process.

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
