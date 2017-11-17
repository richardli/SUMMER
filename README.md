# SUMMER
Spatio-temporal Under-five Mortality Models of Estimation in R

[![Build Status](https://travis-ci.org/martinbryan/SUMMER.svg?branch=master)](https://travis-ci.org/martinbryan/SUMMER)

## Installation

Currently, the package is only available through Github. CRAN submission is waiting on the completion and publication of a manuscript. For now, use the code below to install SUMMER.

``` r
# install.packages("devtools")
devtools::install_github("martinbryan/SUMMER")
```

## Usage
To see example usage of all main functions, build the package vignette. Note that the vignette will take a few minutes to compile.

``` r
# install.packages("devtools")
devtools::install_github("martinbryan/SUMMER", build_vignettes = TRUE)
# Use this to view the vignette in the SUMMER HTML help
help(package = "SUMMER", help_type = "html")
# Use this to view the vignette as an isolated HTML file
utils::browseVignettes(package = "SUMMER")
```

## Bug Reports / Change Requests
If you encounter a bug or would like make a change request, please file it as an issue [here](https://github.com/martinbryan/SUMMER/issues).
