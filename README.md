# SUMMER <img src="man/figures/SUMMER.png" align="right" width="120" />

[![](https://github.com/richardli/SUMMER/workflows/R-CMD-check_INLA_stable/badge.svg)](https://github.com/richardli/SUMMER/actions) [![](https://github.com/richardli/SUMMER/workflows/R-CMD-check_INLA_testing/badge.svg)](https://github.com/richardli/SUMMER/actions) 
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/SUMMER)](https://cran.r-project.org/package=SUMMER) [![](https://cranlogs.r-pkg.org/badges/SUMMER)](https://cran.r-project.org/package=SUMMER) [![](https://cranlogs.r-pkg.org/badges/grand-total/SUMMER?color=orange)](https://cran.r-project.org/package=SUMMER)

SAE Unit/area Models and Methods for Estimation in R

### Overview

**SUMMER** is an R package providing an extensive collection of space-time smoothing and small area estimation methods for prevalence estimation using complex survey data, with a special focus on Demographic Health Surveys (DHS) data, and the estimation of child mortality using full birth history data. The package also provides a collection of plotting functions to visualize the estimates in space and time.



### Citation

To cite the SUMMER package in publications use
```
  @Manual{li2020space,
    title = {Space-Time Smoothing of Demographic and Health Indicators using the R Package SUMMER},
    author = {Zehang R Li and Bryan D Martin and Tracy Q Dong and Geir-Arne Fuglstad and Jessica Godwin and John Paige and Andrea Riebler and Samuel Clark and Jon Wakefield},
    year = {2020},
    journal = {arXiv preprint}
  }
```

To cite specific version of the SUMMER package use
+ v1.0.0 
```
  @Manual{summer2020,
    title = {SUMMER: Spatio-Temporal Under-Five Mortality Methods for Estimation},
    author = {Zehang R Li and Bryan D Martin and Yuan Hsiao and Jessica Godwin and Jon Wakefield and Samuel J Clark and Geir-Arne Fuglstad and Andrea Riebler},
    year = {2020},
    note = {R package version 1.0.0},
  }
```
+ earlier versions (e.g., v0.3.0)
```
  @Manual{summer2019,
    title = {SUMMER: Spatio-Temporal Under-Five Mortality Methods for Estimation},
    author = {Bryan D Martin and Zehang R Li and Yuan Hsiao and Jessica Godwin and Jon Wakefield and Samuel J Clark and Geir-Arne Fuglstad and Andrea Riebler},
    year = {2019},
    note = {R package version 0.3.0},
  }
```


### Installation - CRAN

The package is now available on CRAN. The easiest way to download is to install directly using the code below.

``` r 
install.packages("SUMMER")
```

## Installation - Development Version

To download the development version of the SUMMER package, use the code below.

``` r
# install.packages("devtools")
devtools::install_github("richardli/SUMMER")
```
 
Examples of most of the main functions are described in the several vignettes listed on [https://cran.r-project.org/package=SUMMER](https://cran.r-project.org/package=SUMMER).

### Bug Reports / Change Requests
If you encounter a bug or would like make a change request, please file it as an issue [here](https://github.com/richardli/SUMMER/issues).


