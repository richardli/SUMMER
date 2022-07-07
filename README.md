# SUMMER
[![](https://github.com/richardli/SUMMER/workflows/R-CMD-check_INLA_stable/badge.svg)](https://github.com/richardli/SUMMER/actions) [![](https://github.com/richardli/SUMMER/workflows/R-CMD-check_INLA_testing/badge.svg)](https://github.com/richardli/SUMMER/actions) 
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/SUMMER)](https://cran.r-project.org/package=SUMMER) [![](https://cranlogs.r-pkg.org/badges/SUMMER)](https://cran.r-project.org/package=SUMMER) [![](https://cranlogs.r-pkg.org/badges/grand-total/SUMMER?color=orange)](https://cran.r-project.org/package=SUMMER)

SAE Unit/area Models and Methods for Estimation in R

## Next major update (version 1.4.0)
+ New major functions for simulation and using SPDE spatial model and aggregation over grids.
+ New datasets on Kenya population map.
+ New methods for benchmarking.

## Major update (version 1.3.0)
+ New SAE functions and vignette.

## Major update (version 1.2.0)
+ Major expansion for `smoothSurvey` to implement popular SAE methods. Syntax change to the function. 
+ Allows `smoothDirect` and `smoothCluster` to fit space-only models.
+ S3 class output.

## Major update (version 1.1.0)
+ Major expansion for `smoothSurvey` to implement popular SAE methods. Syntax change to the function. 
+ Allows `smoothDirect` and `smoothCluster` to fit space-only models.
+ `smoothSurvey`, `smoothDirect`, and `smoothCluster` now returns S3 classed objects.

## Major update (version 1.0.0)
Version 1.0.0 contains many major updates from the previous versions:
+ Major updates to functions.
    + ``fitGeneric`` is now ``smoothSurvey``
    + ``fitINLA`` is now ``smoothDirect``
    + ``fitINLA2`` is now ``smoothCluster``
    + More extensions in both smoothed direct and cluster level models.
    + Major changes to how temporal models are specified with `time.model` and `st.time.model`.
    + More interpretable parameterization of slope and random slopes.
    + More visualization options.
    + Note: Previous function name and argument syntax remain to work as before, but may not receive high priority in maintenance in the future.  
+ Many minor improvements in functions.
    + Better model summary message.
    + Removed unnecessary function arguments, e.g., ``geo`` in various functions.
    + Removed the requirement to repeated specifying ``Amat``, ``year_label`` and ``year_range``. Now they are only required in the model fitting stage.
+ New vignettes. 

## Major update (version 0.3.0)
Version 0.3.0 contains some major updates from the previous versions. Some of the substantial changes to existing functions are listed here. For a complete log of changes, see the [News](https://github.com/richardli/SUMMER/blob/master/NEWS.md) section.

+ Function name changes
    * ``countrySummary`` is now ``getDirect``
    * ``cuontrySummary_mult`` is now ``getDirectList``
    * ``fitspace`` is now ``fitGeneric``
    * ``projINLA`` is now ``getSmooth``
+ New functions
    * ``getDiag``: produce diagnostic plots for the fitted model.
    * ``getAdjusted``: produce adjusted estimates for a fitted model.
    * ``getAmat``: automatic extract spatial adjacency matrix from the polygon file.
    * ``hatchPlot``: plot variables on a map with hatching indicating the width of the credible interval.
+New methods
    * ``fitINLA2``: implements new smoothing methods based on binomial models at cluster level. 

## Citation

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


## Installation - CRAN

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

## Bug Reports / Change Requests
If you encounter a bug or would like make a change request, please file it as an issue [here](https://github.com/richardli/SUMMER/issues).


