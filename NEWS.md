# SUMMER - changes

Version 0.3.0 (2019-10-01) 
==========================
+ Version 0.3.0 contains some major updates from the previous versions. 
+ The following functions have been renamed. Most of the function arguments remain the same:
    * ``countrySummary`` is now ``getDirect``
    * ``cuontrySummary_mult`` is now ``getDirectList``
    * ``fitspace`` is now ``fitGeneric``
    * ``projINLA`` is now ``getSmooth``
+ The following new functions are added:
    + ``getDiag``: produce diagnostic plots for the fitted model.
    + ``getAdjusted``: produce adjusted estimates for a fitted model
    + ``getAmat``: automatic extract spatial adjacency matrix from the polygon file.
    + ``hatchPlot``: plot variables on a map with hatching indicating the width of the credible interval.
    + ``fitINLA2``: implements new smoothing methods based on binomial models at cluster level. 
+ Default prior changed to BYM2 + PC priors.
+ Various name changes for function arguments to reduce ambiguity. 
+ Various output column name change to improve consistency across functions