# SUMMER - changes

Version 0.3.0 (2019-10-01) 
==========================
+ Version 0.3.0 contains some major updates from the previous versions. 
+ The following functions have been renamed. Most of the function arguments remain the same:
    * ``countrySummary`` is now ``getDirect``
    * ``cuontrySummary_mult`` is now ``getDirectList``
    * ``fitspace`` is now ``fitGeneric``
    * ``projINLA`` is now ``getSmooth``   
+ Major new functions:
    + ``fitINLA2``: implements new smoothing methods based on binomial models at cluster level. 
    + ``getDiag``: produce diagnostic plots for the fitted model.
    + ``hatchPlot``: plot variables on a map with hatching indicating the width of the credible interval.
+ The following help functions are added:
    + ``getAdjusted``: produce adjusted estimates for a fitted model
    + ``getAmat``: automatic extract spatial adjacency matrix from the polygon file.
    + ``getCounts``: aggregate person-month data into counts and totals by groups.
    + ``mapPoints``: map GPS points to polygon regions.    
+ Default prior changed to BYM2 + PC priors.
+ Various name changes for function arguments to reduce ambiguity. 
+ Various output column name change to improve consistency across functions.
+ Three set of new vignettes are provided to guide users through various modeling workflow.