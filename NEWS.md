# SUMMER - changes
Version 1.0.0 (2020-07-01) 
==========================
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

Version 0.3.1 (2019-10-23) 
==========================
+ Fixed minor typo that lower and upper bounds switched after applying aggregateSurvey().
+ Add option to remove person-month records if last time period contain only k years of data in getBirth().
+ Improved mapPlot().

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