# Package index

## Child mortality estimation

- [`getBirths()`](https://richardli.github.io/SUMMER/reference/getBirths.md)
  : Reformat full birth records into person-month format

- [`getDirect()`](https://richardli.github.io/SUMMER/reference/getDirect.md)
  : Obtain the Horvitz-Thompson direct estimates and standard errors
  using delta method for a single survey.

- [`getDirectList()`](https://richardli.github.io/SUMMER/reference/getDirectList.md)
  : Obtain the Horvitz-Thompson direct estimates and standard errors
  using delta method for multiple surveys.

- [`smoothDirect()`](https://richardli.github.io/SUMMER/reference/smoothDirect.md)
  : Smoothed direct estimates for mortality rates

- [`smoothCluster()`](https://richardli.github.io/SUMMER/reference/smoothCluster.md)
  : Cluster-level space-time smoothing models for mortality rates

- [`getSmoothed()`](https://richardli.github.io/SUMMER/reference/getSmoothed.md)
  : Extract smoothed estimates.

- [`getAdjusted()`](https://richardli.github.io/SUMMER/reference/getAdjusted.md)
  : Adjust direct estimates and their associated variances

- [`getDiag()`](https://richardli.github.io/SUMMER/reference/getDiag.md)
  : Extract posterior summaries of random effects

- [`aggregateSurvey()`](https://richardli.github.io/SUMMER/reference/aggregateSurvey.md)
  : Aggregate estimators from different surveys.

- [`Benchmark()`](https://richardli.github.io/SUMMER/reference/Benchmark.md)
  **\[experimental\]** : Benchmark posterior draws to national estimates

- [`print(`*`<SUMMERmodel>`*`)`](https://richardli.github.io/SUMMER/reference/print.SUMMERmodel.md)
  : Print method for the smoothing models.

- [`summary(`*`<SUMMERmodel>`*`)`](https://richardli.github.io/SUMMER/reference/summary.SUMMERmodel.md)
  : Summary method for the smoothing models.

- [`print(`*`<SUMMERprojlist>`*`)`](https://richardli.github.io/SUMMER/reference/print.SUMMERprojlist.md)
  : Print method for the combined projection output.

- [`summary(`*`<SUMMERprojlist>`*`)`](https://richardli.github.io/SUMMER/reference/summary.SUMMERprojlist.md)
  :

  Summary method for the combined projection output. This function is
  the print method for class `SUMMERprojlist`.

## General SAE models

- [`smoothSurvey()`](https://richardli.github.io/SUMMER/reference/smoothSurvey.md)
  : Fit space-time smoothing models for a binary outcome from complex
  surveys.

- [`smoothArea()`](https://richardli.github.io/SUMMER/reference/smoothArea.md)
  : Small area estimation via basic area level model

- [`smoothUnit()`](https://richardli.github.io/SUMMER/reference/smoothUnit.md)
  : Smooth via basic unit level model

- [`print(`*`<SUMMERmodel.svy>`*`)`](https://richardli.github.io/SUMMER/reference/print.SUMMERmodel.svy.md)
  :

  Print method for the smoothing models from `smoothSurvey`.

- [`summary(`*`<SUMMERmodel.svy>`*`)`](https://richardli.github.io/SUMMER/reference/summary.SUMMERmodel.svy.md)
  :

  Summary method for the smoothing model and output from `smoothSurvey`.

## Visualization

- [`mapPlot()`](https://richardli.github.io/SUMMER/reference/mapPlot.md)
  : Plot region-level variables on a map
- [`ridgePlot()`](https://richardli.github.io/SUMMER/reference/ridgePlot.md)
  : Calculate and plot posterior densities of the projected estimates
- [`hatchPlot()`](https://richardli.github.io/SUMMER/reference/hatchPlot.md)
  : Plot maps with uncertainty hatching.
- [`plot(`*`<SUMMERproj>`*`)`](https://richardli.github.io/SUMMER/reference/plot.SUMMERproj.md)
  : Plot projection output.
- [`tcpPlot()`](https://richardli.github.io/SUMMER/reference/tcpPlot.md)
  : Discrete-color maps based on the True Classification Probabilities
- [`compareEstimates()`](https://richardli.github.io/SUMMER/reference/compareEstimates.md)
  : Plot heatmap comparing pairwise posterior exceedence probabilities
  for svysae object
- [`mapEstimates()`](https://richardli.github.io/SUMMER/reference/mapEstimates.md)
  : Mapping estimates for svysae object

## Utility functions

- [`mapPoints()`](https://richardli.github.io/SUMMER/reference/mapPoints.md)
  : Map GPS points to polygon regions
- [`getAmat()`](https://richardli.github.io/SUMMER/reference/getAmat.md)
  : Extract adjacency matrix from the map
- [`getCounts()`](https://richardli.github.io/SUMMER/reference/getCounts.md)
  : Aggregate person-month data into counts and totals by groups.
- [`changeRegion()`](https://richardli.github.io/SUMMER/reference/ChangeRegion.md)
  : Map region names to a common set.
- [`expit()`](https://richardli.github.io/SUMMER/reference/expit.md) :
  Expit transformation
- [`logit()`](https://richardli.github.io/SUMMER/reference/logit.md) :
  Logit transformation
- [`rst()`](https://richardli.github.io/SUMMER/reference/rst.md) :
  Simulate spatial and temporal random effects

## Pixel-level population simulation

- [`getAreaName()`](https://richardli.github.io/SUMMER/reference/getAreaName.md)
  **\[experimental\]** : Determines which administrative areas contain
  the given points

- [`pixelPopToArea()`](https://richardli.github.io/SUMMER/reference/aggPop.md)
  [`areaPopToArea()`](https://richardli.github.io/SUMMER/reference/aggPop.md)
  : Aggregate populations to the specified areal level

- [`simPopSPDE()`](https://richardli.github.io/SUMMER/reference/simPop.md)
  [`simPopCustom()`](https://richardli.github.io/SUMMER/reference/simPop.md)
  : Simulate populations and areal prevalences

- [`aggPixelPreds()`](https://richardli.github.io/SUMMER/reference/aggPixelPreds.md)
  :

  Helper function of `pixelPopToArea`

- [`poppRegionFromPopMat()`](https://richardli.github.io/SUMMER/reference/poppRegionFromPopMat.md)
  **\[experimental\]** :

  Generate a population frame of a similar format to poppa argument of
  `simPopCustom` with a custom set of regions

- [`simSPDE()`](https://richardli.github.io/SUMMER/reference/simSPDE.md)
  : Simulate from the SPDE spatial model

- [`getExpectedNperEA()`](https://richardli.github.io/SUMMER/reference/simPopInternal.md)
  [`getSortIndices()`](https://richardli.github.io/SUMMER/reference/simPopInternal.md)
  [`rStratifiedMultnomial()`](https://richardli.github.io/SUMMER/reference/simPopInternal.md)
  [`rStratifiedMultnomialBySubarea()`](https://richardli.github.io/SUMMER/reference/simPopInternal.md)
  [`rMyMultinomial()`](https://richardli.github.io/SUMMER/reference/simPopInternal.md)
  [`rMyMultinomialSubarea()`](https://richardli.github.io/SUMMER/reference/simPopInternal.md)
  [`rmultinom1()`](https://richardli.github.io/SUMMER/reference/simPopInternal.md)
  [`sampleNMultilevelMultinomial()`](https://richardli.github.io/SUMMER/reference/simPopInternal.md)
  [`sampleNMultilevelMultinomialFixed()`](https://richardli.github.io/SUMMER/reference/simPopInternal.md)
  : Internal functions for population simulation

- [`makePopIntegrationTab()`](https://richardli.github.io/SUMMER/reference/makePopIntegrationTab.md)
  [`getPoppsub()`](https://richardli.github.io/SUMMER/reference/makePopIntegrationTab.md)
  [`adjustPopMat()`](https://richardli.github.io/SUMMER/reference/makePopIntegrationTab.md)
  **\[experimental\]** : Generating pixellated populations, and
  population frames

- [`calibrateByRegion()`](https://richardli.github.io/SUMMER/reference/calibrateByRegion.md)
  : Calibrate the point level totals so their sum matches the regional
  totals

- [`setThresholdsByRegion()`](https://richardli.github.io/SUMMER/reference/setThresholdsByRegion.md)
  :

  **\[experimental\]**

- [`projKenya()`](https://richardli.github.io/SUMMER/reference/projKenya.md)
  **\[experimental\]** : Map projection for Kenya

## Dataset

- [`BRFSS`](https://richardli.github.io/SUMMER/reference/BRFSS.md) : The
  BRFSS dataset
- [`DemoData`](https://richardli.github.io/SUMMER/reference/DemoData.md)
  : Simulated child mortality person-month dataset.
- [`DemoData2`](https://richardli.github.io/SUMMER/reference/DemoData2.md)
  : Simulated dataset for prevalence mapping.
- [`DemoMap`](https://richardli.github.io/SUMMER/reference/DemoMap.md) :
  Uganda Admin-1 region map for illustration purpose
- [`DemoMap2`](https://richardli.github.io/SUMMER/reference/DemoMap2.md)
  : Kenya Admin-1 region map for illustration purpose
- [`KingCounty`](https://richardli.github.io/SUMMER/reference/KingCounty.md)
  : Map of King County
- [`MalawiData`](https://richardli.github.io/SUMMER/reference/MalawiData.md)
  : Auxiliary data for Malawi 2000, 2004, 2010, and 2015 DHS.
- [`MalawiMap`](https://richardli.github.io/SUMMER/reference/MalawiMap.md)
  : Malawi Admin-2 map
- [`KenData`](https://richardli.github.io/SUMMER/reference/KenData.md) :
  Auxiliary data for Kenya 2014 DHS.
- [`kenyaPopulationData`](https://richardli.github.io/SUMMER/reference/kenyaPopulationData.md)
  [`easpaKenyaNeonatal`](https://richardli.github.io/SUMMER/reference/kenyaPopulationData.md)
  [`poppaKenya`](https://richardli.github.io/SUMMER/reference/kenyaPopulationData.md)
  [`poppsubKenya`](https://richardli.github.io/SUMMER/reference/kenyaPopulationData.md)
  : Kenya 2009 Census Frame and Related Datasets

## Internal functions

- [`rw.new()`](https://richardli.github.io/SUMMER/reference/rw.new.md) :
  New random walk 1 and 2 models for m-year to period random effects
- [`rw.new.pc()`](https://richardli.github.io/SUMMER/reference/rw.new.pc.md)
  : New random walk 1 and 2 models for m-year to period random effects
- [`st.new()`](https://richardli.github.io/SUMMER/reference/st.new.md) :
  New Type I to IV space time interaction models for m-year to period
  random effects
- [`st.new.pc()`](https://richardli.github.io/SUMMER/reference/st.new.pc.md)
  : New Type I to IV space time interaction models for m-year to period
  random effects
- [`iid.new()`](https://richardli.github.io/SUMMER/reference/iid.new.md)
  : New random IID models for m-year to period random effects
- [`iid.new.pc()`](https://richardli.github.io/SUMMER/reference/iid.new.pc.md)
  : New random IID models for m-year to period random effects
- [`logitNormMean()`](https://richardli.github.io/SUMMER/reference/logitNormMean.md)
  : Calculate the mean of a distribution whose logit is Gaussian
- [`simhyper()`](https://richardli.github.io/SUMMER/reference/simhyper.md)
  : Simulate hyperpriors from an GMRF
- [`SUMMER-package`](https://richardli.github.io/SUMMER/reference/SUMMER-package.md)
  [`SUMMER`](https://richardli.github.io/SUMMER/reference/SUMMER-package.md)
  : SUMMER: Small-Area-Estimation Unit/Area Models and Methods for
  Estimation in R
