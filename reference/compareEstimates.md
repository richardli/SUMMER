# Plot heatmap comparing pairwise posterior exceedence probabilities for svysae object

Plot heatmap comparing pairwise posterior exceedence probabilities for
svysae object

## Usage

``` r
compareEstimates(x, posterior.sample = NULL, title = NULL, return.plot = FALSE)
```

## Arguments

- x:

  an object in the S3 class of svysae, fhModel, or clusterModel. Plots
  are created for all models in this object.

- posterior.sample:

  Matrix of posteriors samples of area level quantities with one row for
  each area and one column for each sample. This argument may be
  specified to only provide a heatmap for the desired samples.

- title:

  Optional parameter changing the title of the plot

- return.plot:

  Logical indicator for whether the ggplot object is returned

## Value

ggplot containing heat map of pairwise comparisons

## Examples

``` r
if (FALSE) { # \dontrun{
data(DemoData2)
data(DemoMap2)
library(survey)
des0 <- svydesign(ids = ~clustid+id, strata = ~strata,
                  weights = ~weights, data = DemoData2, nest = TRUE)
Xmat <- aggregate(age~region, data = DemoData2, FUN = mean)

cts.res <- smoothArea(tobacco.use ~ 1,
                      domain = ~region,
                      design = des0,
                      adj.mat = DemoMap2$Amat, 
                      pc.u = 1,
                      pc.alpha = 0.01,
                      pc.u.phi = 0.5,
                      pc.alpha.phi = 2/3,
                      return.samples = TRUE)
compareEstimates(cts.res)
} # }
```
