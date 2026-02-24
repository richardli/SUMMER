# Calculate and plot posterior densities of the projected estimates

The function `ridgePlot` replaces the previous function name
`getSmoothedDensity` (before version 1.0.0).

## Usage

``` r
ridgePlot(
  x = NULL,
  nsim = 1000,
  draws = NULL,
  year.plot = NULL,
  year_plot = deprecated(),
  strata.plot = NULL,
  strata_plot = deprecated(),
  by.year = TRUE,
  ncol = 4,
  scale = 2,
  per1000 = FALSE,
  order = 0,
  direction = 1,
  linewidth = 0.5,
  results = NULL,
  save.density = FALSE,
  ...
)
```

## Arguments

- x:

  output from
  [`smoothDirect`](https://richardli.github.io/SUMMER/reference/smoothDirect.md)
  for the smoothed direct estimates, or
  [`smoothCluster`](https://richardli.github.io/SUMMER/reference/smoothCluster.md)
  for the cluster-level estimates.

- nsim:

  number of posterior draws to take. Only used for cluster-level models
  when `draws` is NULL. Otherwise the posterior draws in `draws` will be
  used instead without resampling.

- draws:

  Output of
  [`getSmoothed`](https://richardli.github.io/SUMMER/reference/getSmoothed.md)
  with `save.draws` set to TRUE. This argument allows the previously
  sampled draws (by setting `save.draws` to be TRUE) be used in new
  aggregation tasks. This argument is only used for cluster-level
  models.

- year.plot:

  A vector indicate which years to plot

- year_plot:

  **\[deprecated\]** replaced by year.plot

- strata.plot:

  Name of the strata to plot. If not specified, the overall is plotted.

- strata_plot:

  **\[deprecated\]** replaced by strata.plot

- by.year:

  logical indicator for whether the output uses years as facets.

- ncol:

  number of columns in the output figure.

- scale:

  numerical value controlling the height of the density plots.

- per1000:

  logical indicator to multiply results by 1000.

- order:

  order of regions when by.year is set to TRUE. Negative values indicate
  regions are ordered from high to low posterior medians from top to
  bottom. Positive values indicate from low to high. 0 indicate
  alphabetic orders.

- direction:

  Direction of the color scheme. It can be either 1 (smaller values are
  darker) or -1 (higher values are darker). Default is set to 1.

- linewidth:

  width of the ridgeline.

- results:

  output from `ridgePlot` returned object with `save.density = TRUE`.
  This argument can be specified to avoid calculating densities again
  when only the visualization changes.

- save.density:

  Logical indicator of whether the densities will be returned with the
  ggplot object. If set to TRUE, the output will be a list consisting
  of (1) a data frame of computed densities and (2) a ggplot object of
  the plot.

- ...:

  additional configurations passed to inla.posterior.sample.

## Value

ridge plot of the density, and if `save.density = TRUE`, also a data
frame of the calculated densities

## See also

[`plot.SUMMERproj`](https://richardli.github.io/SUMMER/reference/plot.SUMMERproj.md)

## Author

Zehang Richard Li

## Examples

``` r
if (FALSE) { # \dontrun{
years <- levels(DemoData[[1]]$time)

data <- getDirectList(births = DemoData, 
years = years,  
regionVar = "region", timeVar = "time", 
clusterVar = "~clustid+id", 
ageVar = "age", weightsVar = "weights", 
geo.recode = NULL)
# obtain direct estimates
data_multi <- getDirectList(births = DemoData, years = years,
  regionVar = "region",  timeVar = "time", clusterVar = "~clustid+id",
  ageVar = "age", weightsVar = "weights", geo.recode = NULL)
data <- aggregateSurvey(data_multi)

#  national model
years.all <- c(years, "15-19")
fit1 <- smoothDirect(data = data, geo = NULL, Amat = NULL, 
  year.label = years.all, year.range = c(1985, 2019), 
  rw = 2, m = 5)
## Plot marginal posterior densities over time
ridgePlot(fit1, year.plot = years.all, 
          ncol = 4, by.year = FALSE)

#  subnational model
fit2 <- smoothDirect(data = data, geo = DemoMap$geo, Amat = DemoMap$Amat, 
  year.label = years.all, year.range = c(1985, 2019), 
  rw = 2, m = 5, type.st = 1)

# Plot marginal posterior densities over time (regions are ordered alphabetically)
ridgePlot(fit2, year.plot = years.all, ncol = 4)

# Re-order the regions and save the density to avoid re-compute later
density <- ridgePlot(fit2, year.plot = years.all,
 ncol = 4, per1000 = TRUE, order = -1, save.density = TRUE)
density$g

# Show each region (instead of each year) in a panel 
## Instead of recalculate the posteriors, we can use previously calculated densities as input 
ridgePlot(results = density, year.plot = years.all, 
ncol = 4, by.year=FALSE, per1000 = TRUE)

# Show more years
ridgePlot(results = density, year.plot = c(1990:2019), 
ncol = 4, by.year=FALSE, per1000 = TRUE)


# Example using surveyPrev package output

library(surveyPrev)
dhsData <- getDHSdata(country = "Rwanda", indicator = "nmr", year = 2019)
data <- getDHSindicator(dhsData, indicator = "nmr")
geo <- getDHSgeo(country = "Rwanda", year = 2019)
poly.adm1 <- geodata::gadm(country="RWA", level=1, path=tempdir())
poly.adm1 <- sf::st_as_sf(poly.adm1)
poly.adm2 <- geodata::gadm(country="RWA", level=2, path=tempdir())
poly.adm2 <- sf::st_as_sf(poly.adm2)
cluster.info <- clusterInfo(geo = geo, 
              poly.adm1 = poly.adm1, 
              poly.adm2 = poly.adm2,
                            by.adm1 = "NAME_1", 
                            by.adm2 = "NAME_2")

fit1 <- directEST(data = data, cluster.info = cluster.info,  admin = 1)
fit2 <- directEST(data = data, cluster.info = cluster.info,  admin = 2) 
ridgePlot(fit1, direction = -1)
ridgePlot(fit2, direction = -1)

} # }
```
