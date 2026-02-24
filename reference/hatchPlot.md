# Plot maps with uncertainty hatching.

This function visualizes the map with different variables. The input
data frame can be either the long or wide format.

## Usage

``` r
hatchPlot(
  data,
  variables,
  values = NULL,
  labels = NULL,
  geo,
  by.data,
  by.geo,
  is.long = FALSE,
  lower,
  upper,
  lim = NULL,
  lim.CI = NULL,
  breaks.CI = NULL,
  ncol = 4,
  hatch = NULL,
  border = NULL,
  size = 1,
  legend.label = NULL,
  per1000 = FALSE,
  direction = 1,
  ...
)
```

## Arguments

- data:

  a data frame with variables to be plotted

- variables:

  vector of variables to be plotted. If long format of data is used,
  only one variable can be selected

- values:

  the column corresponding to the values to be plotted, only used when
  long format of data is used

- labels:

  vector of labels to use for each variable, only used when wide format
  of data is used

- geo:

  SpatialPolygonsDataFrame object for the map

- by.data:

  column name specifying region names in the data

- by.geo:

  variable name specifying region names in the data

- is.long:

  logical indicator of whether the data is in the long format, default
  to FALSE

- lower:

  column name of the lower bound of the CI

- upper:

  column name of the upper bound of the CI

- lim:

  fixed range of values for the variables to plot

- lim.CI:

  fixed range of the CI widths to plot

- breaks.CI:

  a vector of numerical values that decides the breaks in the CI widths
  to be shown

- ncol:

  number of columns for the output tabs

- hatch:

  color of the hatching lines.

- border:

  color of the polygon borders.

- size:

  line width of the polygon borders.

- legend.label:

  Label for the color legend.

- per1000:

  logical indicator to plot mortality rates as rates per 1,000 live
  births. Note that the added comparison data should always be in the
  probability scale.

- direction:

  Direction of the color scheme. It can be either 1 (smaller values are
  darker) or -1 (higher values are darker). Default is set to 1.

- ...:

  unused.

## Author

Zehang Richard Li, Katie Wilson

## Examples

``` r
if (FALSE) { # \dontrun{
years <- levels(DemoData[[1]]$time)

# obtain direct estimates
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

fit2 <- smoothDirect(data = data, geo = geo, Amat = mat, 
  year.label = years.all, year.range = c(1985, 2019), 
  rw = 2, is.yearly=TRUE, m = 5, type.st = 4)
out2 <- getSmoothed(fit2)

plot(out2, is.yearly=TRUE, is.subnational=TRUE)

hatchPlot(data = subset(out2, is.yearly==FALSE), geo = geo,
variables=c("years"), values = c("median"), 
by.data = "region", by.geo = "REGNAME", 
lower = "lower", upper = "upper", is.long=TRUE)

} # }
```
