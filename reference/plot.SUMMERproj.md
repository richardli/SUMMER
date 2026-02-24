# Plot projection output.

Plot projection output.

## Usage

``` r
# S3 method for class 'SUMMERproj'
plot(
  x,
  year.label = c("85-89", "90-94", "95-99", "00-04", "05-09", "10-14", "15-19"),
  year_label = deprecated(),
  year.med = c(1987, 1992, 1997, 2002, 2007, 2012, 2017),
  year_med = deprecated(),
  is.subnational = TRUE,
  year.proj = 2015,
  proj_year = deprecated(),
  data.add = NULL,
  option.add = list(point = NULL, lower = NULL, upper = NULL, by = NULL),
  color.add = "black",
  label.add = NULL,
  dodge.width = 0.5,
  plot.CI = NULL,
  per1000 = FALSE,
  color.CI = NULL,
  alpha.CI = 0.5,
  ...
)
```

## Arguments

- x:

  output from
  [`getSmoothed`](https://richardli.github.io/SUMMER/reference/getSmoothed.md)

- year.label:

  labels for the periods

- year_label:

  **\[deprecated\]** replaced by year.label

- year.med:

  labels for the middle years in each period, only used when both yearly
  and period estimates are plotted. In that case, `year.med` specifies
  where each period estimates are aligned.

- year_med:

  **\[deprecated\]** replaced by year.med

- is.subnational:

  logical indicator of whether the data contains subnational estimates

- year.proj:

  the first year where projections are made, i.e., where no data are
  available.

- proj_year:

  **\[deprecated\]** replaced by year.proj

- data.add:

  data frame for the Comparisons data points to add to the graph. This
  can be, for example, the raw direct estimates. This data frame is
  merged to the projections by column 'region' and 'years'. Except for
  these two columns, this dataset should not have Comparisons columns
  with names overlapping the getSmoothed output.

- option.add:

  list of options specifying the variable names for the points to plot,
  lower and upper bounds, and the grouping variable. This is intended to
  be used to add Comparisons estimates on the same plot as the smoothed
  estimates. See examples for details.

- color.add:

  the color of the Comparisons data points to plot.

- label.add:

  the label of the Comparisons data points in the legend.

- dodge.width:

  the amount to add to data points at the same year to avoid overlap.
  Default to be 0.5.

- plot.CI:

  logical indicator of whether to plot the error bars.

- per1000:

  logical indicator to plot mortality rates as rates per 1,000 live
  births. Note that the added comparison data should always be in the
  probability scale.

- color.CI:

  the color of the error bars of the credible interval.

- alpha.CI:

  the alpha (transparency) of the error bars of the credible interval.

- ...:

  optional arguments, see details

## See also

[`getSmoothed`](https://richardli.github.io/SUMMER/reference/getSmoothed.md)

## Author

Zehang Richard Li

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

#  national model
years.all <- c(years, "15-19")
fit1 <- smoothDirect(data = data, geo = NULL, Amat = NULL, 
  year.label = years.all, year.range = c(1985, 2019), 
  rw = 2, is.yearly=FALSE, m = 5)
out1 <- getSmoothed(fit1)
plot(out1, is.subnational=FALSE)

#  subnational model
fit2 <- smoothDirect(data = data, geo = geo, Amat = mat, 
  year.label = years.all, year.range = c(1985, 2019), 
  rw = 2, is.yearly=TRUE, m = 5, type.st = 4)
out2 <- getSmoothed(fit2)
plot(out2, is.yearly=TRUE, is.subnational=TRUE)


} # }
```
