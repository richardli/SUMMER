# Discrete-color maps based on the True Classification Probabilities

Discrete-color maps based on the True Classification Probabilities

## Usage

``` r
tcpPlot(
  draws,
  geo,
  by.geo = NULL,
  year.plot = NULL,
  year_plot = deprecated(),
  ncol = 4,
  per1000 = FALSE,
  thresholds = NULL,
  intervals = 3,
  size.title = 0.7,
  legend.label = NULL,
  border = "gray20",
  size = 0.5
)
```

## Arguments

- draws:

  a posterior draw object from
  [`getSmoothed`](https://richardli.github.io/SUMMER/reference/getSmoothed.md)

- geo:

  SpatialPolygonsDataFrame object for the map

- by.geo:

  variable name specifying region names in geo

- year.plot:

  vector of year string vector to be plotted.

- year_plot:

  **\[deprecated\]** replaced by year.plot

- ncol:

  number of columns in the output figure.

- per1000:

  logical indicator to multiply results by 1000.

- thresholds:

  a vector of thresholds (on the mortality scale) defining the discrete
  color scale of the maps.

- intervals:

  number of quantile intervals defining the discrete color scale of the
  maps. Required when thresholds are not specified.

- size.title:

  a numerical value giving the amount by which the plot title should be
  magnified relative to the default.

- legend.label:

  Label for the color legend.

- border:

  color of the border

- size:

  size of the border

## Value

a list of True Classification Probability (TCP) tables, a list of
individual spplot maps, and a gridded array of all maps.

## References

Tracy Qi Dong, and Jon Wakefield. (2020) *Modeling and presentation of
vaccination coverage estimates using data from household surveys.* arXiv
preprint arXiv:2004.03127.

## Author

Tracy Qi Dong, Zehang Richard Li

## Examples

``` r
if (FALSE) { # \dontrun{
library(dplyr)
data(DemoData)
# Create dataset of counts, unstratified
counts.all <- NULL
for(i in 1:length(DemoData)){
  counts <- getCounts(DemoData[[i]][, c("clustid", "time", "age", "died",
                                        "region")],
            variables = 'died', by = c("age", "clustid", "region", 
                                         "time"))
  counts <- counts %>% mutate(cluster = clustid, years = time, Y=died)
  counts$strata <- NA
  counts$survey <- names(DemoData)[i] 
  counts.all <- rbind(counts.all, counts)
}

# fit cluster-level model on the periods
periods <- levels(DemoData[[1]]$time)
fit <- smoothCluster(data = counts.all, 
      Amat = DemoMap$Amat, 
      time.model = "rw2", 
      st.time.model = "rw1",
      strata.time.effect =  TRUE, 
      survey.effect = TRUE,
      family = "betabinomial",
      year.label = c(periods, "15-19"))
est <- getSmoothed(fit, nsim = 1000, save.draws=TRUE)

tcp <- tcpPlot(est, DemoMap$geo, by.geo = "REGNAME", interval = 3, year.plot = periods) 
tcp$g
} # }
```
