# Aggregate estimators from different surveys.

Aggregate estimators from different surveys.

## Usage

``` r
aggregateSurvey(data)
```

## Arguments

- data:

  Output from
  [`getDirectList`](https://richardli.github.io/SUMMER/reference/getDirectList.md)

## Value

Estimators aggregated across surveys.

## Author

Zehang Richard Li

## Examples

``` r
if (FALSE) { # \dontrun{
data(DemoData)
data(DemoMap)
years <- levels(DemoData[[1]]$time)

# obtain direct estimates
data <- getDirectList(births = DemoData, 
years = years, 
regionVar = "region", timeVar = "time", 
clusterVar = "~clustid+id", 
ageVar = "age", weightsVar = "weights", 
geo.recode = NULL)

# obtain maps
geo <- DemoMap$geo
mat <- DemoMap$Amat

# Simulate hyper priors
priors <- simhyper(R = 2, nsamp = 1e+05, nsamp.check = 5000, Amat = mat, only.iid = TRUE)

# combine data from multiple surveys
data <- aggregateSurvey(data)
utils::head(data)

} # }
```
