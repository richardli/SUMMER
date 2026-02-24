# Extract posterior summaries of random effects

Extract posterior summaries of random effects

## Usage

``` r
getDiag(
  fitted,
  inla_mod = deprecated(),
  field = c("space", "time", "spacetime")[1],
  CI = 0.95,
  draws = NULL,
  nsim = 1000,
  ...
)
```

## Arguments

- fitted:

  output from
  [`smoothDirect`](https://richardli.github.io/SUMMER/reference/smoothDirect.md)
  or
  [`smoothCluster`](https://richardli.github.io/SUMMER/reference/smoothCluster.md)

- inla_mod:

  **\[deprecated\]** replaced by `fitted`.

- field:

  which random effects to plot. It can be one of the following: space,
  time, and spacetime.

- CI:

  Desired level of credible intervals

- draws:

  Posterior samples drawn from the fitted model. This argument allows
  the previously sampled draws (by setting save.draws to be TRUE) be
  used in new aggregation tasks.

- nsim:

  number of simulations, only applicable for the cluster-level model
  space-time interaction terms when random slopes are included.

- ...:

  Unused arguments, for users with fitted object from the package before
  v1.0.0, arguments including Amat, year.label, and year.range can still
  be specified manually.

## Value

List of diagnostic plots

## Author

Zehang Richard Li

## Examples

``` r
if (FALSE) { # \dontrun{
  data(DemoMap)
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
  fit1 <- smoothDirect(data = data, geo = DemoMap$geo, Amat = DemoMap$Amat, 
    year.label = years.all, year.range = c(1985, 2019), 
    rw = 2, is.yearly=FALSE, m = 5)
random.time <- getDiag(fit1, field = "time")
  random.space <- getDiag(fit1, field = "space")
  random.spacetime <- getDiag(fit1, field = "spacetime")
} # }
```
