# Smoothed direct estimates for mortality rates

The function `smoothDirect` replaces the previous function name
`fitINLA` (before version 1.0.0).

## Usage

``` r
smoothDirect(
  data,
  Amat,
  formula = NULL,
  time.model = c("rw1", "rw2", "ar1")[2],
  st.time.model = NULL,
  year.label,
  year_label = deprecated(),
  year.range = c(1980, 2014),
  year_range = deprecated(),
  is.yearly = TRUE,
  m = 5,
  type.st = 1,
  survey.effect = FALSE,
  hyper = c("pc", "gamma")[1],
  pc.u = 1,
  pc.alpha = 0.01,
  pc.u.phi = 0.5,
  pc.alpha.phi = 2/3,
  pc.u.cor = 0.7,
  pc.alpha.cor = 0.9,
  pc.st.u = NA,
  pc.st.alpha = NA,
  control.compute = list(dic = TRUE, mlik = TRUE, cpo = TRUE, openmp.strategy =
    "default", config = TRUE),
  control.inla = list(strategy = "adaptive", int.strategy = "auto"),
  control.fixed = list(),
  verbose = FALSE,
  geo = NULL,
  rw = NULL,
  ar = NULL,
  options = NULL
)
```

## Arguments

- data:

  Combined dataset

- Amat:

  Adjacency matrix for the regions

- formula:

  INLA formula. See vignette for example of using customized formula.

- time.model:

  Model for the main temporal trend, can be rw1, rw2, or ar1. ar1 is not
  implemented for yearly model with period data input. Default to be
  rw2. For ar1 main effect, a linear slope is also added with time
  scaled to be between -0.5 to 0.5, i.e., the slope coefficient
  represents the total change between the first year and the last year
  in the projection period on the logit scale.

- st.time.model:

  Temporal component model for the interaction term, can be rw1, rw2, or
  ar1. ar1 is not implemented for yearly model with period data input.
  Default to be the same as time.model unless specified otherwise. For
  ar1 interaction model, region-specific random slopes are currently not
  implemented.

- year.label:

  string vector of year names

- year_label:

  **\[deprecated\]** replaced by year.label

- year.range:

  Entire range of the years (inclusive) defined in year.label.

- year_range:

  **\[deprecated\]** replaced by year.range

- is.yearly:

  Logical indicator for fitting yearly or period model.

- m:

  Number of years in each period.

- type.st:

  type for space-time interaction

- survey.effect:

  logical indicator whether to include a survey iid random effect. If
  this is set to TRUE, there needs to be a column named 'survey' in the
  input data frame. In prediction, this random effect term will be set
  to 0. Notice this survey effect is implemented according to the Merter
  et al. (2015) model, and differently compared to the smoothCluster()
  function.

- hyper:

  which hyperpriors to use. Default to be using the PC prior ("pc").

- pc.u:

  hyperparameter U for the PC prior on precisions.

- pc.alpha:

  hyperparameter alpha for the PC prior on precisions.

- pc.u.phi:

  hyperparameter U for the PC prior on the mixture probability phi in
  BYM2 model.

- pc.alpha.phi:

  hyperparameter alpha for the PC prior on the mixture probability phi
  in BYM2 model.

- pc.u.cor:

  hyperparameter U for the PC prior on the autocorrelation parameter in
  the AR prior, i.e. Prob(cor \> pc.u.cor) = pc.alpha.cor.

- pc.alpha.cor:

  hyperparameter alpha for the PC prior on the autocorrelation parameter
  in the AR prior.

- pc.st.u:

  hyperparameter U for the PC prior on precisions for the interaction
  term.

- pc.st.alpha:

  hyperparameter alpha for the PC prior on precisions for the
  interaction term.

- control.compute:

  list of options to be passed to control.compute() in the inla()
  function. The default argument saves the internal objects created by
  INLA for posterior sampling later. If the fitted object is too large
  in size and there is no need to perform joint posterior sampling from
  the model (only used in benchmarking), this argument can be set to
  `control.compute = list(config = FALSE)` to reduce the size of the
  fitted object.

- control.inla:

  list of options to be passed to control.inla() in the inla() function.
  Default to the "adaptive" integration strategy.

- control.fixed:

  list of options to be passed to control.fixed() in the inla()
  function.

- verbose:

  logical indicator to print out detailed inla() intermediate steps.

- geo:

  Deprecated.

- rw:

  Deprecated.

- ar:

  Deprecated.

- options:

  Deprecated.

## Value

List of fitted object

## References

Li, Z., Hsiao, Y., Godwin, J., Martin, B. D., Wakefield, J., Clark, S.
J., & with support from the United Nations Inter-agency Group for Child
Mortality Estimation and its technical advisory group. (2019). *Changes
in the spatial distribution of the under-five mortality rate: Small-area
analysis of 122 DHS surveys in 262 subregions of 35 countries in
Africa.* PloS one, 14(1), e0210645.

Mercer, L. D., Wakefield, J., Pantazis, A., Lutambi, A. M., Masanja, H.,
& Clark, S. (2015). *Space-time smoothing of complex survey data: small
area estimation for child mortality.* The annals of applied statistics,
9(4), 1889.

## See also

[`getDirect`](https://richardli.github.io/SUMMER/reference/getDirect.md)

## Author

Zehang Richard Li

## Examples

``` r
if (FALSE) { # \dontrun{
  years <- levels(DemoData[[1]]$time)
  # obtain direct estimates
  data_multi <- getDirectList(births = DemoData, years = years,
  regionVar = "region",  timeVar = "time", clusterVar = "~clustid+id",
  ageVar = "age", weightsVar = "weights", geo.recode = NULL)
  data <- aggregateSurvey(data_multi)
  
  #  national model
  years.all <- c(years, "15-19")
  fit1 <- smoothDirect(data = data, Amat = NULL, 
  year.label = years.all, year.range = c(1985, 2019), 
  time.model = 'rw2', m = 5, control.compute = list(config =TRUE))
  out1 <- getSmoothed(fit1)
  plot(out1)
  
  #  subnational model
  fit2 <- smoothDirect(data = data, Amat = DemoMap$Amat, 
  year.label = years.all, year.range = c(1985, 2019), 
  time.model = 'rw2', m = 5, type.st = 4)
  out2 <- getSmoothed(fit2)
  plot(out2)
  
  #  subnational space-only model for one period
  fit3 <- smoothDirect(data = subset(data, years == "10-14"), 
           time.model = NULL, Amat = DemoMap$Amat)
  out3 <- getSmoothed(fit3)
  plot(out3, plot.CI = TRUE)
} # }
```
