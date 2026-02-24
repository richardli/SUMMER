# Cluster-level space-time smoothing models for mortality rates

The function `smoothCluster` replace the previous function name
`fitINLA2` (before version 1.0.0).

## Usage

``` r
smoothCluster(
  data,
  X = NULL,
  family = c("betabinomial", "binomial")[1],
  age.group = c("0", "1-11", "12-23", "24-35", "36-47", "48-59"),
  age.groups = deprecated(),
  age.n = c(1, 11, 12, 12, 12, 12),
  age.time.group = c(1, 2, 3, 3, 3, 3),
  age.strata.fixed.group = c(1, 2, 3, 4, 5, 6),
  time.model = c("rw1", "rw2", "ar1")[2],
  st.time.model = NULL,
  Amat,
  bias.adj = NULL,
  bias.adj.by = NULL,
  formula = NULL,
  year.label,
  year_label = deprecated(),
  type.st = 4,
  survey.effect = FALSE,
  linear.trend = TRUE,
  common.trend = FALSE,
  strata.time.effect = FALSE,
  hyper = "pc",
  pc.u = 1,
  pc.alpha = 0.01,
  pc.u.phi = 0.5,
  pc.alpha.phi = 2/3,
  pc.u.cor = 0.7,
  pc.alpha.cor = 0.9,
  pc.st.u = NA,
  pc.st.alpha = NA,
  pc.st.slope.u = NA,
  pc.st.slope.alpha = NA,
  overdisp.mean = 0,
  overdisp.prec = 0.4,
  options = list(config = TRUE),
  control.inla = list(strategy = "adaptive", int.strategy = "auto"),
  control.fixed = list(),
  verbose = FALSE,
  geo = NULL,
  rw = NULL,
  ar = NULL,
  st.rw = NULL,
  age.rw.group = NULL,
  ...
)
```

## Arguments

- data:

  count data of person-months with the following columns

  - cluster: cluster ID

  - years: time period

  - region: region of the cluster

  - strata: stratum of the cluster

  - age: age group corresponding to the row

  - total: total number of person-month in this age group, stratum,
    cluster, and period

  - Y: total number of deaths in this age group, stratum, cluster, and
    period

- X:

  Covariate matrix. It must contain either a column with name "region",
  or a column with name "years", or both. The covariates must not have
  missing values for all regions (if varying in space) and all time
  periods (if varying in time). The rest of the columns are treated as
  covariates in the mean model.

- family:

  family of the model. This can be either binomial (with logistic normal
  prior), betabiniomial.

- age.group:

  a character vector of age groups in increasing order.

- age.groups:

  **\[deprecated\]** replaced by `age.group`

- age.n:

  number of months in each age groups in the same order.

- age.time.group:

  vector indicating grouping of the ages groups in the temporal model.
  For example, if each age group is assigned a different temporal
  component, then set age.rw.group to c(1:length(age.group)); if all age
  groups share the same random walk component, then set age.rw.group to
  a rep(1, length(age.group)). The default for 6 age groups is
  c(1,2,3,3,3,3), which assigns a separate temporal trend to the first
  two groups and a common random walk for the rest of the age groups.
  The vector should contain values starting from 1. This argument
  replaces the previous `age.rw.group` argument.

- age.strata.fixed.group:

  vector indicating grouping of the ages groups for different strata in
  the intercept. The default is c(1:length(age.group)), which correspond
  to each age group within each stratum receives a separate intercept.
  If several age groups are specified to be the same value in this
  vector, the stratum specific deviation from the baseline is assumed to
  be the same for these age groups. For example, if
  `age.strata.fixed.group = c(1, 2, 3, 3, 3, 3)`, then the intercept
  part of the linear predictor consists of 6 overall age-specific
  intercepts and 3 set of strata effects (where a base stratum is chosen
  internally), for age groups 1, 2, and the rest respectively. Notice
  that this argument does not control the linear temporal trends (which
  is also parameterized as fixed effect, but determined by the
  `age.rw.group` argument). The vector should contain values starting
  from 1.

  More specific examples: (1) if each age group is assigned a different
  intercept, then set age.strata.fixed.group to
  c(1:length(age.group)) (2) if all age groups share the same intercept,
  then set age.strata.fixed.group to a rep(1, length(age.group)). The
  default for 6 age groups is the former. (3) If each temporal trend is
  associated with its own intercept, set it to be the same as
  `age.rw.group`.

- time.model:

  Model for the main temporal trend, can be rw1, rw2, ar1, or NULL (for
  spatial-only smoothing). Default to be rw2. For ar1 main effect, a
  linear slope is also added with time scaled to be between -0.5 to 0.5,
  i.e., the slope coefficient represents the total change between the
  first year and the last year in the projection period on the logit
  scale.

- st.time.model:

  Temporal component model for the interaction term, can be rw1, rw2, or
  ar1. Default to be the same as time.model unless specified otherwise.
  The default does not include region-specific random slopes. They can
  be added to the interaction term by specifying `pc.st.slope.u` and
  `pc.st.slope.alpha`.

- Amat:

  Adjacency matrix for the regions

- bias.adj:

  the ratio of unadjusted mortality rates or age-group-specific hazards
  to the true rates or hazards. It needs to be a data frame that can be
  merged to thee outcome, i.e., with the same column names for time
  periods (for national adjustment), or time periods and region (for
  subnational adjustment). The column specifying the adjustment ratio
  should be named "ratio".

- bias.adj.by:

  vector of the column names specifying how to merge the bias adjustment
  to the count data. For example, if bias adjustment factor is provided
  in bias.adj for each region and time, then bias.adj.by should be
  `c("region", "time")`.

- formula:

  INLA formula. See vignette for example of using customized formula.

- year.label:

  string vector of year names

- year_label:

  **\[deprecated\]** replaced by year.label

- type.st:

  type for space-time interaction

- survey.effect:

  logical indicator whether to include a survey fixed effect. If this is
  set to TRUE, there needs to be a column named 'survey' in the input
  data frame. In prediction, this effect term will be set to 0.

- linear.trend:

  logical indicator whether a linear trend is added to the temporal main
  effect. If the temporal main effect is RW2, then it will be forced to
  FALSE. Default is TRUE.

- common.trend:

  logical indicator whether all age groups and strata share the same
  linear trend in the temporal main effect.

- strata.time.effect:

  logical indicator whether to include strata specific temporal trends.

- hyper:

  Deprecated. which hyperpriors to use. Only supports PC prior ("pc").

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

- pc.st.slope.u:

  hyperparameter U for the PC prior on precisions for the area-level
  random slope. If both pc.st.slope.u and pc.st.slope.alpha are not NA,
  an area-level random slope with iid prior will be added to the model.
  The parameterization of the random slope is so that Prob(\|beta\| \>
  pc.st.slope.u) = pc.st.slope.alpha, where time covariate is rescaled
  to be -0.5 to 0.5, so that the random slope can be interpreted as the
  total deviation from the main trend from the first year to the last
  year to be projected, on the logit scale.

- pc.st.slope.alpha:

  hyperparameter alpha for the PC prior on precisions for the area-level
  random slope. See above for the parameterization.

- overdisp.mean:

  hyperparameter for the betabinomial likelihood. Mean of the
  over-dispersion parameter on the logit scale.

- overdisp.prec:

  hyperparameter for the betabinomial likelihood. Precision of the
  over-dispersion parameter on the logit scale.

- options:

  list of options to be passed to control.compute() in the inla()
  function.

- control.inla:

  list of options to be passed to control.inla() in the inla() function.
  Default to the "adaptive" integration strategy.

- control.fixed:

  list of options to be passed to control.fixed() in the inla()
  function.

- verbose:

  logical indicator to print out detailed inla() intermediate steps.

- geo:

  Deprecated. Spatial polygon file, legacy parameter from previous
  versions of the package.

- rw:

  Deprecated. Take values 0, 1 or 2, indicating the order of random
  walk. If rw = 0, the autoregressive process is used instead of the
  random walk in the main trend. See the description of the argument ar
  for details.

- ar:

  Deprecated. Order of the autoregressive component. If ar is specified
  to be positive integer, the random walk components will be replaced by
  AR(p) terms in the interaction part. The main temporal trend remains
  to be random walk of order rw unless rw = 0.

- st.rw:

  Deprecated. Take values 1 or 2, indicating the order of random walk
  for the interaction term. If not specified, it will take the same
  order as the argument rw in the main effect. Notice that this argument
  is only used if ar is set to 0.

- age.rw.group:

  Deprecated. Legacy parameter replaced by `age.time.group`.

- ...:

  arguments to be passed to the inla() function call.

## Value

INLA model fit using the provided formula, country summary data, and
geographic data

## See also

[`getDirect`](https://richardli.github.io/SUMMER/reference/getDirect.md)

## Author

Zehang Richard Li

## Examples

``` r
if (FALSE) { # \dontrun{
library(dplyr)
data(DemoData)
# Create dataset of counts
counts.all <- NULL
for(i in 1:length(DemoData)){
  counts <- getCounts(DemoData[[i]][, c("clustid", "time", "age", "died",
                                        "region", "strata")],
            variables = 'died', by = c("age", "clustid", "region", 
                                         "time", "strata"))
  counts <- counts %>% mutate(cluster = clustid, years = time, Y=died)
  counts$strata <- gsub(".*\\.","",counts$strata)
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
summary(fit)
est <- getSmoothed(fit, nsim = 1000)
plot(est$stratified, plot.CI=TRUE) + ggplot2::facet_wrap(~strata) 

# fit cluster-level space-time model with covariate
# notice without projected covariates, we use periods up to 10-14 only
# construct a random covariate matrix for illustration
periods <- levels(DemoData[[1]]$time)
X <- expand.grid(years = periods, 
       region = unique(counts.all$region))
X$X1 <- rnorm(dim(X)[1])
X$X2 <- rnorm(dim(X)[1])
fit.covariate <- smoothCluster(data = counts.all, 
   X = X,
     Amat = DemoMap$Amat, 
     time.model = "rw2", 
     st.time.model = "rw1",
     strata.time.effect =  TRUE, 
     survey.effect = TRUE,
     family = "betabinomial",
     year.label = c(periods))
est <- getSmoothed(fit.covariate, nsim = 1000)

# fit cluster-level model for one time point only
# i.e., space-only model
fit.sp <- smoothCluster(data = subset(counts.all, time == "10-14"), 
     Amat = DemoMap$Amat, 
     time.model = NULL, 
     survey.effect = TRUE,
     family = "betabinomial")
summary(fit.sp)
est <- getSmoothed(fit.sp, nsim = 1000)
plot(est$stratified, plot.CI = TRUE) + ggplot2::facet_wrap(~strata) 

# fit cluster-level model for one time point and covariate
# construct a random covariate matrix for illustration
X <- data.frame(region = unique(counts.all$region),
      X1 = c(1, 2, 2, 1), 
      X2 = c(1, 1, 1, 2))
fit.sp.covariate <- smoothCluster(data = subset(counts.all, time == "10-14"), 
     X = X, 
     Amat = DemoMap$Amat, 
     time.model = NULL, 
     survey.effect = TRUE,
     family = "betabinomial")
summary(fit.sp.covariate)
est <- getSmoothed(fit.sp.covariate, nsim = 1000)
} # }
```
