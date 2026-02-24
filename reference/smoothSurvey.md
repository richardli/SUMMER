# Fit space-time smoothing models for a binary outcome from complex surveys.

This function calculates the direct estimates by region and fit a simple
spatial smoothing model to the direct estimates adjusting for survey
design. Normal or binary variables are currently supported. For binary
variables, the logit transformation is performed on the direct estimates
of probabilities, and a Gaussian additive model is fitted on the logit
scale using INLA.

## Usage

``` r
smoothSurvey(
  data,
  geo = NULL,
  Amat = NULL,
  region.list = NULL,
  X = NULL,
  X.unit = NULL,
  responseType = deprecated(),
  response.type = c("binary", "gaussian")[1],
  responseVar,
  strataVar = "strata",
  weightVar = "weights",
  regionVar = "region",
  clusterVar = "~v001+v002",
  pc.u = 1,
  pc.alpha = 0.01,
  pc.u.phi = 0.5,
  pc.alpha.phi = 2/3,
  CI = 0.95,
  formula = NULL,
  timeVar = NULL,
  time.model = c("rw1", "rw2")[1],
  include_time_unstruct = deprecated(),
  include.time.unstruct = FALSE,
  type.st = 1,
  direct.est = NULL,
  direct.est.var = NULL,
  is.unit.level = FALSE,
  is.agg = FALSE,
  strataVar.within = NULL,
  totalVar = NULL,
  weight.strata = NULL,
  nsim = 1000,
  save.draws = FALSE,
  smooth = TRUE,
  ...
)
```

## Arguments

- data:

  The input data frame. The input data with column of the response
  variable (`responseVar`), region ID (`regionVar`), stratification
  within region (`strataVar`), and cluster ID (`clusterVar`).

  - For area-level model, the data frame consist of survey observations
    and corresponding survey weights (`weightVar`).

  - For unit-level model and `is.agg = FALSE`, the data frame should
    consist of aggregated counts by clusters (for binary responses), or
    any cluster-level response (for continuous response). For binary
    response (`response.type = 'binary'`), the beta-binomial model will
    be fitted for cluster-level counts. For continuous response
    (`response.type = 'gaussian'`), a Gaussian smoothing model will be
    fitted on the cluster-level response.

  - For unit-level model and `is.agg = TRUE`, the data frame should be
    the same as in the area-level model. For binary response
    (`response.type = 'binary'`), the beta-binomial model will be fitted
    for cluster-level counts aggregated internally. For continuous
    response (`response.type = 'gaussian'`), the nested error model will
    be fitted on unit-level response.

- geo:

  Deprecated argument from early versions.

- Amat:

  Adjacency matrix for the regions. If set to NULL, the IID spatial
  effect will be used.

- region.list:

  a vector of region names. Only used when IID model is used and the
  adjacency matrix not specified. This allows the output to include
  regions with no sample in the data. When the spatial adjacency matrix
  is specified, the column names of the adjacency matrix will be used to
  determine region.list. If set to NULL, all regions in the data are
  used.

- X:

  Areal covariates data frame. One of the column name needs to match the
  `regionVar` specified in the function call, in order to be linked to
  the data input. Currently only supporting time-invariant region-level
  covariates.

- X.unit:

  Column names of unit-level covariates. When `X.unit` is specified, a
  nested error model will be fitted with unit-level IID noise, and
  area-level predictions are produced by plugging in the covariate
  specified in the `X` argument. When `X` is not specified, the
  empirical mean of each covariate will be used. This is only
  implemented for continuous response with the Gaussian likelihood model
  and unit-level model.

- responseType:

  **\[deprecated\]** The argument has been renamed into `response.type`.

- response.type:

  Type of the response variable, currently supports 'binary' (default
  with logit link function) or 'gaussian'.

- responseVar:

  the response variable

- strataVar:

  the strata variable used in the area-level model.

- weightVar:

  the weights variable

- regionVar:

  Variable name for region.

- clusterVar:

  Variable name for cluster. For area-level model, this should be a
  formula for cluster in survey design object, e.g., '~clusterID +
  householdID'. For unit-level model, this should be the variable name
  for cluster unit.

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

- CI:

  the desired posterior credible interval to calculate

- formula:

  a string of user-specified random effects model to be used in the INLA
  call

- timeVar:

  The variable indicating time period. If set to NULL then the temporal
  model and space-time interaction model are ignored.

- time.model:

  the model for temporal trends and interactions. It can be either "rw1"
  or "rw2".

- include_time_unstruct:

  **\[deprecated\]** The argument has been renamed into
  `include.time.unstruct`.

- include.time.unstruct:

  Indicator whether to include the temporal unstructured effects (i.e.,
  shocks) in the smoothed estimates from cluster-level model. The
  argument only applies to the unit-level models. Default is FALSE which
  excludes all unstructured temporal components. If set to TRUE all the
  unstructured temporal random effects will be included.

- type.st:

  can take values 0 (no interaction), or 1 to 4, corresponding to the
  type I to IV space-time interaction.

- direct.est:

  data frame of direct estimates, with column names of response and
  region specified by `responseVar`, `regionVar`, and `timeVar`. When
  `direct.est` is specified, it overwrites the `data` input.

- direct.est.var:

  the column name corresponding to the variance of direct estimates.

- is.unit.level:

  logical indicator of whether unit-level model is fitted instead of
  area-level model.

- is.agg:

  logical indicator of whether the input is at the aggregated counts by
  cluster. Only used for unit-level model and binary response variable.

- strataVar.within:

  the variable specifying within region stratification variable. This is
  only used for the unit-level model.

- totalVar:

  the variable specifying total observations in `counts`. This is only
  used for the unit-level model when `counts` is specified.

- weight.strata:

  a data frame with one column corresponding to `regionVar`, and columns
  specifying proportion of each strata for each region. This argument
  specifies the weights for strata-specific estimates. This is only used
  for the unit-level model.

- nsim:

  number of posterior draws to take. This is only used for the
  unit-level model when `weight.strata` is provided.

- save.draws:

  logical indicator of whether to save the full posterior draws.

- smooth:

  logical indicator of whether to perform smoothing. If set to FALSE, a
  data frame of direct estimate is returned. Only used when
  `is.unit.level` is FALSE.

- ...:

  additional arguments passed to `svydesign` function.

## Value

- HT:

  Direct estimates

- smooth:

  Smoothed direct estimates

- fit:

  a fitted INLA object

- CI:

  input argument

- Amat:

  input argument

- response.type:

  input argument

- formula:

  INLA formula

## Details

The function `smoothSurvey` replaces the previous function name
`fitGeneric` (before version 1.0.0).

## See also

[`getDirectList`](https://richardli.github.io/SUMMER/reference/getDirectList.md),
[`smoothDirect`](https://richardli.github.io/SUMMER/reference/smoothDirect.md)

## Author

Zehang Richard Li

## Examples

``` r
if (FALSE) { # \dontrun{
##
## 1. Area-level model with binary response
##

data(DemoData2)
data(DemoMap2)
fit0 <- smoothSurvey(data=DemoData2,  
Amat=DemoMap2$Amat, response.type="binary", 
responseVar="tobacco.use", strataVar="strata", 
weightVar="weights", regionVar="region", 
clusterVar = "~clustid+id", CI = 0.95)
summary(fit0)

# if only direct estimates without smoothing is of interest
fit0.dir <- smoothSurvey(data=DemoData2,  
Amat=DemoMap2$Amat, response.type="binary", 
responseVar="tobacco.use", strataVar="strata", 
weightVar="weights", regionVar="region", 
clusterVar = "~clustid+id", CI = 0.95, smooth = FALSE)

# posterior draws can be returned with save.draws = TRUE
fit0.draws <- smoothSurvey(data=DemoData2,  
Amat=DemoMap2$Amat, response.type="binary", 
responseVar="tobacco.use", strataVar="strata", 
weightVar="weights", regionVar="region", 
clusterVar = "~clustid+id", CI = 0.95, save.draws = TRUE)
# notice the posterior draws are on the latent scale
head(fit0.draws$draws.est[, 1:10]) 

# Example with region-level covariates
 Xmat <- aggregate(age~region, data = DemoData2, 
            FUN = function(x) mean(x))
 fit1 <- smoothSurvey(data=DemoData2,  
  Amat=DemoMap2$Amat, response.type="binary", 
  X = Xmat,
  responseVar="tobacco.use", strataVar="strata", 
  weightVar="weights", regionVar="region", 
  clusterVar = "~clustid+id", CI = 0.95)

# Example with using only direct estimates as input instead of the full data
direct <- fit0$direct[, c("region", "direct.est", "direct.var")]
fit2 <- smoothSurvey(data=NULL, direct.est = direct, 
                    Amat=DemoMap2$Amat, regionVar="region",
                    responseVar="direct.est", direct.est.var = "direct.var", 
                    response.type = "binary")
# Check it is the same as fit0
plot(fit2$smooth$mean, fit0$smooth$mean)

# Example with using only direct estimates as input, 
#   and after transformation into a Gaussian smoothing model
# Notice: the output are on the same scale as the input 
#   and in this case, the logit estimates.    
direct.logit <- fit0$direct[, c("region", "direct.logit.est", "direct.logit.var")]
fit3 <- smoothSurvey(data=NULL, direct.est = direct.logit, 
               Amat=DemoMap2$Amat, regionVar="region",
               responseVar="direct.logit.est", direct.est.var = "direct.logit.var",
               response.type = "gaussian")
# Check it is the same as fit0
plot(fit3$smooth$mean, fit0$smooth$logit.mean)

# Example with non-spatial smoothing using IID random effects
fit4 <- smoothSurvey(data=DemoData2, response.type="binary", 
       responseVar="tobacco.use", strataVar="strata", 
       weightVar="weights", regionVar="region", 
       clusterVar = "~clustid+id", CI = 0.95)

# Example with missing regions in the raw input
DemoData2.sub <- subset(DemoData2, region != "central")
fit.without.central <- smoothSurvey(data=DemoData2.sub,  
                         Amat=NULL, response.type="binary", 
                         responseVar="tobacco.use", strataVar="strata", 
                         weightVar="weights", regionVar="region", 
                         clusterVar = "~clustid+id", CI = 0.95)
fit.without.central$direct
fit.without.central$smooth

fit.with.central <- smoothSurvey(data=DemoData2.sub,  
                         Amat=NULL, region.list = unique(DemoData2$region),
                         response.type="binary", 
                         responseVar="tobacco.use", strataVar="strata", 
                         weightVar="weights", regionVar="region", 
                         clusterVar = "~clustid+id", CI = 0.95)
fit.with.central$direct
fit.with.central$smooth

# Using the formula argument, further customizations can be added to the 
#  model fitted. For example, we can fit the Fay-Harriot model with 
#  IID effect instead of the BYM2 random effect as follows.
#  The "region.struct" and "hyperpc1" are picked to match internal object 
#  names. Other object names can be inspected from the source of smoothSurvey.
fit5 <- smoothSurvey(data=DemoData2,  
       Amat=DemoMap2$Amat, response.type="binary", 
       formula = "f(region.struct, model = 'iid', hyper = hyperpc1)",
       pc.u = 1, pc.alpha = 0.01,
       responseVar="tobacco.use", strataVar="strata", 
       weightVar="weights", regionVar="region", 
       clusterVar = "~clustid+id", CI = 0.95)
# Check it is the same as fit4, notice the region order may be different
regions <- fit5$smooth$region
plot(fit4$smooth[match(regions, fit4$smooth$region),]$logit.mean, fit5$smooth$logit.mean)

##
## 2. Unit-level model with binary response  
##

# For unit-level models, we need to create stratification variable within regions
data <- DemoData2
data$urbanicity <- "rural"
data$urbanicity[grep("urban", data$strata)] <- "urban"

# Beta-binomial likelihood is used in this model
fit6 <- smoothSurvey(data=data, 
  Amat=DemoMap2$Amat, response.type="binary", 
  X = Xmat, is.unit.level = TRUE,
  responseVar="tobacco.use", strataVar.within = "urbanicity", 
  regionVar="region", clusterVar = "clustid", CI = 0.95)

# We may use aggregated PSU-level counts as input as well
#    in the case of modeling a binary outcome 
data.agg <- aggregate(tobacco.use~region + urbanicity + clustid, 
                      data = data, FUN = sum)
data.agg.total <- aggregate(tobacco.use~region + urbanicity + clustid, 
                      data = data, FUN = length)
colnames(data.agg.total)[4] <- "total"
data.agg <- merge(data.agg, data.agg.total)
head(data.agg)

fit7 <- smoothSurvey(data=data.agg, 
  Amat=DemoMap2$Amat, response.type="binary", 
  X = Xmat, is.unit.level = TRUE, is.agg = TRUE,
  responseVar = "tobacco.use", strataVar.within = "urbanicity", 
  totalVar = "total", regionVar="region", clusterVar = "clustid", CI = 0.95)

# Check it is the same as fit6
plot(fit6$smooth$mean, fit7$smooth$mean)  

##
## 3. Area-level model with continuous response
##

# The smoothing model is the same as area-level model with binary response
#  the continuous direct estimates are smoothed instead of 
#  their logit-transformed versions for binary response.
fit8 <- smoothSurvey(data=DemoData2, Amat=DemoMap2$Amat, 
       response.type="gaussian", responseVar="age", strataVar="strata", 
       weightVar="weights", regionVar="region", 
       pc.u.phi = 0.5, pc.alpha.phi = 0.5,
       clusterVar = "~clustid+id", CI = 0.95)

##
## 4. Unit-level model with continuous response  
##    (or nested error models)

# The unit-level model assumes for each of the i-th unit,
#    Y_{i} ~ intercept + region_effect + IID_i
#    where IID_i is the error term specific to i-th unit

# When more than one level of cluster sampling is carried out, 
#   they are ignored here. Only the input unit is considered.
#   So here we do not need to specify clusterVar any more. 
fit9 <- smoothSurvey(data= data, 
  Amat=DemoMap2$Amat, response.type="gaussian", 
  is.unit.level = TRUE, responseVar="age", strataVar.within = NULL,
  regionVar="region", clusterVar = NULL, CI = 0.95)

# To compare, we may also model PSU-level responses. As an illustration, 
data.median <- aggregate(age~region + urbanicity + clustid, 
                      data = data, FUN = median)

fit10 <- smoothSurvey(data= data.median, 
  Amat=DemoMap2$Amat, response.type="gaussian", 
  is.unit.level = TRUE, responseVar="age", strataVar.within = NULL,
  regionVar="region", clusterVar = "clustid", CI = 0.95)


# To further incorporate within-area stratification

fit11 <- smoothSurvey(data = data, 
  Amat = DemoMap2$Amat, response.type = "gaussian", 
  is.unit.level = TRUE, responseVar="age", strataVar.within = "urbanicity",
  regionVar = "region", clusterVar = NULL, CI = 0.95)  

# Notice the usual output is now stratified within each region
# The aggregated estimates require strata proportions for each region
# For illustration, we set strata population proportions below
prop <- data.frame(region = unique(data$region), 
                            urban = 0.3, 
                            rural = 0.7)
fit12 <- smoothSurvey(data=data, 
  Amat=DemoMap2$Amat, response.type="gaussian", 
  is.unit.level = TRUE, responseVar="age", strataVar.within = "urbanicity",
  regionVar="region", clusterVar = NULL, CI = 0.95,
  weight.strata = prop)  

# aggregated outcome
head(fit12$smooth.overall)

# Compare aggregated outcome with direct aggregating posterior means. 
# There could be small differences if only 1000 posterior draws are taken.
est.urb <- subset(fit11$smooth, strata == "urban")
est.rural <- subset(fit11$smooth, strata == "rural")
est.mean.post <- est.urb$mean * 0.3 + est.rural$mean * 0.7
plot(fit12$smooth.overall$mean, est.mean.post)


##
## 6. Unit-level model with continuous response and unit-level covariate 
## 

# For area-level prediction, area-level covariate mean needs to be  
#   specified in X argument. And unit-level covariate names are specified
#   in X.unit argument.

set.seed(1)
sim <- data.frame(region = rep(c(1, 2, 3, 4), 1000),
                   X1 = rnorm(4000), X2 = rnorm(4000))
Xmean <- aggregate(.~region, data = sim, FUN = sum)
sim$Y <- rnorm(4000, mean = sim$X1 + 0.3 * sim$X2 + sim$region)
samp <- sim[sample(1:4000, 20), ]
fit.sim <- smoothSurvey(data=samp , 
                  X.unit = c("X1", "X2"),
                  X = Xmean, Amat=NULL, response.type="gaussian", 
                  is.unit.level = TRUE, responseVar="Y", regionVar = "region",  
                  pc.u = 1, pc.alpha = 0.01, CI = 0.95) 

} # }
```
