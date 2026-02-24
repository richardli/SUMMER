# Small area estimation via basic area level model

Generates small area estimates by smoothing direct estimates using an
area level model. For continuous outcome, standardizing the outcome
variable so that it has unit standard deviation is recommended when
using the default priors.

## Usage

``` r
smoothArea(
  formula,
  domain,
  design = NULL,
  adj.mat = NULL,
  X.domain = NULL,
  direct.est = NULL,
  domain.size = NULL,
  transform = c("identity", "logit", "log"),
  pc.u = 1,
  pc.alpha = 0.01,
  pc.u.phi = 0.5,
  pc.alpha.phi = 2/3,
  level = 0.95,
  n.sample = 250,
  var.tol = 1e-10,
  return.samples = F
)
```

## Arguments

- formula:

  An object of class 'formula' describing the model to be fitted. If
  direct.est is specified, the right hand side of the formula is not
  necessary.

- domain:

  One-sided formula specifying factors containing domain labels

- design:

  An object of class "svydesign" containing the data for the model

- adj.mat:

  Adjacency matrix with rownames matching the domain labels. If set to
  NULL, the IID spatial effect will be used.

- X.domain:

  Data frame of areal covariates. One of the column names needs to match
  the name of the domain variable, in order to be linked to the data
  input. Currently only supporting time-invariant covariates.

- direct.est:

  Data frame of direct estimates, with first column containing the
  domain variable, second column containing direct estimate, and third
  column containing the variance of direct estimate.

- domain.size:

  Data frame of domain sizes. One of the column names needs to match the
  name of the domain variable, in order to be linked to the data input
  and there must be a column names 'size' containing domain sizes.

- transform:

  Optional transformation applied to the direct estimates before fitting
  area level model. The default option is no transformation, but logit
  and log are implemented.

- pc.u:

  Hyperparameter U for the PC prior on precisions. See the INLA
  documentation for more details on the parameterization.

- pc.alpha:

  Hyperparameter alpha for the PC prior on precisions.

- pc.u.phi:

  Hyperparameter U for the PC prior on the mixture probability phi in
  BYM2 model.

- pc.alpha.phi:

  Hyperparameter alpha for the PC prior on the mixture probability phi
  in BYM2 model.

- level:

  The specified level for the posterior credible intervals

- n.sample:

  Number of draws from posterior used to compute summaries

- var.tol:

  Tolerance parameter; if variance of an area's direct estimator is
  below this value, that direct estimator is dropped from model

- return.samples:

  If TRUE, return matrix of posterior samples of area level quantities

## Value

A svysae object

## Examples

``` r
if (FALSE) { # \dontrun{
data(DemoData2)
data(DemoMap2)
library(survey)
des0 <- svydesign(ids = ~clustid+id, strata = ~strata,
                  weights = ~weights, data = DemoData2, nest = TRUE)
Xmat <- aggregate(age~region, data = DemoData2, FUN = mean)

# EXAMPLE 1: Continuous response model
cts.res <- smoothArea(tobacco.use ~ 1,
                      domain = ~region,
                      design = des0,
                      adj.mat = DemoMap2$Amat, 
                      pc.u = 1,
                      pc.alpha = 0.01,
                      pc.u.phi = 0.5,
                      pc.alpha.phi = 2/3)

# EXAMPLE 2: Including area level covariates
cts.cov.res <- smoothArea(tobacco.use ~ age, 
                          domain = ~region,
                          design = des0,
                          adj.mat = DemoMap2$Amat, 
                          X.domain = Xmat,
                          pc.u = 1,
                          pc.alpha = 0.01,
                          pc.u.phi = 0.5,
                          pc.alpha.phi = 2/3)

# EXAMPLE 3: Binary response model
bin.res <- smoothArea(tobacco.use ~ 1, 
                      domain = ~region,
                      design = des0,
                      adj.mat = DemoMap2$Amat, 
                      transform = "logit",
                      pc.u = 1,
                      pc.alpha = 0.01,
                      pc.u.phi = 0.5,
                      pc.alpha.phi = 2/3)

# EXAMPLE 4: Including area level covariates in binary response model
bin.cov.res <- smoothArea(tobacco.use ~ age, 
                          domain = ~region,
                          design = des0,
                          adj.mat = DemoMap2$Amat, 
                          transform = "logit",
                          X.domain = Xmat,
                          pc.u = 1,
                          pc.alpha = 0.01,
                          pc.u.phi = 0.5,
                          pc.alpha.phi = 2/3)
} # }
```
