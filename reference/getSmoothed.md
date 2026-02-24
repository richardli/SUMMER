# Extract smoothed estimates.

Extract smoothed estimates.

## Usage

``` r
getSmoothed(
  fitted,
  inla_mod = deprecated(),
  nsim = 1000,
  weight.strata = NULL,
  weight.frame = NULL,
  verbose = FALSE,
  mc = 0,
  include.time.unstruct = FALSE,
  include_time_unstruct = deprecated(),
  CI = 0.95,
  draws = NULL,
  save.draws = FALSE,
  include.subnational = TRUE,
  include_subnational = deprecated(),
  only.age = NULL,
  joint = FALSE,
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

  **\[deprecated\]** replaced by `fitted`

- nsim:

  number of simulations, only applicable for the cluster-level model or
  direct model when `joint = TRUE` is specified. The smooth direct model
  will draw 1e5 samples from the marginal distribution when
  `joint = FALSE` since the computation is faster.

- weight.strata:

  a data frame with two columns specifying time and region, followed by
  columns specifying proportion of each strata for each region. This
  argument specifies the weights for strata-specific estimates on the
  probability scale.

- weight.frame:

  a data frame with three columns, years, region, and the weight of each
  frame for the corresponding time period and region. This argument
  specifies the weights for frame-specific estimates on the logit scale.
  Notice this is different from weight.strata argument.

- verbose:

  logical indicator whether to print progress messages from
  inla.posterior.sample.

- mc:

  number of monte carlo draws to approximate the marginal
  prevalence/hazards for binomial model. If mc = 0, analytical
  approximation is used. The analytical approximation is invalid for
  hazard modeling with more than one age groups.

- include.time.unstruct:

  Indicator whether to include the temporal unstructured effects (i.e.,
  shocks) in the smoothed estimates from cluster-level model. The
  argument only applies to the cluster-level models (from
  [`smoothCluster`](https://richardli.github.io/SUMMER/reference/smoothCluster.md)).
  Default is FALSE which excludes all unstructured temporal components.
  If set to TRUE all the unstructured temporal random effects will be
  included. Alternatively, if this is specified as a vector of subset of
  year labels (as in the year.label argument), only the unstructured
  terms in the corresponding time periods will be added to the
  prediction.

- include_time_unstruct:

  **\[deprecated\]** replaced by `include.time.unstruct`

- CI:

  Desired level of credible intervals

- draws:

  Posterior samples drawn from the fitted model. This argument allows
  the previously sampled draws (by setting save.draws to be TRUE) be
  used in new aggregation tasks.

- save.draws:

  Logical indicator whether the raw posterior draws will be saved. Saved
  draws can be used to accelerate aggregations with different weights.

- include.subnational:

  logical indicator whether to include the spatial and space-time
  interaction components in the smoothed estimates. If set to FALSE,
  only the main temporal trends are returned.

- include_subnational:

  **\[deprecated\]** replaced by `include.subnational`

- only.age:

  a vector of age groups used to compute the final estimates. Default to
  be NULL, which includes all age groups in the model. This argument can
  be used to extract mortality rates of finer age groups when fitting
  multiple age groups jointly.

- joint:

  Logical indicator whether the posterior draws should be drawn from the
  joint posterior or marginal distributions. Only releveant for the
  smooth direct model. Default from the marginal distribution (joint =
  FALSE).

- ...:

  Unused arguments, for users with fitted object from the package before
  v1.0.0, arguments including Amat, year.label, and year.range can still
  be specified manually.

## Value

A data frame or a list of data frames of S3 class SUMMERproj, which
contains the smoothed estimates.

## See also

[`plot.SUMMERproj`](https://richardli.github.io/SUMMER/reference/plot.SUMMERproj.md)

## Author

Zehang Richard Li

## Examples
