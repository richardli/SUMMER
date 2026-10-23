# Obtain the Horvitz-Thompson direct estimates and standard errors using delta method for a single survey.

Obtain the Horvitz-Thompson direct estimates and standard errors using
delta method for a single survey.

## Usage

``` r
getDirect(
  births,
  years,
  regionVar = "region",
  timeVar = "time",
  clusterVar = "~v001+v002",
  ageVar = "age",
  weightsVar = "v005",
  Ntrials = NULL,
  geo.recode = NULL,
  national.only = FALSE,
  CI = 0.95,
  var.fix = FALSE,
  refStrataVar = NULL
)
```

## Arguments

- births:

  A matrix child-month data from
  [`getBirths`](https://richardli.github.io/SUMMER/reference/getBirths.md)

- years:

  String vector of the year intervals used

- regionVar:

  Variable name for region in the input births data.

- timeVar:

  Variable name for the time period indicator in the input births data.

- clusterVar:

  Variable name for cluster, typically '~v001 + v002'

- ageVar:

  Variable name for age group. This variable need to be in the form of
  "a-b" where a and b are both ages in months. For example, "1-11" means
  age between 1 and 11 months, including both end points. An exception
  is age less than one month can be represented by "0" or "0-0".

- weightsVar:

  Variable name for sampling weights, typically 'v005'

- Ntrials:

  Variable for the total number of person-months if the input data
  (births) is in the compact form.

- geo.recode:

  The recode matrix to be used if region name is not consistent across
  different surveys. See `changeRegion`.

- national.only:

  Logical indicator to obtain only the national estimates

- CI:

  the desired confidence interval to calculate

- var.fix:

  Logical; if `TRUE`, we attempt a target-specific synthetic-PSU
  correction when Rule 2 detects no between-PSU hazard contrast. The
  default preserves the usual estimates and attributes.

- refStrataVar:

  Optional name of a reference-group column. With `NULL`, we pool all
  real observations across areas and sampling strata separately within
  each period and age band. An explicit column must have one nonmissing
  value per design stratum within each period. We never infer this
  column from `births$strata`.

## Value

A data frame of period-region direct estimates, confidence limits, logit
estimates, logit variances and precisions, and a survey placeholder.
When `var.fix = TRUE`, we also attach `fixed_areas`, a character vector
of triggered regions (including `"All"` when triggered), and
`variance.fix`, a data frame with one row per triggered region-period.
These attributes include unresolved attempts; `fixed_areas` does not
imply successful correction.

## Details

With `var.fix = TRUE`, we evaluate Rule 2 separately for each area and
period. A single real PSU triggers correction. Otherwise, correction
requires every exposed PSU's weighted hazard to equal the area's
weighted hazard in every exposed age band, within \\10^{-12}\\. PSU
identities include the sampling stratum, so repeated PSU labels in
different strata remain distinct.

We construct references before converting compact death counts to
proportions. Within each reference group and period, \\K\\ counts all
real PSUs with retained exposure, including PSUs without exposure in a
particular age band. We average each PSU's exposure-weighted mean
sampling weight to obtain \\w^\*\\. For each age band, synthetic
exposure is \\E^\*\_a = \sum_i w_i E\_{ia}/(K w^\*)\\, and synthetic
deaths are \\D^\*\_a = \sum_i w_i D\_{ia}/(K w^\*)\\. Fractional counts
are supported. Observations with missing area labels still contribute to
references and national estimates. Zero-exposure rows do not contribute
to the correction.

We add one synthetic PSU in each design stratum present in a triggered
target, using that stratum's reference group. We augment a separate
full-survey design for each target before subsetting it and applying the
existing estimator. We retain the original design for untriggered
targets and treat the national estimate as its own target. Thus, another
area's synthetic PSUs cannot change an untriggered estimate or its
variance. The correction can change both the point estimate and
variance. Simple named multistage cluster formulas, such as
`~psu + household`, are supported in the correction path.

We retain the zero-death guard after augmentation. Targets without data
remain missing, and degenerate references need not produce positive
variance. We do not impose a variance floor or pool reference periods.
The audit has character columns `region`, `years`, `status`, `reason`,
`trigger`, and `reference.groups`; integer counts `real.psus`,
`phantom.psus`, and `phantom.rows`; and numeric totals
`phantom.weighted.deaths` and `phantom.weighted.exposure`. We label an
attempt `"fixed"` only when the augmented target no longer satisfies
Rule 2 and has a finite logit estimate and finite positive variance.
Otherwise, we label it `"unresolved"` and record the reason. Empty
trigger and audit attributes retain the same types and columns.

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

[`getDirectList`](https://richardli.github.io/SUMMER/reference/getDirectList.md)

## Author

Zehang Richard Li, Bryan Martin, Laina Mercer

## Examples

``` r
if (FALSE) { # \dontrun{
data(DemoData)
years <- c("85-89", "90-94", "95-99", "00-04", "05-09", "10-14")
mean <- getDirect(births = DemoData[[1]],  years = years, 
regionVar = "region", timeVar = "time", clusterVar = "~clustid+id", 
ageVar = "age", weightsVar = "weights", geo.recode = NULL)
} # }
```
