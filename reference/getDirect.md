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
  CI = 0.95
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

## Value

a matrix of period-region summary of the Horvitz-Thompson direct
estimates by region and time period specified in the argument, the
standard errors using delta method for a single survey, the 95%
confidence interval, and the logit of the estimates.

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
