# Print method for the smoothing models from `smoothSurvey`.

This function is the print method for class `SUMMERmodel.svy`.

## Usage

``` r
# S3 method for class 'SUMMERmodel.svy'
print(x, ...)
```

## Arguments

- x:

  output from
  [`smoothSurvey`](https://richardli.github.io/SUMMER/reference/smoothSurvey.md).

- ...:

  not used

## See also

[`summary.SUMMERmodel.svy`](https://richardli.github.io/SUMMER/reference/summary.SUMMERmodel.svy.md)

## Author

Zehang Li

## Examples

``` r
if (FALSE) { # \dontrun{
data(DemoData2)
data(DemoMap2)
fit0 <- smoothSurvey(data=DemoData2,  
Amat=DemoMap2$Amat, responseType="binary", 
responseVar="tobacco.use", strataVar="strata", 
weightVar="weights", regionVar="region", 
clusterVar = "~clustid+id", CI = 0.95)
fit0
} # }
```
