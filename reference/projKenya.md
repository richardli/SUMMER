# Map projection for Kenya

**\[experimental\]**

## Usage

``` r
projKenya(lon, lat = NULL, inverse = FALSE)
```

## Arguments

- lon:

  either longitude or, if inverse == TRUE, easting in km

- lat:

  either latitude or, if inverse == TRUE, northing in km

- inverse:

  if FALSE, projects from lon/lat to easting/northing. Else from
  easting/northing to lon/lat

## Value

A 2 column matrix of easting/northing coordinates in km if inverse ==
FALSE. Otherwise, a 2 column matrix of longitude/latitude coordinates.

## Details

Projection specifically chosen for Kenya. Project from lat/lon to
northing/easting in kilometers. Uses epsg=21097 with km units. May not
work on all systems due to differences in the behavior between different
PROJ and GDAL versions.

## Author

John Paige

## Examples

``` r
eastLim = c(-110.6405, 832.4544)
northLim = c(-555.1739, 608.7130)
coordMatrixEN = cbind(eastLim, northLim)
coordMatrixLL = projKenya(coordMatrixEN, inverse=TRUE)

coordMatrixLL
#>           lon       lat
#> [1,] 33.50075 -5.002281
#> [2,] 42.00093  5.496810
# if the coordMatrixLL isn't the following, projKenya may not support 
# your installation of GDAL and/or PROJ:
#      east north
# [1,] 33.5  -5.0
# [2,] 42.0   5.5

projKenya(coordMatrixLL, inverse=FALSE)
#>           east     north
#> [1,] -110.6405 -555.1739
#> [2,]  832.4544  608.7130
# regardless of your PROJ/GDAL installations, the result of the 
# above line of could should be:
#            lon       lat
# [1,] -110.6405 -555.1739
# [2,]  832.4544  608.7130
```
