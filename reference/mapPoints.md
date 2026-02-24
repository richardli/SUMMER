# Map GPS points to polygon regions

Map GPS points to polygon regions

## Usage

``` r
mapPoints(data, geo, long, lat, names)
```

## Arguments

- data:

  point data with two columns of GPS locations.

- geo:

  SpatialPolygonsDataFrame of the map

- long:

  column name for longitudinal coordinate in the data

- lat:

  column name for latitude coordinate in the data

- names:

  character vector of region ids to be added to the neighbours list

## Value

Spatial djacency matrix.

## Author

Zehang Richard Li

## Examples

``` r
data(DemoMap) 
dat <- data.frame(ID = c(1,2,3), lon = c(32.2, 33.7, 33), lat = c(0.1, 0.9, 2.8))
dat2 <- mapPoints(dat, DemoMap$geo, long = "lon", lat = "lat", names = "REGNAME")
dat2
#>   ID  lon lat  REGNAME
#> 1  1 32.2 0.1  central
#> 2  2 33.7 0.9  eastern
#> 3  3 33.0 2.8 northern
 
```
