# Determines which administrative areas contain the given points

**\[experimental\]**

## Usage

``` r
getAreaName(
  pts,
  shapefile,
  areaNameVar = "NAME_1",
  delta = 0.05,
  mean.neighbor = 50,
  max.bytes = 3 * 2^30
)
```

## Arguments

- pts:

  2 column matrix of lon/lat coordinates

- shapefile:

  A SpatialPolygonsDataFrame object

- areaNameVar:

  The column name in `slot(shapefile, "data")` corresponding to the area
  level of interest

- delta:

  Argument passed to fields::fields.rdist.near in fields package

- mean.neighbor:

  Argument passed to fields::fields.rdist.near in fields package

- max.bytes:

  Maximum allowed memory in bytes (default is 3Gb). Determines whether
  to call fields::fields.rdist.near which saves memory but requires
  delta and mean.neighbor inputs to be specified for
  fields::fields.rdist.near

## Value

A list of area IDs, area names, whether or not points are in multiple
areas, and whether or not points are in no areas and assigned to the
nearest one.

## Details

For any points not in an area, they are assigned the nearest area using
fields::fields.rdist.near or fields::rdist depending on the number of
points and the maximum memory in bytes with a warning.

delta and mean.neighbor arguments only used when some points are not in
areas, perhaps due to inconsistencies in shapefiles.

## See also

[`projKenya`](https://richardli.github.io/SUMMER/reference/projKenya.md),
[`fields.rdist.near`](https://rdrr.io/pkg/fields/man/rdist.html)

## Author

John Paige

## Examples

``` r
if (FALSE) { # \dontrun{
# download Kenya GADM shapefiles from SUMMERdata github repository
githubURL <- "https://github.com/paigejo/SUMMERdata/blob/main/data/kenyaMaps.rda?raw=true"
download.file(githubURL,"kenyaMaps.rda")

# load it in
load("kenyaMaps.rda")

# use the shapefile data to see what Admin1 and 2 areas the 
# points (0, 37) and (0.5, 38) are in
# (these are longitude/latitude coordinates)
pts = cbind(c(37, 38), c(0, .5))
head(slot(adm1, "data"))
admin1Areas = getAreaName(pts, adm1, "NAME_1")
admin2Areas = getAreaName(pts, adm2, "NAME_2")
} # }
```
