# Extract adjacency matrix from the map

Extract adjacency matrix from the map

## Usage

``` r
getAmat(geo, names)
```

## Arguments

- geo:

  SpatialPolygonsDataFrame of the map

- names:

  character vector of region ids to be added to the neighbours list

## Value

Spatial djacency matrix.

## Author

Zehang Richard Li

## Examples

``` r
if (FALSE) { # \dontrun{
data(DemoMap) 
mat <- getAmat(geo = DemoMap$geo, names = DemoMap$geo$REGNAME)
mat
DemoMap$Amat
} # } 
```
