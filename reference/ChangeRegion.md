# Map region names to a common set.

Map region names to a common set.

## Usage

``` r
changeRegion(data, Bmat, regionVar = "region")
```

## Arguments

- data:

  Preprocessed data

- Bmat:

  Matrix of changes. Each row corresponds to a region name possibly in
  the data files, and each column corresponds to a region after mapping.
  The values in the matrix are binary. The row names and column names
  need to be specified to the region names.

- regionVar:

  String indicating the region variable. Defaults to 'region'.

## Value

Data after changing region names

## Author

Zehang Richard Li

## Examples

``` r
# Construct a small test data
testdata <- data.frame(region = c("north", "south", "east",
 "south", "east"), index = c(1:5))

# Construct a changing rule: combining south and east
Bmat <- matrix(c(1, 0, 0, 0, 1, 1), 3, 2)
colnames(Bmat) <- c("north", "south and east")
rownames(Bmat) <- c("north", "south", "east")
print(Bmat)
#>       north south and east
#> north     1              0
#> south     0              1
#> east      0              1

# New data after transformation
test <- changeRegion(testdata, Bmat, "region")
#> 2 names changed, in total 4 rows in data changed
print(test)
#>           region index
#> 1          north     1
#> 2 south and east     2
#> 3 south and east     3
#> 4 south and east     4
#> 5 south and east     5
```
