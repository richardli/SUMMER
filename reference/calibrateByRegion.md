# Calibrate the point level totals so their sum matches the regional totals

Calibrate/normalize the point level totals so their sum matches the
regional totals. Technically, the totals can be at any level smaller
than the region level specified.

## Usage

``` r
calibrateByRegion(point.totals, point.regions, regions, region.totals)
```

## Arguments

- point.totals:

  Vector of point level totals that will be calibrated/normalized

- point.regions:

  Vector of regions associated with each point

- regions:

  Vector of region names

- region.totals:

  Vector of desired region level totals associated with `regions`

## Value

A vector of same length as point.totals and point.regions containing the
calibrated/normalized point totals that sum to the correct regional
totals

Vector of updated point level totals, calibrated to match region totals

## Details

**\[experimental\]**

## Author

John Paige

## Examples

``` r
point.totals = c(1, 1, 1, 2)
point.regions = c("a", "a", "b", "b")
region.totals = c(10, 20)
regions = c("a", "b")
calibrateByRegion(point.totals, point.regions, regions, region.totals)
#> [1]  5.000000  5.000000  6.666667 13.333333
```
