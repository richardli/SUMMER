# Aggregate person-month data into counts and totals by groups.

Aggregate person-month data into counts and totals by groups.

## Usage

``` r
getCounts(data, variables, by, ignore = NULL, addtotal = TRUE, drop = TRUE)
```

## Arguments

- data:

  dataset in person-month format

- variables:

  a character vector of the variables to aggregate

- by:

  a character vector of columns that specifies which groups to aggregate
  by.

- ignore:

  list of conditions not to impute 0. If left unspecified, any group
  levels not in the data will be imputed to have 0 counts.

- addtotal:

  logical indicator of whether to add a column of group total counts.

- drop:

  logical indicator of whether to drop all rows with total = 0.

## Value

data.frame of the ggregated counts.

## Author

Zehang Richard Li

## Examples

``` r
 
# a toy dataset with 4 time periods but one missing in data
timelist <- factor(1:4)
data = data.frame(died = c(0,0,0,1,1,0,0), 
          area = c(rep(c("A", "B"), 3), "A"), 
          time = timelist[c(1,1,2,3,3,3,3)])
data
#>   died area time
#> 1    0    A    1
#> 2    0    B    1
#> 3    0    A    2
#> 4    1    B    3
#> 5    1    A    3
#> 6    0    B    3
#> 7    0    A    3
# without ignore argument, all levels will be imputed
getCounts(data, variables = "died", by = c("area", "time"))
#>   area time died total
#> 1    A    1    0     1
#> 2    B    1    0     1
#> 3    A    2    0     1
#> 4    A    3    1     2
#> 5    B    3    1     2

# ignoring time = 4, the ignored level will not be imputed (but still in the output)
getCounts(data, variables = "died", by = c("area", "time"), 
      ignore = list("time"=c(4)) )
#>   area time died total
#> 1    A    1    0     1
#> 2    B    1    0     1
#> 3    A    2    0     1
#> 4    A    3    1     2
#> 5    B    3    1     2

 
```
