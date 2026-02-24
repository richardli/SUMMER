# Simulate hyperpriors from an GMRF

Simulate hyperpriors from an GMRF

## Usage

``` r
simhyper(
  R = 2,
  nsamp = 1e+05,
  nsamp.check = 5000,
  Amat = NULL,
  nperiod = 6,
  only.iid = TRUE
)
```

## Arguments

- R:

  Desired prior odds ratio. Default to 2, i.e., a 95% prior interval for
  the residual odds ratios lies in the interval (R, 1/R).

- nsamp:

  Sample to simulate for scaling factor

- nsamp.check:

  Sample to simulate for checking range

- Amat:

  Adjacency matrix of the areas in the data.

- nperiod:

  numerical value of how many time periods in the data

- only.iid:

  Indicator for whether or not only IID hyperpriors are simulated

## References

Wakefield, J. Multi-level modelling, the ecologic fallacy, and hybrid
study designs. *International Journal of Epidemiology*, 2009, vol. 38
(pg. 330-336).

## Author

Zehang Richard Li, Laina Mercer

## Examples

``` r
if (FALSE) { # \dontrun{
data(DemoMap)
mat <- DemoMap$Amat
priors <- simhyper(R = 2, nsamp = 1e+05, nsamp.check = 5000, Amat = mat)
} # }
```
