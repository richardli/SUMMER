# Simulate from the SPDE spatial model

Generates nCoords x nsim matrix of simulated values of the SPDE spatial
process

## Usage

``` r
simSPDE(
  coords,
  nsim = 1,
  mesh,
  eff.range = (max(coords[, 1]) - min(coords[, 1]))/3,
  marg.var = 1,
  inla.seed = 0L
)
```

## Arguments

- coords:

  2 column matrix of spatial coordinates at which to simulate the
  spatial process

- nsim:

  number of draws from the SPDE model

- mesh:

  SPDE mesh

- eff.range:

  effective spatial range

- marg.var:

  marginal variance of the spatial process

- inla.seed:

  seed input to inla.qsample. 0L sets seed intelligently, positive value
  sets a specific seed, negative value keeps existing RNG

## Details

**\[experimental\]**

## References

Lindgren, F., Rue, H., Lindström, J., 2011. An explicit link between
Gaussian fields and Gaussian Markov random fields: the stochastic
differential equation approach (with discussion). Journal of the Royal
Statistical Society, Series B 73, 423–498.

## Author

John Paige

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
require(INLA)
coords = matrix(runif(10*2), ncol=2)
mesh = inla.mesh.2d(loc.domain=cbind(c(0, 0, 1, 1), c(0, 1, 0, 1)), 
  n=3000, max.n=5000, max.edge=c(.01, .05), offset=-.1)
simVals = simSPDE(coords, nsim=1, mesh, eff.range=.2, inla.seed=1L)
} # }
```
