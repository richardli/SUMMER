# Simulate spatial and temporal random effects

This function simulates spatial and temporal random effects with mean
zero. The method is described in Algorithm 3.1 of Rue & Held 2005.

## Usage

``` r
rst(
  n = 1,
  type = c("s", "t", "st")[1],
  type.s = "ICAR",
  type.t = c("RW1", "RW2")[2],
  Amat = NULL,
  n.t = NULL,
  scale.model = TRUE
)
```

## Arguments

- n:

  sample size

- type:

  type of random effects: temporal (t), spatial (s), or spatial-temporal
  (st)

- type.s:

  type of spatial random effect, currently only ICAR is available

- type.t:

  type of temporal random effect, currently only RW1 and RW2 are
  available

- Amat:

  adjacency matrix for the spatial regions

- n.t:

  number of time points for the temporal random effect

- scale.model:

  logical indicator of whether to scale the random effects to have unit
  generalized variance. See Sørbye 2013 for more details

## Value

a matrix (for spatial or temporal) or a three-dimensional array (for
spatial-temporal) of the random effects.

## References

Rue, H., & Held, L. (2005). *Gaussian Markov random fields: theory and
applications*. CRC press.

Sørbye, S. H. (2013). *Tutorial: Scaling IGMRF-models in R-INLA*.
Department of Mathematics and Statistics, University of Tromsø.

## Author

Zehang Richard Li

## Examples

``` r
if (FALSE) { # \dontrun{
data(DemoMap)
## Spatial random effects 
out <- rst(n=10000, type = "s", Amat = DemoMap$Amat)
# To verify the mean under the conditional specification
mean(out[,1] - apply(out[,c(2,3,4)], 1, mean))  
mean(out[,2] - apply(out[,c(1,3)], 1, mean)) 
mean(out[,3] - apply(out[,c(1,2,4)], 1, mean))  
mean(out[,4] - apply(out[,c(1,3)], 1, mean)) 

## Temporal random effects (RW1)
out <- rst(n=1, type = "t", type.t = "RW1", n.t = 200, scale.model = FALSE)
par(mfrow = c(1,2))
plot(1:dim(out)[2], out, col = 1, type = "l", xlab = "Time", ylab = "Random effects")
# verify the first order difference is normally distributed
first_diff <- diff(as.numeric(out[1,]))
qqnorm(first_diff )  
abline(c(0,1))

## Temporal random effects (RW2)
out <- rst(n=1, type = "t", type.t = "RW2", n.t = 200, scale.model = FALSE)
par(mfrow = c(1,2))
plot(1:dim(out)[2], out, col = 1, type = "l", xlab = "Time", ylab = "Random effects")
# verify the second order difference is normally distributed
first_diff <- diff(as.numeric(out[1,]))
second_diff <- diff(first_diff)
qqnorm(second_diff)  
abline(c(0,1))

## Spatial-temporal random effects
out <- rst(n=1, type = "st", type.t = "RW2", Amat = DemoMap$Amat, n.t = 50)
dimnames(out)
par(mfrow = c(1,1))
plot(1:dim(out)[3], out[1,1,], col = 1,
 type = "l", ylim = range(out), xlab = "Time", ylab = "Random effects")
for(i in 2:4) lines(1:dim(out)[3], out[1,i,], col = i)
legend("bottomright", colnames(DemoMap$Amat), col = c(1:4), lty = rep(1,4))
} # }
```
