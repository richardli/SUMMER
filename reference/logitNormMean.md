# Calculate the mean of a distribution whose logit is Gaussian

Adapted from logitnorm package. Calculates the mean of a distribution
whose logit is Gaussian. Each row of muSigmaMat is a mean and standard
deviation on the logit scale.

## Usage

``` r
logitNormMean(muSigmaMat, logisticApprox = FALSE, ...)
```

## Arguments

- muSigmaMat:

  An n x 2 matrix where each row is \\\mu\\ and \\\sigma\\ on the logit
  scale for an independent random variable.

- logisticApprox:

  Whether or not to use logistic approximation to speed up computation.
  See details for more information.

- ...:

  More arguments, passed to `integrate` function

## Value

A vector of expectations of the specified random variables

## Details

If \\\mbox{logit}(Y) \sim N(\mu, \sigma^2)\\, This function calculates
\\E\[Y\]\\ via either numerical integration or by assuming that Y
follows a logistic distribution. Under this approximation, setting \\k =
16 \sqrt{3} / (15 \pi)\\, we approximate the expectation as: \$\$E\[Y\]
= expit(\mu / \sqrt{1 + k^2 \sigma^2})\$\$ The above logistic
approximation speeds up the computation, but also sacrifices some
accuracy.

## Author

John Paige

## Examples

``` r
mus = c(-5, 0, 5)
sigmas = rep(1, 3)
logitNormMean(cbind(mus, sigmas))
#> [1] 0.01079679 0.50000000 0.98920321
logitNormMean(cbind(mus, sigmas), TRUE)
#> [1] 0.01325606 0.50000000 0.98674394
```
