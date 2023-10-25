#' Calculate the mean of a distribution whose 
#' logit is Gaussian
#' 
#' Adapted from logitnorm package.  Calculates the mean of a distribution whose 
#' logit is Gaussian. Each row of muSigmaMat is a mean and standard deviation 
#' on the logit scale. 
#' 
#' @param muSigmaMat An n x 2 matrix where each row is \eqn{\mu}{mu} and \eqn{\sigma}{sigma} 
#' on the logit scale for an independent random variable.
#' @param logisticApprox Whether or not to use logistic approximation to speed 
#' up computation. See details for more information.
#' @param ... More arguments, passed to \code{integrate} function
#' 
#' @return A vector of expectations of the specified random variables
#' 
#' @author John Paige
#' 
#' @details If \eqn{\mbox{logit}(Y) \sim N(\mu, \sigma^2)}{logit(Y) ~ N(mu, sigma^2)}, 
#' This function calculates \eqn{E[Y]}{E[Y]} via either numerical integration or by 
#' assuming that Y follows a logistic distribution. Under this approximation, setting 
#' \eqn{k = 16 \sqrt(3) / (15 \pi)}{k = 16 * sqrt(3) / (15 * pi)}, we approximate 
#' the expectation as:
#' \deqn{E[Y] = expit(\mu / \sqrt(1 + k^2 \sigma^2))}{E[Y] = expit(mu / sqrt(1 + k^2 * sigma^2))}.
#' The above logistic approximation speeds up the computation, but also sacrifices 
#' some accuracy.
#' 
#' @examples
#' mus = c(-5, 0, 5)
#' sigmas = rep(1, 3)
#' logitNormMean(cbind(mus, sigmas))
#' logitNormMean(cbind(mus, sigmas), TRUE)
#' 
#' @export
logitNormMean = function(muSigmaMat, logisticApprox=FALSE, ...) {
  if(length(muSigmaMat) > 2) {
    apply(muSigmaMat, 1, logitNormMean, logisticApprox=logisticApprox, ...)
  }
  else {
    mu = muSigmaMat[1]
    sigma = muSigmaMat[2]
    if(sigma == 0)
      SUMMER::expit(mu)
    else {
      if(any(is.na(c(mu, sigma))))
        NA
      else if(!logisticApprox) {
        # numerically calculate the mean
        fExp <- function(x) exp(stats::plogis(x, log.p=TRUE) + stats::dnorm(x, mean = mu, sd = sigma, log=TRUE))
        stats::integrate(fExp, mu-10*sigma, mu+10*sigma, abs.tol = 0, ...)$value
      } else {
        # use logistic approximation
        k = 16 * sqrt(3) / (15 * pi)
        SUMMER::expit(mu / sqrt(1 + k^2 * sigma^2))
      }
    }
  }
}