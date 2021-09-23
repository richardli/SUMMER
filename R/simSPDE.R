#' Simulate from the SPDE spatial model
#' 
#' Generates nCoords x nsim matrix of simulated 
#' values of the SPDE spatial process
#'
#' @param coords 2 column matrix of spatial coordinates at which to simulate the spatial process
#' @param nsim number of draws from the SPDE model
#' @param mesh SPDE mesh
#' @param effRange effective spatial range
#' @param margVar marginal variance of the spatial process
#' @param inla.seed seed input to inla.qsample. 0L sets seed intelligently, 
#' > 0 sets a specific seed, < 0 keeps existing RNG
#' 
#' @author John Paige
#' 
#' @references Lindgren, F., Rue, H., Lindström, J., 2011. An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic differential equation approach (with discussion). Journal of the Royal Statistical Society, Series B 73, 423–498.
#' 
#' @examples
#' \dontrun{
#' set.seed(123)
#' require(INLA)
#' coords = matrix(runif(10*2), ncol=2)
#' mesh = inla.mesh.2d(loc.domain=cbind(c(0, 0, 1, 1), c(0, 1, 0, 1)), 
#'   n=3000, max.n=5000, max.edge=c(.01, .05), offset=-.1)
#' simVals = simSPDE(coords, nsim=1, mesh, effRange=.2, inla.seed=1L)
#' }
#' 
#' @export
simSPDE = function(coords, nsim=1, mesh, effRange=(max(coords[,1])-min(coords[,1]))/3, margVar=1, inla.seed=0L) {
  
  # calculate SPDE model parameters based on Lindgren Rue (2015) "Bayesian Spatial Modelling with R-INLA"
  meshSize <- min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2]))))
  # it is easier to use theta and set sigma0 to 1 then to set sigma0 and the effective range directly
  # kappa0 <- sqrt(8)/effRange * meshSize # since nu = 1
  # kappa0 <- sqrt(8)/effRange # since nu = 1
  # kappa0 = sqrt(8) / 5
  # logKappa = log(kappa0)
  sigma0 = 1
  # tau0 <- 1/(sqrt(4 * pi) * kappa0 * sigma0)
  # logTau = log(tau0)
  
  # from page 5 of the paper listed above:
  logKappa = 0.5 * log(8)
  logTau = 0.5 * (lgamma(1) - (lgamma(2) + log(4*pi))) - logKappa
  theta = c(log(sqrt(margVar)), log(effRange))
  spde <- INLA::inla.spde2.matern(mesh, B.tau = cbind(logTau, -1, +1),
                            B.kappa = cbind(logKappa, 0, -1), theta.prior.mean = theta,
                            theta.prior.prec = c(0.1, 1))
  
  # generate A and Q precision matrix
  Q = INLA::inla.spde2.precision(spde, theta = theta)
  A = INLA::inla.spde.make.A(mesh, coords)
  
  # generate simulations
  simField = INLA::inla.qsample(nsim, Q, seed=inla.seed)
  simDat = as.matrix(A %*% simField)
  
  simDat
}