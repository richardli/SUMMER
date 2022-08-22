# This can be a function within Benchmark() to get the fixed effect draws, 
#    this only applies to the fixed effects specified in control.fixed
# input: from getSmoothed()
# output: 
#   1. draws: K by Nsim matrix of draws of K named intercepts
#   2. priors: K by 2 matrix of prior mean and precision of the K named intercepts 
get_intercepts <- function(fitted){
  control <- fitted$control.fixed
  K <- length(control$mean)
  priors <- matrix(NA, K, 2)
  colnames(priors) <- c("mean", "prec")
  names <- names(control$mean)
  rownames(priors) <- names
  for(i in 1:K){
    priors[i, 1] <- control$mean[[i]]
    priors[i, 2] <- control$prec[[names[i]]]
  }	
  Nsim <-length(fitted$draws)
  draws <- matrix(NA, K, Nsim)
  rownames(draws) <- names
  for(i in 1:Nsim){
    draws[, i] <- fitted$draws[[i]]$latent[names, ]
  }
  colnames(draws) <- paste0("sample:", 1:Nsim)
  return(list(draws = draws, priors = priors))
  
}

## A_full: the function that calculates the acceptance probability given:
# old_thetas: current fitted values, on probability scale (not logit scale)
# new_thetas: proposed fitted values, on probability scale (not logit scale)
# old_betas: current intercepts
# new_betas: proposed intercepts
# z: national benchmarks, in order arrange(time)
# intercept_means: the means for the intercept terms that were input to inla()
# weights: population size weights, in order arrange(region, time), where they sum to 1 at each time point
# var_z: variances for the national benchmarks
# var_plus: the variances for the intercept terms that were input to inla() (VARIANCES, not precisions or sds)
# nregion: number of regions
# ntime: number of time points
# assumes thetas are in order arrange(region, time)
A_full <- function(old_thetas, 
                   new_thetas, 
                   old_betas,
                   new_betas,
                   z, 
                   intercept_means,
                   var_z, 
                   var_plus,
                   nregion,
                   ntime) {
  
  # initialize numerator and denominator
  num <- 1
  denom <- 1
  
  # add likelihoods for z
  for (i in 1:ntime) {
    num <- num * dnorm(z[i], new_thetas[i], sqrt(var_z[i]))
    denom <- denom * dnorm(z[i], old_thetas[i], sqrt(var_z[i]))
  }
  
  # add likelihoods for intercepts
  for (i in 1:length(old_betas)) {
    num <- num * dnorm(old_betas[i], intercept_means[i], sqrt(var_plus[i]))
    denom <- denom * dnorm(new_betas[i], intercept_means[i], sqrt(var_plus[i]))
  }
  # print(num/denom)
  return(min(1, num/denom))
}
