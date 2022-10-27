#' Benchmark posterior draws to national estimates
#' 
#' @param fitted output from \code{\link{getSmoothed}} to be benchmarked. 
#' @param national a data frame of national level estimates that is benchmarked against, with at least two columns indicating national estimates (probability scale) and the associated standard error. If benchmarking over multiple time period, a third column indicating time period is needed.  
#' @param estVar column name in \code{national} that indicates national estimates.
#' @param sdVar column name in \code{national} that indicates standard errors of national estimates.
#' @param timeVar column name in \code{national} that indicates time periods.
#' @param weight.region a data frame with a column `region` specifying subnational regions, a column `proportion` that specifies the proportion of population in each region. When multiple time periods exist, a third column `years` is required and the population proportions are the population proportions of each region in the corresponding time period.  
#' @param method a string denoting the algorithm to use for benchmarking. Options include `MH` for Metropolis-Hastings, and `Rejection` for rejection sampler. Defaults to `Rejection`.
#' 
#' @return Benchmarked object in S3 class SUMMERproj or SUMMERprojlist in the same format as the input object \code{fitted}.
#' @author Taylor Okonek, Zehang Richard Li 
#' @importFrom stats runif
#' @importFrom stats dnorm
#' @importFrom stats rbinom
#' @examples
#' \dontrun{
#' ##  ------------------------------------------ ##
#' ##     Benchmarking with smoothCluster output
#' ##  ------------------------------------------ ##
#' 
#' data(DemoData)
#' # fit unstratified cluster-level model
#' counts.all <- NULL
#' for(i in 1:length(DemoData)){
#' vars <- c("clustid", "region", "time", "age")
#' counts <- getCounts(DemoData[[i]][, c(vars, "died")], 
#' 						variables = 'died',
#' 	 			    by = vars, drop=TRUE)
#' counts$cluster <- counts$clustid
#' counts$years <- counts$time
#' counts$Y <- counts$died
#' counts$survey <- names(DemoData)[i]	
#' counts.all <- rbind(counts.all, counts)
#' }
#' periods <- c("85-89", "90-94", "95-99", "00-04", "05-09", "10-14", "15-19")
#' fit.bb  <- smoothCluster(data = counts.all, Amat = DemoMap$Amat, 
#' 				family = "betabinomial",
#' 				year_label = periods, 
#' 				survey.effect = TRUE)
#' est.bb <- getSmoothed(fit.bb, nsim = 1e4, CI = 0.95, save.draws=TRUE)
#' 
#' # construct a simple population weight data frame with equal weights
#' weight.region <- expand.grid(region = unique(counts.all$region), 
#' 						 years = periods)
#' weight.region$proportion <- 0.25
#' 
#' # construct a simple national estimates
#' national <- data.frame(years = periods, 
#' 				   est = seq(0.27, 0.1, length = 7), 
#' 				   sd = runif(7, 0.01, 0.03))
#' 
#'  # benchmarking
#' est.bb.bench <- Benchmark(est.bb, national, weight.region = weight.region, 
#' 						estVar = "est", sdVar = "sd", timeVar = "years")
#' 
#' # Sanity check: Benchmarking comparison
#' compare <- national
#' compare$before <- NA
#' compare$after <- NA
#' for(i in 1:dim(compare)[1]){
#' 	weights <- subset(weight.region, years == national$years[i])
#' 	sub <- subset(est.bb$overall, years == national$years[i])
#' 	sub <- merge(sub, weights)
#' 	sub.bench <- subset(est.bb.bench$overall, years == national$years[i])
#' 	sub.bench <- merge(sub.bench, weights)
#' 	compare$before[i] <- sum(sub$proportion * sub$median)
#' 	compare$after[i] <- sum(sub.bench$proportion * sub.bench$median)
#' }
#' plot(compare$est, compare$after, col = 2, pch = 10,
#' 		 xlim = range(c(compare$est, compare$before, compare$after)),
#' 		 ylim = range(c(compare$est, compare$before, compare$after)),
#' 		 xlab = "External national estimates", 
#' 		 ylab = "Weighted posterior median after benchmarking",
#'     main = "Sanity check: weighted average of area medians")
#' points(compare$est, compare$before)
#' abline(c(0, 1))
#' legend("topleft", c("Before benchmarking", "After benchmarking"), pch = c(1, 10), col = c(1, 2))
#' 
#' #  construct a simple national estimates
#' national <- data.frame(years = periods, 
#' 					   est = seq(0.22, 0.1, length = 7), 
#' 					   sd = runif(7, 0.01, 0.03))
#' # national does not need to have all years
#' national_sub <- national[1:3,]
#' 
#' # benchmarking
#' est.bb.bench <- Benchmark(est.bb, national_sub, 
#' 						weight.region = weight.region, 
#' 						estVar = "est", sdVar = "sd", timeVar = "years")
#' 
#' # Sanity check: only benchmarked for three periods
#' compare <- national
#' compare$before <- NA
#' compare$after <- NA
#' for(i in 1:dim(compare)[1]){
#' 	weights <- subset(weight.region, years == national$years[i])
#' 	sub <- subset(est.bb$overall, years == national$years[i])
#' 	sub <- merge(sub, weights)
#' 	sub.bench <- subset(est.bb.bench$overall, years == national$years[i])
#' 	sub.bench <- merge(sub.bench, weights)
#' 	compare$before[i] <- sum(sub$proportion * sub$median)
#' 	compare$after[i] <- sum(sub.bench$proportion * sub.bench$median)
#' }
#' plot(compare$est, compare$after, col = 2, pch = 10,
#' 		 xlim = range(c(compare$est, compare$before, compare$after)),
#' 		 ylim = range(c(compare$est, compare$before, compare$after)),
#' 		 xlab = "External national estimates", 
#' 		 ylab = "Weighted posterior median after benchmarking",
#'     main = "Sanity check: weighted average of area medians")
#' points(compare$est, compare$before)
#' abline(c(0, 1))
#' legend("topleft", c("Before benchmarking", "After benchmarking"), pch = c(1, 10), col = c(1, 2))
#' 
#' #  Another extreme benchmarking example, where almost all weights in central region
#' weight.region$proportion <- 0.01
#' weight.region$proportion[weight.region$region == "central"] <- 0.97
#' # benchmarking
#' est.bb.bench <- Benchmark(est.bb, national, weight.region = weight.region, 
#' 					estVar = "est", sdVar = "sd", timeVar = "years")
#' # It can be seen the central region are pulled to the national benchmark
#' plot(national$est, 
#' 	 subset(est.bb.bench$overall, region == "central")$mean,
#' 	 col = 2, pch = 10, xlab = "External national estimates", 
#' 	 ylab = "Central region estimates") 
#' points(national$est, 
#' 	 subset(est.bb$overall, region == "central")$mean) 
#' legend("topleft", c("Before benchmarking", "After benchmarking"), pch = c(1, 10),  col = c(1, 2))
#' abline(c(0, 1))
#' 
#' # Example with the MH method
#' # Benchmarking with MH should be applied when customized priors are 
#' #  specified for fixed effects when fitting the model
#' fit.bb.new  <- smoothCluster(data = counts.all, Amat = DemoMap$Amat, 
#' 				family = "betabinomial",
#' 				year_label = periods, 
#' 				survey.effect = TRUE, 
#' 				control.fixed = list(
#' 					mean=list(`age.intercept0:1`=-4, 
#' 						       `age.intercept1-11:1`=-5,
#' 						       `age.intercept12-23:1`=-8,
#' 						       `age.intercept24-35:1`=-9,
#' 						       `age.intercept36-47:1`=-10,
#' 						       `age.intercept48-59:1`=-11), 
#' 					prec=list(`age.intercept0:1`=10, 
#' 						       `age.intercept1-11:1`=10,
#' 						       `age.intercept12-23:1`=10,
#' 						       `age.intercept24-35:1`=10,
#' 						       `age.intercept36-47:1`=10,
#' 						       `age.intercept48-59:1`=10)))
#' est.bb.new <- getSmoothed(fit.bb.new, nsim = 10000, CI = 0.95, save.draws=TRUE)
#' 
#' #  construct a simple national estimates
#' national <- data.frame(years = periods, 
#' 					   est = seq(0.22, 0.1, length = 7), 
#' 					   sd = runif(7, 0.01, 0.03))
#' weight.region <- expand.grid(region = unique(counts.all$region), 
#' 						 years = periods)
#' weight.region$proportion <- 0.25					   
#' est.bb.bench.MH <- Benchmark(est.bb.new, national, 
#' 	weight.region = weight.region, 
#' 	estVar = "est", sdVar = "sd", timeVar = "years",
#' 	method = "MH")
#' 
#' compare <- national
#' compare$before <- NA
#' compare$after <- NA
#' for(i in 1:dim(compare)[1]){
#' 	weights <- subset(weight.region, years == national$years[i])
#' 	sub <- subset(est.bb.new$overall, years == national$years[i])
#' 	sub <- merge(sub, weights)
#' 	sub.bench <- subset(est.bb.bench.MH$overall, years == national$years[i])
#' 	sub.bench <- merge(sub.bench, weights)
#' 	compare$before[i] <- sum(sub$proportion * sub$median)
#' 	compare$after[i] <- sum(sub.bench$proportion * sub.bench$median)
#' }
#' plot(compare$est, compare$after, col = 2, pch = 10,
#' 		 xlim = range(c(compare$est, compare$before, compare$after)),
#' 		 ylim = range(c(compare$est, compare$before, compare$after)),
#' 		 xlab = "External national estimates", 
#' 		 ylab = "Weighted posterior median after benchmarking",
#'     main = "Sanity check: weighted average of area medians")
#' points(compare$est, compare$before)
#' abline(c(0, 1))
#' legend("topleft", c("Before benchmarking", "After benchmarking"), pch = c(1, 10), col = c(1, 2))
#' 
#' ##  ------------------------------------------ ##
#' ##     Benchmarking with smoothDirect output
#' ##  ------------------------------------------ ##
#' years <- levels(DemoData[[1]]$time)
#' # obtain direct estimates
#' data_multi <- getDirectList(births = DemoData, years = years,
#'                         regionVar = "region",  timeVar = "time", clusterVar = "~clustid+id",
#'                         ageVar = "age", weightsVar = "weights", geo.recode = NULL)
#' data <- aggregateSurvey(data_multi)
#' #  subnational model
#' years.all <- c(years, "15-19")
#' fit2 <- smoothDirect(data = data, Amat = DemoMap$Amat,
#'                  year_label = years.all, year_range = c(1985, 2019),
#'                  time.model = "rw2", m = 5, type.st = 4)
#' out2a <- getSmoothed(fit2, joint = TRUE, nsim = 1e5, save.draws = TRUE)
#' 
#' ##
#' ## Benchmarking for yearly estimates
#' ##
#' weight.region <- expand.grid(region = unique(data$region[data$region != "All"]),
#'                              years = 1985:2019)
#' weight.region$proportion <- 0.25
#' # construct a simple national estimates
#' national <- data.frame(years = 1985:2019,
#'                        est = seq(0.25, 0.15, length = 35),
#'                        sd = runif(35, 0.03, 0.05))
#' # Benchmarking to national estimates on the yearly scale
#' out2b <- Benchmark(out2a, national, weight.region = weight.region,
#'                           estVar = "est", sdVar = "sd", timeVar = "years")
#' plot(out2a$overall)  
#' plot(out2b$overall) 
#' 
#' # combine the point estimate and compare with the benchmark values
#' national.est <- aggregate(mean ~ years, 
#'    data = out2a$overall[out2a$overall$is.yearly, ], FUN = mean)
#' national.est.bench <- aggregate(mean ~ years, 
#'    data = out2b$overall[out2b$overall$is.yearly, ], FUN = mean)
#' 
#' plot(national$est, national.est$mean,  
#' 		 xlim = range(c(national$est, national.est$mean, national.est.bench$mean)),
#' 		 ylim = range(c(national$est, national.est$mean, national.est.bench$mean)),
#' 		 xlab = "External national estimates", 
#' 		 ylab = "Weighted posterior median after benchmarking",
#'     main = "Sanity check: weighted average of area means")
#' points(national$est, national.est.bench$mean, col = 2, pch = 10)
#' abline(c(0, 1))
#' legend("topleft", c("Before benchmarking", "After benchmarking"), pch = c(1, 10), col = c(1, 2))
#' 
#' 
#' 
#' ##
#' ## Benchmarking for period estimates
#' ##
#' weight.region <- expand.grid(region = unique(data$region[data$region != "All"]),
#'                              years = years.all)
#' weight.region$proportion <- 0.25
#' # construct a simple national estimates
#' national <- data.frame(years = years.all,
#'                        est = seq(0.25, 0.15, len = 7),
#'                        sd = runif(7, 0.01, 0.03))
#' # Benchmarking to national estimates on the period scale
#' out2c <- Benchmark(out2a, national, weight.region = weight.region,
#'                           estVar = "est", sdVar = "sd", timeVar = "years")
#' plot(out2a$overall)
#' plot(out2c$overall)
#' 
#' # combine the point estimate and compare with the benchmark values
#' national.est <- aggregate(mean ~ years, 
#' 			data = out2a$overall[!out2a$overall$is.yearly, ], FUN = mean)
#' national.est.bench <- aggregate(mean ~ years, 
#' 			data = out2c$overall[!out2b$overall$is.yearly, ], FUN = mean)
#' 
#' plot(national$est, national.est$mean,  
#' 		 xlim = range(c(national$est, national.est$mean, national.est.bench$mean)),
#' 		 ylim = range(c(national$est, national.est$mean, national.est.bench$mean)),
#' 		 xlab = "External national estimates", 
#' 		 ylab = "Weighted posterior median after benchmarking",
#'     main = "Sanity check: weighted average of area means")
#' points(national$est, national.est.bench$mean, col = 2, pch = 10)
#' abline(c(0, 1))
#' legend("topleft", c("Before benchmarking", "After benchmarking"), pch = c(1, 10), col = c(1, 2))
#' 
#' 
#'  }
#' @export
#' 

Benchmark <- function(fitted, national, estVar, sdVar, timeVar = NULL, weight.region = NULL, method = c("MH","Rejection")[2]) {

	############################################################################
	#  Helper functions
	############################################################################

	# This function grabs the fixed effect draws, 
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
	                   nregion) {

	  ntime <- length(old_thetas)

	  # initialize numerator and denominator
	  num <- 1
	  denom <- 1
	  
	  # add likelihoods for z
	  for (i in 1:ntime) {
	    num <- num + dnorm(z[i], new_thetas[i], sqrt(var_z[i]), log = TRUE)
	    denom <- denom + dnorm(z[i], old_thetas[i], sqrt(var_z[i]), log = TRUE)
	  }
	  
	  # add likelihoods for intercepts
	  for (i in 1:length(old_betas)) {
	    num <- num + dnorm(old_betas[i], intercept_means[i], sqrt(var_plus[i]), log = TRUE)
	    denom <- denom + dnorm(new_betas[i], intercept_means[i], sqrt(var_plus[i]), log = TRUE)
	  }
	  # print(num/denom)
	  return(min(1, exp(num - denom)))
	}

	############################################################################
	#  Check input format
	############################################################################

	# Check input types
	if(!is(fitted, "SUMMERprojlist")){
		if(is(fitted, "SUMMERproj")){
			stop("Please use the full returned object from getSmoothed() with 'save.draws = TRUE'")
		}else{
			stop("The argument 'fitted' needs to be a returned object from getSmoothed() with 'save.draws = TRUE'")
		}
	}
	if(!is.null(fitted$final)){
		stop("Benchmarking on multiple-frame model not implemented yet.")
	}
	if(is.null(fitted$CI)){
		stop("Please rerun the getSmoothed() function with a valid CI.")
	}else{
		lowerCI <- (1 - fitted$CI) / 2
		upperCI <- 1 - (1 - fitted$CI) / 2
	}
  if(!(method %in% c("MH","Rejection"))) {
    stop(paste0("Method ", method, " is not available. Choose one of 'MH', 'Rejection'."))
  }
  message(paste0("Method ", method, " being used."))

  years <- fitted$overall$years[match(1:max(fitted$overall$time), fitted$overall$time)]

	# Message which is used if both yearly and periods exist
  # check which part matches the national estimates years
  # if m = 1, is.yearly = FALSE for all rows here so it won't enter the if statement
	skip <- rep(FALSE, dim(fitted$overall)[1])
	fitted_overall_orig <- fitted$overall
	if(length(unique(fitted$overall$is.yearly)) == 2){
		years.national <- national$years 
		years.exist <- fitted$overall$years %in% years.national
		if(length(unique(fitted$overall$is.yearly[years.exist])) == 2){
			stop("The national estimates and model estimates contain both yearly and period estimates. Only one type of estimates can be used for benchmarking.")
		}else if(unique(fitted$overall$is.yearly[years.exist])[1]){
			message("'yearly' estimates are benchmarked")
		}else{
			message("'period' estimates are benchmarked")
		}
		skip <- fitted$over$is.yearly != unique(fitted$overall$is.yearly[years.exist])[1]
		fitted$overall <- fitted$overall[!skip, ]
	}
	# Check the population weight matrix
	is.time <- length(unique(fitted$overall$years)) > 1
	ntime <- max(fitted$overall$time)
	nregion <- max(fitted$overall$area)
	if(is.null(weight.region)){
		stop("weight.region argument is required.")
	}
	if("region" %in% colnames(weight.region) == FALSE){
		stop("There is no 'region' column in weight.region.")
	}
	if("years" %in% colnames(weight.region) == FALSE && is.time){
		stop("There is no 'years' column in weight.region.")
	}
	if("proportion" %in% colnames(weight.region) == FALSE){
		if(is.time && dim(weight.region)[2] == 3){
			tmp <- colnames(weight.region)[!colnames(weight.region) %in% c("region", "years")]			
			weight.region$proportion <- weight.region[, tmp] 
			message(paste0("Use ", tmp, " column as the population proportion."))
		}else if(!is.time && dim(weight.region)[2] == 2){
			tmp <- colnames(weight.region)[!colnames(weight.region) %in% c("region")]			
			weight.region$proportion <- weight.region[, tmp] 
			message(paste0("Use ", tmp, " column as the population proportion."))
		}else{
			stop("Cannot determine the column indicating population proportion. Specify the column name to be 'proportion'.")
		}
	}
	if(sum(!fitted$overall$region %in% weight.region$region) > 0 && nregion > 1){
		stop("weight.region$region does not contain all the regions in the fitted model.")
	}
	if(is.time){
		if(sum(!weight.region$years %in% national[, timeVar]) > 0){
			weight.region <- weight.region[weight.region$years %in% national[, timeVar], ]
		}
		if(sum(!fitted$overall$years %in% weight.region$years) > 0){
			tmp <- unique(fitted$overall$years[fitted$overall$years %in% weight.region$years])
			warning(paste0("weight.region$years does not contain all the time periods in the fitted model. Benchmarking only performed for the following period:\n", paste(tmp, collapse = ", ")))
		}
	}
	if(is.time){
		for(j in unique(weight.region$years)){
			sub <- which(weight.region$years == j)
			if(abs(1 - sum(weight.region$proportion[sub])) > 0.001){
				stop(paste0("Population proportion in ", j, " does not sum to 1."))
			}
		}
	}else{
		if(abs(1 - sum(weight.region$proportion)) > 0.001){
				stop(paste0("Population proportion does not sum to 1."))
		}	
	}
	control <- fitted$control.fixed
	if(length(control) == 0 && method == "MH"){
		method <- "Rejection"
		warning("The model is fitted with default priors for the intercept, so benchmarking using the rejection algorithm is used instead.")
	}


	############################################################################
	#  Arrange input data
	############################################################################

	# Match benchmarked values
	t_sub <- 1:ntime
	nat <- data.frame(est = national[, estVar], 
				sd = national[, sdVar])
	if(!is.null(timeVar)){
		nat$years = national[, timeVar]
		order <- match(years, nat$years)
		nat <- nat[order[which(!is.na(order))], ]
		t_sub <- match(nat$years, years)
	}

	# Sampling and save message
	n0 <- length(fitted$draws)
	q_mat <- matrix(0, nrow = ntime, ncol = n0)

	# fill in q_mat
	if(nregion == 1){
		for(i in 1:length(fitted$draws.est.overall)){
			if(!skip[i]){
				index <- which(years == fitted$draws.est.overall[[i]]$years)
				q_mat[index, ] <- fitted$draws.est.overall[[i]]$draws				
			}
		}
	}else{
		for(i in 1:length(fitted$draws.est.overall)){
			if(!skip[i]){
				index <- which(years == fitted$draws.est.overall[[i]]$years)
				tmp <- which(weight.region$years == fitted$draws.est.overall[[i]]$years &
							 weight.region$region == fitted$draws.est.overall[[i]]$region)
				if(length(tmp) == 0) next
				tmp_pop <- weight.region$proportion[tmp]
				q_mat[index, ] <- q_mat[index, ] + fitted$draws.est.overall[[i]]$draws * tmp_pop
			}
		}
	}

	# if benchmarking only happens in subset of years
	q_mat <- q_mat[t_sub, ]


	############################################################################
	#  Rejection algorithm
	############################################################################
  if (method == "Rejection") {
	  # generate Us
	  U <- runif(n0, 0, 1)
	  
	  # set up fitted_list and prop_accepted for return
	  fitted_list <- list()
	  prop_accepted <- 0
	  accept_ratio <- matrix(0, nrow = nrow(q_mat), ncol = ncol(q_mat))
	  for (i in 1:dim(q_mat)[1]) {
	    for (j in 1:ncol(q_mat)) {
	      accept_ratio[i, j] <- exp((-1 / (2 * nat$sd[i] ^ 2)) * (q_mat[i, j] - nat$est[i]) ^ 2)
	    }
	  }
	  
	  # get accepted samples if we take all time points at once
	  multi_accepted_samps <- rep(FALSE, ncol(accept_ratio))
	  for (i in 1:ncol(accept_ratio)) {
	    multi_accepted_samps[i] <- ( U[i] < prod(accept_ratio[,i]) )
	  }
	  acc <- which(multi_accepted_samps == TRUE)
	  
	  if(length(acc) == 0){
	    stop("All posterior samples have been rejected. Please rerun getSmoothed() with a larger 'nsim' argument.")
	  } 
	  
	  fitted$msg <- paste0(fitted$msg, 
	                       "\nThe posterior draws have been benchmarked to external information. The acceptance ratio is ", 
	                       round(length(acc) / n0, 3),
	                       "\n")
	  message("The posterior draws have been benchmarked to external information. The acceptance ratio is ", 
	          round(length(acc) / n0, 3),
	          "\n")
	  
	############################################################################
	#  Metropolis-Hastings algorithm
	############################################################################
	} else {
	  # get intercept draws and priors for fixed effects
	  ints <- get_intercepts(fitted)
	  if (any(ints$priors[,2] <= 0)) {
	    stop("fitted model cannot have fixed effect precisions <= 0")
	  }
	  
	  # Initialize (beta^0, theta^0)
	  old_thetas <- q_mat[,1]
	  old_betas <- ints$draws[,1]
	  
	  accepted_ids <- c()
	  num_accepted <- 0
	  for (i in 2:ncol(q_mat)) {
	    # Sample new (beta', theta')
	    new_thetas <- q_mat[,i]
	    new_betas <- ints$draws[,i]
	    
	    # compute A
	    accept_prob <- A_full(old_thetas = old_thetas,
	                          old_betas = old_betas,
	                          new_thetas = new_thetas,
	                          new_betas = new_betas,
	                          z = nat$est,
	                          intercept_means = ints$priors[,1], 
	                          var_z = nat$sd^2,
	                          var_plus = 1/ints$priors[,2],
	                          nregion = nregion)
	    
	    accept_yn <- rbinom(n = 1, size = 1, prob = accept_prob)
	    
	    if (accept_yn) {
	      accepted_ids <- c(accepted_ids, i)
	      old_betas <- new_betas
	      old_thetas <- new_thetas
	      num_accepted <- num_accepted + 1
	      # print(i)
	    } else {
	      # if we're not at the initialization value...
	      if (i > 2) {
	        accepted_ids <- c(accepted_ids, tail(accepted_ids,1))
	      }
	    }
	  }
	  
	  # messaging
	  if(num_accepted == 0){
	    stop("All posterior samples have been rejected. Please rerun getSmoothed() with a larger 'nsim' argument.")
	  } 
	  
	  fitted$msg <- paste0(fitted$msg, 
	                       "\nThe posterior draws have been benchmarked to external information. The acceptance ratio is ", 
	                       round(num_accepted / n0, 3),
	                       "\n")
	  message("The posterior draws have been benchmarked to external information. The acceptance ratio is ", 
	          round(num_accepted / n0, 3),
	          "\n")
	  
	  acc <- accepted_ids
	  
	}
	

	# Update fitted
	for(i in 1:length(fitted$draws.est )){
		fitted$draws.est[[i]]$draws <- fitted$draws.est[[i]]$draws[acc]  
	}
	for(i in 1:length(fitted$draws.est.overall )){
		fitted$draws.est.overall[[i]]$draws <- fitted$draws.est.overall[[i]]$draws[acc]  
	}
	draws.raw <- NULL
	counter <- 1
	for(i in acc){
		draws.raw[[counter]] <- fitted$draws[[i]]
		counter <- counter + 1
	}
	fitted$draws <- draws.raw

	if(!is.null(fitted$stratified)){
	fitted$stratified$variance <- fitted$stratified$median <- fitted$stratified$mean <-  fitted$stratified$lower <-  fitted$stratified$upper <- NA
	for(i in 1:length(fitted$draws.est)){
		y <- fitted$draws.est[[i]]$years
		r <- fitted$draws.est[[i]]$region
		s <- fitted$draws.est[[i]]$strata
		j <- which(fitted$stratified$years == y & fitted$stratified$region == r & fitted$stratified$strata == s)
		fitted$stratified[j, c("lower", "median", "upper")] <- stats::quantile(fitted$draws.est[[j]]$draws, c(lowerCI, 0.5, upperCI))
        fitted$stratified[j, "mean"] <- mean(fitted$draws.est[[j]]$draws)
        fitted$stratified[j, "variance"] <- var(fitted$draws.est[[j]]$draws)
	}
	}

  fitted$overall <- fitted_overall_orig
	fitted$overall$variance <- fitted$overall$median <- fitted$overall$mean <-  fitted$overall$lower <-  fitted$overall$upper <- NA
	for(i in 1:length(fitted$draws.est.overall)){
		y <- fitted$draws.est.overall[[i]]$years
		r <- fitted$draws.est.overall[[i]]$region
		j <- which(fitted_overall_orig$years == y & fitted_overall_orig$region == r)
		fitted$overall[j, c("lower", "median", "upper")] <- stats::quantile(fitted$draws.est.overall[[j]]$draws, c(lowerCI, 0.5, upperCI))
        fitted$overall[j, "mean"] <- mean(fitted$draws.est.overall[[j]]$draws)
        fitted$overall[j, "variance"] <- var(fitted$draws.est.overall[[j]]$draws)
	}

	fitted$nsim <- length(acc)
	if (method == "Rejection") {
	  fitted$accept.ratio <- length(acc) / n0
	} else {
	  fitted$accept.ratio <- num_accepted / n0
	}
	fitted$benchmarked <- TRUE

	return(fitted)
}

