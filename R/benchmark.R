#' Benchmark posterior draws to national estimates
#' 
#' @param fitted output from \code{\link{getSmoothed}} to be benchmarked. 
#' @param national a data frame of national level estimates that is benchmarked against, with at least two columns indicating national estimates (probability scale) and the associated standard error. If benchmarking over multiple time period, a third column indicating time period is needed.  
#' @param estVar column name in \code{national} that indicates national estimates.
#' @param sdVar column name in \code{national} that indicates standard errors of national estimates.
#' @param timeVar column name in \code{national} that indicates time periods.
#' @param weight.region a data frame with a column `region` specifying subnational regions, a column `proportion` that specifies the proportion of population in each region. When multiple time periods exist, a third column `years` is required and the population proportions are the population proportions of each region in the corresponding time period.  
#' 
#' @return Benchmarked object in S3 class SUMMERproj or SUMMERprojlist in the same format as the input object \code{fitted}.
#' @author Taylor Okonek, Zehang Richard Li 
#' @importFrom stats runif
#' @examples
#' \dontrun{
#' data(DemoData)
# fit unstratified cluster-level model
#'counts.all <- NULL
#'for(i in 1:length(DemoData)){
#'	vars <- c("clustid", "region", "time", "age")
#'  counts <- getCounts(DemoData[[i]][, c(vars, "died")], 
#' 						variables = 'died',
#'		 			    by = vars, drop=TRUE)
#'	counts$cluster <- counts$clustid
#'	counts$years <- counts$time
#'	counts$Y <- counts$died
#'	counts$survey <- names(DemoData)[i]	
#'	counts.all <- rbind(counts.all, counts)
#'}
#'periods <- c("85-89", "90-94", "95-99", "00-04", "05-09", "10-14", "15-19")
#'fit.bb  <- smoothCluster(data = counts.all, Amat = DemoMap$Amat, 
#'					family = "betabinomial",
#'					year_label = periods, 
#'					survey.effect = TRUE)
#'est.bb <- getSmoothed(fit.bb, nsim = 1000, CI = 0.95, save.draws=TRUE)
#'
#'# construct a simple population weight data frame with equal weights
#'weight.region <- expand.grid(region = unique(counts.all$region), 
#'							 years = periods)
#'weight.region$proportion <- 0.25
#'
#'# construct a simple national estimates
#'national <- data.frame(years = periods, 
#'					   est = seq(0.25, 0.15, length = 7), 
#'					   sd = runif(7, 0.01, 0.03))
#'# benchmarking
#'est.bb.bench <- Benchmark(est.bb, national, weight.region = weight.region, 
#'						estVar = "est", sdVar = "sd", timeVar = "years")
#'
#'# Sanity check
#'compare <- national
#'compare$before <- NA
#'compare$after <- NA
#'for(i in 1:dim(compare)[1]){
#'	weights <- subset(weight.region, years == national$years[i])
#'	sub <- subset(est.bb$overall, years == national$years[i])
#'	sub <- merge(sub, weights)
#'	sub.bench <- subset(est.bb.bench$overall, years == national$years[i])
#'	sub.bench <- merge(sub.bench, weights)
#'	compare$before[i] <- sum(sub$proportion * sub$median)
#'	compare$after[i] <- sum(sub.bench$proportion * sub.bench$median)
#'}
#'plot(compare$est, compare$after, col = 2, pch = 10,
#'		 xlim = range(c(compare$est, compare$before, compare$after)),
#'		 ylim = range(c(compare$est, compare$before, compare$after)),
#'		 xlab = "External national estimates", 
#'		 ylab = "Weighted posterior median after benchmarking",
#'       main = "Sanity check: weighted average of area medians")
#'points(compare$est, compare$before)
#'abline(c(0, 1))
#' legend("topleft", c("Before benchmarking", "After benchmarking"), pch = c(1, #'10), col = c(1, 2))
#'
# Another extreme benchmarking example, where almost all weights in central #'region
#'weight.region$proportion <- 0.01
#'weight.region$proportion[weight.region$region == "central"] <- 0.97
#'# benchmarking
#'est.bb.bench <- Benchmark(est.bb, national, weight.region = weight.region, 
#'						estVar = "est", sdVar = "sd", timeVar = "years")
#'# It can be seen the central region are pulled to the national benchmark
#'subset(est.bb.bench$overall, region == "central")
#'subset(est.bb$overall, region == "central") 
#'  }
#' @export
#' 

Benchmark <- function(fitted, national, estVar, sdVar, timeVar = NULL, weight.region = NULL){
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

	# Check the population weight matrix
	is.time <- length(unique(fitted$stratified$years)) > 1
	ntime <- max(fitted$overall$time)
	nregion <- max(fitted$overall$region)
	years <- fitted$overall$years[match(1:max(fitted$overall$time), fitted$overall$time)]
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
	if(sum(!fitted$stratified$region %in% weight.region$region) > 0 && nregion > 1){
		stop("weight.region$region does not contain all the regions in the fitted model.")
	}
	if(is.time){
		if(sum(!fitted$stratified$years %in% weight.region$years) > 0){
			stop("weight.region$years does not contain all the time periods in the fitted model.")
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


	# Match benchmarked values
	nat <- list(est = national[, estVar], 
				sd = national[, sdVar])
	if(!is.null(timeVar)){
		nat$years = national[, timeVar]
	}



	# Rejection sampling and save message
	n0 <- length(fitted$draws)
	q_mat <- matrix(0, nrow = ntime, ncol = n0)

	# fill in q_mat
	if(nregion == 1){
		for(i in 1:length(fitted$draws.est.overall)){
			index <- which(years == fitted$draws.est.overall[[i]]$years)
			q_mat[index, ] <- fitted$draws.est.overall[[i]]$draws
		}
	}else{
		for(i in 1:length(fitted$draws.est.overall)){
			index <- which(years == fitted$draws.est.overall[[i]]$years)
			tmp <- which(weight.region$years == fitted$draws.est.overall[[i]]$years &
							 weight.region$region == fitted$draws.est.overall[[i]]$region)
			tmp_pop <- weight.region$proportion[tmp]
			q_mat[index, ] <- q_mat[index, ] + fitted$draws.est.overall[[i]]$draws * tmp_pop
		}
	}

	# generate Us
	U <- runif(n0, 0, 1)

	# set up fitted_list and prop_accepted for return
	fitted_list <- list()
	prop_accepted <- 0
	accept_ratio <- matrix(0, nrow = nrow(q_mat), ncol = ncol(q_mat))
	for (i in 1:ntime) {
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

	fitted$msg <- paste0(fitted$msg, 
						"\nThe posterior draws have been benchmarked to external information. The acceptance ratio is ", 
						round(length(acc) / n0, 3),
						"\n")
	message("The posterior draws have been benchmarked to external information. The acceptance ratio is ", 
						round(length(acc) / n0, 3),
						"\n")

	if(length(acc) == 0){
		stop("All posterior samples have been rejected. Please rerun getSmoothed() with a larger 'nsim' argument.")
	} 

	# Update fitted
	for(i in 1:length(fitted$draws.est )){
		fitted$draws.est[[i]]$draws <- fitted$draws.est[[i]]$draws[acc]  
	}
	for(i in 1:length(fitted$draws.est.overall )){
		fitted$draws.est.overall[[i]]$draws <- fitted$draws.est.overall[[i]]$draws[acc]  
	}


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

	fitted$overall$variance <- fitted$overall$median <- fitted$overall$mean <-  fitted$overall$lower <-  fitted$overall$upper <- NA
	for(i in 1:length(fitted$draws.est.overall)){
		y <- fitted$draws.est.overall[[i]]$years
		r <- fitted$draws.est.overall[[i]]$region
		j <- which(fitted$overall$years == y & fitted$overall$region == r)
		fitted$overall[j, c("lower", "median", "upper")] <- stats::quantile(fitted$draws.est.overall[[j]]$draws, c(lowerCI, 0.5, upperCI))
        fitted$overall[j, "mean"] <- mean(fitted$draws.est.overall[[j]]$draws)
        fitted$overall[j, "variance"] <- var(fitted$draws.est.overall[[j]]$draws)
	}

	fitted$nsim <- length(acc)
	fitted$accept.ratio <- length(acc) / n0
	fitted$benchmarked <- TRUE

	return(fitted)
}

