#' Adjust direct estimates and their associated variances
#' 
#' @param data data frame of the adjusted estimates and the associated uncertainties, see the arguments below for specific columns.
#' @param ratio the ratio of unadjusted mortality rates to the true mortality rates. It can be either a data frame with the following three columns (region, time, and adj) if adjustment factor differ by region; or a data frame with the following two columns (time and adj) if adjustment factor only varies over time. The column names specifying region, time, and adjustment are specified by the arguments in the function call.
#' @param time the column name for time in the data and adjustment ratio.
#' @param region the column name for region in the data  and adjustment ratio.
#' @param est the column name for unadjusted mortality rates in the data 
#' @param logit the column name for the logit of the unadjusted mortality rates in the data 
#' @param logit.var the column name for the variance of the logit of the unadjusted mortality rates in the data
#' @param logit.prec the column name for the precision of the logit of the unadjusted mortality rates in the data
#' @param logit.lower the column name for the 95\% lower bound of the logit of the unadjusted mortality rates in the data
#' @param logit.upper the column name for the 95\% lower bound of the logit of the unadjusted mortality rates in the data
#' @param prob.lower the column name for the 95\% lower bound of the unadjusted mortality rates in the data. If this is provided instead of logit.lower, the logit scale lower bound will be created.
#' @param prob.upper the column name for the 95\% lower bound of the unadjusted mortality rates in the data. if this is provided instead of logit.upper, the logit scale upper bound will be created.
#' @param adj the column name for the adjustment ratio
#' @param lower previous argument name for prob.lower. Will be removed in the next update
#' @param upper previous argument name for prob.upper. Will be removed in the next update
#' @param verbose logical indicator for whether to print out unadjusted row index
#' 
#' @return adjusted dataset of the same columns.
#' @author Zehang Richard Li 
#' @examples
#' \dontrun{
#' years <- levels(DemoData[[1]]$time)
#' 
#' # obtain direct estimates
#' data <- getDirectList(births = DemoData, 
#' years = years,  
#' regionVar = "region", timeVar = "time", 
#' clusterVar = "~clustid+id", 
#' ageVar = "age", weightsVar = "weights", 
#' geo.recode = NULL)
#' # obtain direct estimates
#' data_multi <- getDirectList(births = DemoData, years = years,
#'   regionVar = "region",  timeVar = "time", clusterVar = "~clustid+id",
#'   ageVar = "age", weightsVar = "weights", geo.recode = NULL)
#' data <- aggregateSurvey(data_multi)
#' 
#' # randomly simulate adjustment factor
#' adj <- expand.grid(region = unique(data$region), years = years)
#' adj$ratio <- runif(dim(adj)[1], min = 0.5, max = 0.8)
#' data.adj <- getAdjusted(data = data, ratio = adj)
#'  }
#' @export
#' 

getAdjusted <- function(data, ratio, time = "years", region = "region",  est = "mean", logit = "logit.est", logit.var = "var.est", logit.prec = "logit.prec", logit.lower = "lower", logit.upper = "upper", prob.lower = NULL, prob.upper = NULL, adj = "ratio", verbose = FALSE, lower = NULL, upper = NULL){

	adjust <- function(p, v, c){
		f.prime <- 1 - (c - 1) * (p/(1-p)) / (c + (c-1) * (p/(1-p)))
		p <- p / c
		v <- v * f.prime^2
		return(c(p, log(p/(1-p)), v, 1/v))
	}
	if(!is.null(lower)){
		if(is.null(prob.lower)){
			message("Arguments 'lower' and 'upper' has been renamed to 'prob.lower' and 'prob.upper' to avoid confusion. 'lower' and 'upper' will be removed in the next major update.")
			prob.lower <- lower
			prob.upper <- upper			
		}else{
			message("Arguments 'lower' and 'upper' has been renamed to 'prob.lower' and 'prob.upper' to avoid confusion. 'lower' and 'upper' arguments will be ignored.")
		}
	}

	if(identical(prob.lower, logit.lower)){
		stop("The same column has been specified for logit.lower and prob.lower. Please specify only one of them.")
	}
	if(identical(prob.upper, logit.upper)){
		stop("The same column has been specified for logit.upper and prob.upper. Please specify only one of them.")
	}

	if(!is.null(prob.lower) && (is.null(logit.lower) || is.na(logit.lower))){
		data$logit.lower <- logit(data[, prob.lower])
		data$logit.upper <- logit(data[, prob.upper])
		logit.lower <- "logit.lower"
		logit.upper <- "logit.upper"
	}

	## Check all regions are adjusted for each time period
	if(region %in% colnames(ratio)){
		region.names <- unique(data[, region])
		time.names <- unique(ratio[, time])
		time.names <- time.names[time.names %in% data[, time]]
		for(t in time.names){
			sub <- ratio[ratio[, time] == t, ]
			all.there <- sum(region.names %in% sub[, region]) == length(region.names)
			if(!all.there){
				stop("For some time periods, there are missing adjustment ratios for some regions. Check that all regions have adjustment ratios in each time period, including national estimates (region = 'All'), if they exist in the data.")
			}
		}
	}

	count <- 0
	warn_national <- FALSE
	unadj <- NULL
	for(i in 1:dim(data)[1]){
		t <- data[i, time]
		s <- data[i, region]
		if(region %in% colnames(ratio)){
			matched <- intersect(which(ratio[, region] == s), which(ratio[, time] == t))
		}else{
			matched <- which(ratio[, time] == t)
		}
		if(length(matched) == 0){
			warn_national <- (s == "All") && (region %in% colnames(ratio))
			unadj <- c(unadj, i)
			next
		}else if(length(matched) == 1){
			count <- count + 1
			new <- adjust(data[i, est], data[i, logit.var], ratio[matched, adj])
			data[i, est] <- new[1]
			data[i, logit] <- new[2]
			data[i, logit.var] <- new[3]
			data[i, logit.prec] <- new[4]
			data[i, logit.lower] <- new[2] + stats::qnorm(0.025)*sqrt(new[3])
			data[i, logit.upper] <- new[2] - stats::qnorm(0.025)*sqrt(new[3])

		}else{
			stop("Multiple adjustment factors found for some estimates.")
		}
	}
	message(paste0(count, " estimates were adjusted."))
	if(warn_national){
		warning("National estimates are not adjusted. If you will benchmark smoothed direct estimates,  make sure you have national adjustment factors as well. ")
	}
	if(!is.null(unadj) && verbose){
		warning(paste0("The following rows of the data are not adjusted: ", paste(unadj, collapse = ", ")))
	}
	if(!is.null(prob.lower)){
		data[, prob.lower] <- expit(data[, logit.lower])
		data[, prob.upper] <- expit(data[, logit.upper])
	}
	return(data)
}

