#' Adjust direct estimates and their associated variances
#' 
#' @param data data frame of the adjusted estimates and the associated uncertainties, see the arguments below for specific columns.
#' @param ratio the ratio of unadjusted mortality rates to the true mortality rates. It can be either a data frame with the following three columns (region, time, and ratio) if adjustment factor differ by region; or a data frame with the following two columns (time and ratio) if adjustment factor only varies over time.
#' @param time the column name for time in the data 
#' @param region the column name for region in the data 
#' @param est the column name for unadjusted mortality rates in the data 
#' @param logit the column name for the logit of the unadjusted mortality rates in the data 
#' @param logit.var the column name for the variance of the logit of the unadjusted mortality rates in the data
#' @param logit.prec the column name for the precision of the logit of the unadjusted mortality rates in the data
#' @param logit.lower the column name for the 95\% lower bound of the logit of the unadjusted mortality rates in the data
#' @param logit.upper the column name for the 95\% lower bound of the logit of the unadjusted mortality rates in the data
#' 
#' @return adjusted dataset of the same columns.
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
#' adj <- expand.grid(region = unique(data$region), time = years)
#' adj$ratio <- runif(dim(adj)[1], min = 0.5, max = 0.8)
#' data.adj <- getAdjusted(data = data, ratio = adj)
#'  }
#' @export
#' 

getAdjusted <- function(data, ratio, time = "years", region = "region",  est = "mean", logit = "logit.est", logit.var = "var.est", logit.prec = "logit.prec", logit.lower = "lower", logit.upper = "upper"){

	adjust <- function(p, v, c){
		f.prime <- 1 - (c - 1) * (p/(1-p)) / (c + (c-1) * (p/(1-p)))
		p <- p / c
		v <- v * f.prime^2
		return(c(p, log(p/(1-p)), v, 1/v))
	}

	count <- 0
	for(i in 1:dim(data)[1]){
		t <- data[i, time]
		s <- data[i, region]
		if("region" %in% colnames(ratio)){
			matched <- intersect(which(ratio$region == s), which(ratio$time == t))
		}else{
			matched <- which(ratio$time == t)
		}
		if(length(matched) == 0){
			next
		}else if(length(matched) == 1){
			count <- count + 1
			new <- adjust(data[i, est], data[i, logit.var], ratio[matched, "ratio"])
			data[i, est] <- new[1]
			data[i, logit] <- new[2]
			data[i, logit.var] <- new[3]
			data[i, logit.prec] <- new[4]
			data[i, logit.lower] <- new[2] - stats::qnorm(0.025)*sqrt(new[3])
			data[i, logit.upper] <- new[2] + stats::qnorm(0.025)*sqrt(new[3])

		}else{
			stop("Multiple adjustment factors found for some estimates.")
		}
	}
	message(paste0(count, " estimates were adjusted."))
	return(data)
}

