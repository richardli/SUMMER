#' Summary method for the smoothing models.
#' 
#' This function is the summary method for class \code{SUMMERmodel}.
#' 
#' 
#' @param object output from \code{\link{smoothDirect}} or \code{\link{smoothCluster}}
#' @param ... not used
#' 
#' @author Zehang Li 
#' 
#' @seealso \code{\link{summary.SUMMERmodel}} 
#' @method summary SUMMERmodel
#' @examples
#' \dontrun{
#'   library(SUMMER)
#'   library(dplyr)
#'   data(DemoData)
#' 
#'   # Smooth Direct Model
#'   years <- levels(DemoData[[1]]$time)
#'   # obtain direct estimates
#'   data_multi <- getDirectList(births = DemoData, years = years,
#'   regionVar = "region",  timeVar = "time", clusterVar = "~clustid+id",
#'   ageVar = "age", weightsVar = "weights", geo.recode = NULL)
#'   data <- aggregateSurvey(data_multi)
#'   
#'   years.all <- c(years, "15-19")
#'   fit <- smoothDirect(data = data, Amat = NULL, 
#'   year_label = years.all, year_range = c(1985, 2019), 
#'   time.model = 'rw2', is.yearly=FALSE, m = 5)
#'   summary(fit)
#' 
#'   # Cluster-level Model
#'   counts.all <- NULL
#'   for(i in 1:length(DemoData)){
#'   counts <- getCounts(DemoData[[i]][, c("clustid", "time", "age", "died",
#'                                        "region", "strata")],
#'            variables = 'died', by = c("age", "clustid", "region", 
#'                                         "time", "strata"))
#'   counts <- counts %>% mutate(cluster = clustid, years = time, Y=died)
#'   counts$strata <- gsub(".*\\.","",counts$strata)
#'   counts$survey <- names(DemoData)[i] 
#'   counts.all <- rbind(counts.all, counts)
#'   }
#'   
#'   # fit cluster-level model on the periods
#'   periods <- levels(DemoData[[1]]$time)
#'   fit <- smoothCluster(data = counts.all, 
#'      Amat = DemoMap$Amat, 
#'      time.model = "rw2", 
#'      st.time.model = "rw1",
#'      strata.time.effect =  TRUE, 
#'      survey.effect = TRUE,
#'      family = "betabinomial",
#'      year_label = c(periods, "15-19"))
#'   summary(fit) 
#' }
#' @export 

summary.SUMMERmodel <- function(object,...){
	cat("----------------------------------")
	config <- object$msg
	cat(config)  
	cat("----------------------------------\n")
	cat("Fixed Effects\n")
	fixed <- summary(object$fit)$fixed
	print(fixed)
	cat("----------------------------------\n")
	cat("Random Effects\n")
	random <- data.frame(Name = summary(object$fit)$random.names,
						 Model = summary(object$fit)$random.model) 
	print(random)
	cat("----------------------------------\n")
	cat("Model hyperparameters\n")
	hyperpar <- summary(object$fit)$hyperpar
	print(hyperpar)

	neffp <- summary(object$fit)$neffp
	print(neffp)
	mlik <- summary(object$fit)$mlik
	print(mlik)
}
 

#' Print method for the smoothing models.
#' 
#' This function is the print method for class \code{SUMMERmodel}.
#' 
#' 
#' @param x output from \code{\link{smoothDirect}} or \code{\link{smoothCluster}}
#' @param ... not used
#' @method print SUMMERmodel
#' @author Zehang Li 
#' 
#' @seealso \code{\link{summary.SUMMERmodel}} 
#' 
#' @examples
#' \dontrun{
#'   library(SUMMER)
#'   library(dplyr)
#'   data(DemoData)
#' 
#'   # Smooth Direct Model
#'   years <- levels(DemoData[[1]]$time)
#'   # obtain direct estimates
#'   data_multi <- getDirectList(births = DemoData, years = years,
#'   regionVar = "region",  timeVar = "time", clusterVar = "~clustid+id",
#'   ageVar = "age", weightsVar = "weights", geo.recode = NULL)
#'   data <- aggregateSurvey(data_multi)
#'   
#'   years.all <- c(years, "15-19")
#'   fit <- smoothDirect(data = data, Amat = NULL, 
#'   year_label = years.all, year_range = c(1985, 2019), 
#'   time.model = 'rw2', is.yearly=FALSE, m = 5)
#'   fit
#' 
#'   # Cluster-level Model
#'   counts.all <- NULL
#'   for(i in 1:length(DemoData)){
#'   counts <- getCounts(DemoData[[i]][, c("clustid", "time", "age", "died",
#'                                        "region", "strata")],
#'            variables = 'died', by = c("age", "clustid", "region", 
#'                                         "time", "strata"))
#'   counts <- counts %>% mutate(cluster = clustid, years = time, Y=died)
#'   counts$strata <- gsub(".*\\.","",counts$strata)
#'   counts$survey <- names(DemoData)[i] 
#'   counts.all <- rbind(counts.all, counts)
#'   }
#'   
#'   # fit cluster-level model on the periods
#'   periods <- levels(DemoData[[1]]$time)
#'   fit <- smoothCluster(data = counts.all, 
#'      Amat = DemoMap$Amat, 
#'      time.model = "rw2", 
#'      st.time.model = "rw1",
#'      strata.time.effect =  TRUE, 
#'      survey.effect = TRUE,
#'      family = "betabinomial",
#'      year_label = c(periods, "15-19"))
#'   fit
#' }
#' @export 

print.SUMMERmodel <- function(x,...){
	cat("----------------------------------")
	cat(x$msg)  
	cat("----------------------------------\n INLA ")
	print(x$fit)
}


#' Print method for the combined projection output.
#' 
#' This function is the print method for class \code{SUMMERprojlist}.
#' 
#' 
#' @param x output from \code{\link{getSmoothed}}
#' @param ... not used
#' @method print SUMMERprojlist
#' @author Zehang Li 
#' 
#' @examples
#' \dontrun{
#'  library(SUMMER)
#'  library(dplyr)
#'  data(DemoData)
#'  # Create dataset of counts
#'  counts.all <- NULL
#'  for(i in 1:length(DemoData)){
#'  counts <- getCounts(DemoData[[i]][, c("clustid", "time", "age", "died",
#'                                       "region", "strata")],
#'           variables = 'died', by = c("age", "clustid", "region", 
#'                                        "time", "strata"))
#'  counts <- counts %>% mutate(cluster = clustid, years = time, Y=died)
#'  counts$strata <- gsub(".*\\.","",counts$strata)
#'  counts$survey <- names(DemoData)[i] 
#'  counts.all <- rbind(counts.all, counts)
#'  }
#'  
#'  # fit cluster-level model on the periods
#'  periods <- levels(DemoData[[1]]$time)
#'  fit <- smoothCluster(data = counts.all, 
#'     Amat = DemoMap$Amat, 
#'     time.model = "rw2", 
#'     st.time.model = "rw1",
#'     strata.time.effect =  TRUE, 
#'     survey.effect = TRUE,
#'     family = "betabinomial",
#'     year_label = c(periods, "15-19"))
#'  summary(fit)
#'  est <- getSmoothed(fit, nsim = 1000)
#' }
#' @export 

print.SUMMERprojlist <- function(x, ...){
	cat("---------------------------------------------\n")
	cat("Stratified estimates stored in ...$stratified\n")
	cat("Aggregated estimates stored in ...$overall\n")
	cat("---------------------------------------------\n")
	cat(paste0("Estimates computed for ", max(x$overall$time), " time period(s) and ", max(x$overall$area), " area(s)"))
	cat(x$msg)  
	cat("\n")
	cat(paste0(x$nsim, " posterior draws taken.\n"))
}


#' Summary method for the smoothing model and output from \code{smoothSurvey}.
#' 
#' This function is the summary method for class \code{SUMMERmodel.svy}.
#' 
#' 
#' @param object output from \code{\link{smoothSurvey}} 
#' @param ... not used
#' 
#' @author Zehang Li 
#' 
#' @seealso \code{\link{summary.SUMMERmodel.svy}} 
#' @method summary SUMMERmodel.svy
#' @examples
#' \dontrun{
#' data(DemoData2)
#' data(DemoMap2)
#' fit0 <- smoothSurvey(data=DemoData2,  
#' Amat=DemoMap2$Amat, responseType="binary", 
#' responseVar="tobacco.use", strataVar="strata", 
#' weightVar="weights", regionVar="region", 
#' clusterVar = "~clustid+id", CI = 0.95)
#' summary(fit0)
#' }
#' @export 

summary.SUMMERmodel.svy <- function(object,...){
	cat("----------------------------------\n")
	config <- object$msg
	cat(config)  
	cat("\n----------------------------------\n")
	cat("Fixed Effects\n")
	fixed <- summary(object$fit)$fixed
	print(fixed)
	cat("----------------------------------\n")
	cat("Random Effects\n")
	random <- data.frame(Name = summary(object$fit)$random.names,
						 Model = summary(object$fit)$random.model) 
	print(random)
	cat("----------------------------------\n")
	cat("Model hyperparameters\n")
	hyperpar <- summary(object$fit)$hyperpar
	print(hyperpar)

	neffp <- summary(object$fit)$neffp
	print(neffp)
	mlik <- summary(object$fit)$mlik
	print(mlik)
}
 

#' Print method for the smoothing models from \code{smoothSurvey}.
#' 
#' This function is the print method for class \code{SUMMERmodel.svy}.
#' 
#' 
#' @param x output from \code{\link{smoothSurvey}}.
#' @param ... not used
#' @method print SUMMERmodel.svy
#' @author Zehang Li 
#' 
#' @seealso \code{\link{summary.SUMMERmodel.svy}} 
#' 
#' @examples
#' \dontrun{
#' data(DemoData2)
#' data(DemoMap2)
#' fit0 <- smoothSurvey(data=DemoData2,  
#' Amat=DemoMap2$Amat, responseType="binary", 
#' responseVar="tobacco.use", strataVar="strata", 
#' weightVar="weights", regionVar="region", 
#' clusterVar = "~clustid+id", CI = 0.95)
#' fit0
#' }
#' @export 

print.SUMMERmodel.svy <- function(x,...){
	cat("----------------------------------\n")
	cat(x$msg)  
	cat("\n\n$formula\n")
	print(paste(as.character(x$formula)[c(2,1,3)], collapse = " "))
	if(!is.null(x$HT)){
		cat("----------------------------------\n")
		cat("$HT\n")
		print(head(x$HT))
		cat("...\n")
	}
	if(!is.null(x$smooth)){
		cat("----------------------------------\n")
		cat("$smooth\n")
		print(head(x$smooth))
		cat("...\n")
	}
	if(!is.null(x$smooth.overall)){
		cat("----------------------------------\n")
		cat("$smooth\n")
		print(head(x$smooth.overall))
		cat("...\n")
	}

}
