#' Deprecated function for smoothed direct estimates for mortality rates 
#' 
#' This is an wrapper function for smoothSurvey function for the purpose of backward compatibility. 
#' 
#' @param ... arguments of the new \code{\link{smoothSurvey}} function
#' @seealso  \code{\link{smoothSurvey}}
#' @examples
#' \dontrun{
#' data(DemoData2)
#' data(DemoMap2)
#' fit <- smoothSurvey(data=DemoData2,  
#' Amat=DemoMap2$Amat, responseType="binary", 
#' responseVar="tobacco.use", strataVar="strata", 
#' weightVar="weights", regionVar="region", 
#' clusterVar = "~clustid+id", CI = 0.95)
#' 
#' # Example with region-level covariates
#'  Xmat <- aggregate(age~region, data = DemoData2, FUN = mean)
#'  fit <- smoothSurvey(data=DemoData2,  
#'   Amat=DemoMap2$Amat, responseType="binary", 
#'   X = Xmat,
#'   responseVar="tobacco.use", strataVar="strata", 
#'   weightVar="weights", regionVar="region", 
#'   clusterVar = "~clustid+id", CI = 0.95)
#' }
#' @export
fitGeneric <- function(...){
	smoothSurvey(...)
}


#' Deprecated function for smoothed direct estimates for mortality rates 
#' 
#' This is an wrapper function for smoothDirect function for the purpose of backward compatibility. 
#' 
#' @param ... arguments of the new \code{\link{smoothDirect}} function
#' @seealso  \code{\link{smoothDirect}}
#' @examples
#' \dontrun{
#' years <- levels(DemoData[[1]]$time)
#' # obtain direct estimates
#' data_multi <- getDirectList(births = DemoData, years = years,
#'   regionVar = "region",  timeVar = "time", clusterVar = "~clustid+id",
#'   ageVar = "age", weightsVar = "weights", geo.recode = NULL)
#' data <- aggregateSurvey(data_multi)
#' 
#' #  national model
#' years.all <- c(years, "15-19")
#' fit1 <- smoothDirect(data = data, Amat = NULL, 
#'   year_label = years.all, year_range = c(1985, 2019), 
#'   rw = 2, is.yearly=FALSE, m = 5)
#' out1 <- getSmoothed(fit1)
#' plot(out1, is.subnational=FALSE)
#' 
#' #  subnational model
#' fit2 <- smoothDirect(data = data, Amat = DemoMap$Amat, 
#'   year_label = years.all, year_range = c(1985, 2019), 
#'   rw = 2, is.yearly=TRUE, m = 5, type.st = 4)
#' out2 <- getSmoothed(fit2, Amat = mat)
#' plot(out2, is.yearly=TRUE, is.subnational=TRUE)
#' }
#' @export
fitINLA <- function(...){
	smoothDirect(...)
}



#' Deprecated function for cluster-level space-time smoothing models for mortality rates 
#' 
#' This is an wrapper function for smoothCluster function for the purpose of backward compatibility. 
#' 
#' @param ... arguments of the new \code{\link{smoothCluster}} function
#' @seealso  \code{\link{smoothCluster}}
#' @examples
#' \dontrun{
#' library(dplyr)
#' data(DemoData)
#' # Create dataset of counts
#' counts.all <- NULL
#' for(i in 1:length(DemoData)){
#'   counts <- getCounts(DemoData[[i]][, c("clustid", "time", "age", "died",
#'                                         "region", "strata")],
#'             variables = 'died', by = c("age", "clustid", "region", 
#'                                          "time", "strata"))
#'   counts <- counts %>% mutate(cluster = clustid, years = time, Y=died)
#'   counts$strata <- gsub(".*\\.","",counts$strata)
#'   counts$survey <- names(DemoData)[i] 
#'   counts.all <- rbind(counts.all, counts)
#' }
#' 
#' # fit cluster-level model on the periods
#' periods <- levels(DemoData[[1]]$time)
#' fit <- smoothCluster(data = counts.all, 
#'       Amat = DemoMap$Amat, 
#'       time.model = "rw2", 
#'       st.time.model = "rw1",
#'       strata.time.effect =  TRUE, 
#'       survey.effect = TRUE,
#'       family = "betabinomial",
#'       year_label = c(periods, "15-19"))
#' est <- getSmoothed(fit, Amat = DemoMap$Amat, 
#'       year_label =c(periods, "15-19"), nsim = 1000)
#' plot(est$stratified, plot.CI=TRUE) + ggplot2::facet_wrap(~strata) 
#' }
#' @export
fitINLA2 <- function(...){
	smoothCluster(...)
}