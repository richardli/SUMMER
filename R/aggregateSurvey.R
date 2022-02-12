#' Aggregate estimators from different surveys. 
#' 
#' @param data Output from \code{\link{getDirectList}}
#' 
#' @return Estimators aggregated across surveys.
#' @author Zehang Richard Li
#' @examples
#' \dontrun{
#' data(DemoData)
#' data(DemoMap)
#' years <- levels(DemoData[[1]]$time)
#' 
#' # obtain direct estimates
#' data <- getDirectList(births = DemoData, 
#' years = years, 
#' regionVar = "region", timeVar = "time", 
#' clusterVar = "~clustid+id", 
#' ageVar = "age", weightsVar = "weights", 
#' geo.recode = NULL)
#' 
#' # obtain maps
#' geo <- DemoMap$geo
#' mat <- DemoMap$Amat
#' 
#' # Simulate hyper priors
#' priors <- simhyper(R = 2, nsamp = 1e+05, nsamp.check = 5000, Amat = mat, only.iid = TRUE)
#' 
#' # combine data from multiple surveys
#' data <- aggregateSurvey(data)
#' utils::head(data)
#' 
#' }
#' 
#' @export
aggregateSurvey <- function(data) {
  data0 <- data
  data0$logit.prec <- 1/data0$var.est
  time_region <- unique(data0[, c("region", "years")])
  
  data <- data.frame(region = time_region$region, years = time_region$years, mean = NA, lower=NA, upper=NA, logit.est=NA, var.est=NA, region_num = NA, survey = NA, logit.prec = NA)
  for(i in 1:dim(data)[1]){
    tmp <- intersect(which(data0$region == data$region[i]),
                     which(data0$years == data$years[i]))
    # Version adjusting for HIV
    data[i, "logit.prec"] <- sum(data0[tmp, "logit.prec"], na.rm = TRUE)
    if(data[i, "logit.prec"] == 0){
      data[i, "var.est"] <- NA
      data[i, "logit.prec"] <- NA
    }else{
      data[i, "var.est"] <- 1 / data[i, "logit.prec"]
      weights <- data0[tmp, "logit.prec"] / data[i, "logit.prec"]
      data[i, "logit.est"] <- sum(weights * data0[tmp, "logit.est"], na.rm = TRUE)
      data[i, "mean"] <- expit(data[i, "logit.est"])
      
      data[i, "lower"] <- expit(data[i, "logit.est"] + stats::qnorm(0.025)*sqrt(data[i, "var.est"]))
      data[i, "upper"] <- expit(data[i, "logit.est"] + stats::qnorm(0.975)*sqrt(data[i, "var.est"]))
    }
  }
  return(data)
}
