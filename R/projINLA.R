#' Function to obtain projected estimates from INLA for each time and region.
#' 
#' 
#'
#' @param inla_mod output from \code{\link{fitINLA}}
#' @param year_range range corresponding to year label
#' @param year_label vector of year string vector
#' @param Amat adjacency matrix
#' @param nsim number of simulations
#' @param ... not used
#' 
#' @return Results from RW2 model fit, including projection.
#' 
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
#' 
#' # Model fitting with INLA
#' years.all <- c(years, "15-19")
#' fit <- fitINLA(data = data, geo = geo, Amat = mat, 
#' year_names = years.all, year_range = c(1985, 2019), 
#' priors = priors, rw = 2, is.yearly=TRUE, 
#' m = 5, type.st = 4)
#' # Projection
#' out <- getSmoothed(fit, Amat = mat)
#' plot(out, is.subnational=TRUE) + ggplot2::ggtitle("Subnational yearly model")
#' 
#' }
#' 
#' 
#' @export
getSmoothed <- function(inla_mod, year_range = c(1985, 2019), year_label = c("85-89", "90-94", "95-99", "00-04", "05-09", "10-14", "15-19"), 
                            Amat = NULL, nsim = 1000, ...){

  if (!isTRUE(requireNamespace("INLA", quietly = TRUE))) {
    stop("You need to install the packages 'INLA'. Please run in your R terminal:\n install.packages('INLA', repos='https://www.math.ntnu.no/inla/R/stable')")
  }
  # If INLA is installed, then attach the Namespace (so that all the relevant functions are available)
  if (isTRUE(requireNamespace("INLA", quietly = TRUE))) {
    if (!is.element("INLA", (.packages()))) {
      attachNamespace("INLA")
    }
  if(is.null(Amat)){
    region_names <- "All"
    region_nums <- 0
  }else{
    region_names <- colnames(Amat)
    region_nums <- 1:length(region_names)
  }
  is.yearly = inla_mod$is.yearly
  if(is.yearly){
    timelabel.yearly <- c(year_range[1] : year_range[2], year_label)
  }else{
    timelabel.yearly <- year_label
  }
  results <- expand.grid(District = region_nums, Year = timelabel.yearly)
  results$median <- results$lower <- results$upper <- results$logit.median <- results$logit.lower <- results$logit.upper <- NA
  mod <- inla_mod$fit
  lincombs.info <- inla_mod$lincombs.info

  for(i in 1:length(timelabel.yearly)){
    for(j in 1:length(region_names)){
        index <- lincombs.info$Index[lincombs.info$District == region_nums[j] & lincombs.info$Year == i]
        tmp.logit <- INLA::inla.rmarginal(nsim, mod$marginals.lincomb.derived[[index]])
        marg <- INLA::inla.tmarginal(expit, mod$marginals.lincomb.derived[[index]])
        tmp <- INLA::inla.rmarginal(nsim, marg)

        results$median[results$District == region_nums[j] & results$Year == timelabel.yearly[i]] <- stats::median(tmp)
        results$upper[results$District == region_nums[j] & results$Year == timelabel.yearly[i]] <- stats::quantile(tmp, .975)
        results$lower[results$District == region_nums[j] & results$Year == timelabel.yearly[i]] <- stats::quantile(tmp, .025)
        results$logit.median[results$District == region_nums[j] & results$Year == timelabel.yearly[i]] <- stats::median(tmp.logit)
        results$logit.upper[results$District == region_nums[j] & results$Year == timelabel.yearly[i]] <- stats::quantile(tmp.logit, .975)
        results$logit.lower[results$District == region_nums[j] & results$Year == timelabel.yearly[i]] <- stats::quantile(tmp.logit, .025)

    }
  }
  results$is.yearly <- !(results$Year %in% year_label)
  results$Year.num <- suppressWarnings(as.numeric(as.character(results$Year)))
  if(region_names[1] != "All"){
    results$District <- region_names[results$District]
  }else{
    results$District <- "All"
  }
  colnames(results)[which(colnames(results) == "District")] <- "region"
  colnames(results)[which(colnames(results) == "Year")] <- "years"
  # Add S3 method
  class(results) <- c("SUMMERproj", "data.frame")
  return(results)
  }
}


