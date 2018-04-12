#' Plot projection output.
#' @param x output from \code{\link{projINLA}}
#' @param years_label labels for the periods
#' @param years_med labels for the middle years in each period
#' @param is.yearly logical indicator of whether the data contains yearly estimates
#' @param is.subnational logical indicator of whether the data contains subnational estimates
#' @param proj_year The first year where projections are made, i.e., where no data are available. 
#' @param ... optional arguments, see details
#' 
#' @details 
#' Note that arguments after \code{...} must match exactly.
#' \itemize{
#'  \item{\code{years_label}}{string of year labels, defaults to \code{c("85-89", "90-94", "95-99", "00-04", "05-09", "10-14", "15-19")}}
#'  \item{\code{proj_year}}{projection year as numeric, defaults to \code{2015}}
#'  \item{\code{years_med}}{ median of year intervals, defaults to \code{c(1987, 1992, 1997, 2002, 2007, 2012, 2017)}}
#'  \item{\code{is.yearly}}{indicator for yearly model, defaults to \code{TRUE}}
#'  \item{\code{is.subnational}}{indicator for subnational model, defaults to \code{TRUE}}
#' }
#' @examples
#' \dontrun{
#' data(DemoData)
#' deta(DemoMap)
#' years <- levels(DemoData[[1]]$time)
#' 
#' # obtain direct estimates
#' data <- countrySummary_mult(births = DemoData, 
#' years = years, idVar = "id", 
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
#' out <- projINLA(fit, Amat = mat, is.yearly = TRUE)
#' plot(out, is.yearly=TRUE, is.subnational=TRUE) + ggplot2::ggtitle("Subnational yearly model")
#' 
#' }
#' @export
plot.projINLA <- function(x, years_label = c("85-89", "90-94", "95-99", "00-04", "05-09", "10-14", "15-19"), years_med = c(1987, 1992, 1997, 2002, 2007, 2012, 2017), is.yearly=TRUE, is.subnational = TRUE, proj_year = 2015, ...){

  if(is.null(proj_year)) {
    proj_year = 0
  }
  
  is.periods <- x$Year %in% years_label
  x$Year.num[is.periods] <- years_med[match(x$Year[is.periods], years_label)]
  x$project <- FALSE
  x$project[x$Year.num >= proj_year] <- TRUE
  # fix for global variable issue
  Year.num <- NULL; District <- NULL; med <- NULL; q025 <- NULL; q975 <- NULL; project <- NULL
  
  if(is.subnational){
    g <- ggplot2::ggplot(ggplot2::aes(x = Year.num, y = med, ymin = q025, ymax = q975, color = District), data = x)
    my.dodge <- ggplot2::position_dodge(width = 1)
  }else{
    g <- ggplot2::ggplot(ggplot2::aes(x = Year.num, y = med, ymin = q025, ymax = q975), data = x)
    my.dodge <- ggplot2::position_dodge(width = 0.2)
  }
  
  if(!is.yearly){
    g <- g + ggplot2::geom_point(position = my.dodge)
    g <- g + ggplot2::geom_line(position = my.dodge)
    g <- g + ggplot2::geom_errorbar(ggplot2::aes(linetype=project), size = .7, width = .05, position = my.dodge)
    g <- g + ggplot2::theme_bw() + ggplot2::xlab("Year") + ggplot2::ylab("U5MR")
    g <- g + ggplot2::scale_x_continuous(breaks=years_med, labels=years_label)
  }else if(!is.subnational){
    g <- g + ggplot2::geom_point(position = my.dodge, data=subset(x, is.periods==FALSE), alpha = 0.5, color = 1)
    g <- g + ggplot2::geom_line(position = my.dodge, data=subset(x, is.periods==FALSE), alpha = 0.5, color = 1)
    g <- g + ggplot2::geom_errorbar(ggplot2::aes(linetype=project), size = .5, width = .05, position = my.dodge, data=subset(x, is.periods==FALSE), alpha = 0.3, color = 1)
    g <- g + ggplot2::geom_point(shape = 17, size = 2.5, position = my.dodge, data=subset(x, is.periods==TRUE), color = 2)
    g <- g + ggplot2::geom_errorbar(ggplot2::aes(linetype=project), size = .7, width = .05, position = my.dodge, data=subset(x, is.periods==TRUE), color = 2)
    g <- g + ggplot2::theme_bw() + ggplot2::xlab("Year") + ggplot2::ylab("U5MR")
  }else if(is.subnational){
    g <- g + ggplot2::geom_point(position = my.dodge, data=subset(x, is.periods==FALSE), alpha = 0.5)
    g <- g + ggplot2::geom_line(position = my.dodge, data=subset(x, is.periods==FALSE), alpha = 0.5)
    g <- g + ggplot2::geom_point(shape = 17, size = 2.5, position = my.dodge, data=subset(x, is.periods==TRUE))
    g <- g + ggplot2::geom_errorbar(ggplot2::aes(linetype=project), size = .7, width = .05, position = my.dodge, data=subset(x, is.periods==TRUE))
    g <- g + ggplot2::theme_bw() + ggplot2::xlab("Year") + ggplot2::ylab("U5MR")
  }
  g
}