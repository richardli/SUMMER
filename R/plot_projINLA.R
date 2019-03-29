#' Plot projection output.
#' @param x output from \code{\link{projINLA}}
#' @param year_label labels for the periods
#' @param year_med labels for the middle years in each period
#' @param is.yearly logical indicator of whether the data contains yearly estimates
#' @param is.subnational logical indicator of whether the data contains subnational estimates
#' @param proj_year the first year where projections are made, i.e., where no data are available. 
#' @param data.add data frame for the Comparisons data points to add to the graph. This can be, for example, the raw direct estimates. This data frame is merged to the projections by column 'region' and 'years'. Except for these two columns, this dataset should not have Comparisons columns with names overlapping the projINLA output.
#' @param option.add list of options specifying the variable names for the points to plot, lower and upper bounds, and the grouping variable. This is intended to be used to add Comparisons estimates on the same plot as the smoothed estimates. See examples for details.  
#' @param color.add the color of the Comparisons data points to plot.
#' @param dodge.width the amount to add to data points at the same year to avoid overlap. Default to be 1.
#' @param ... optional arguments, see details
#' 
#' @details 
#' Note that arguments after \code{...} must match exactly.
#' \itemize{
#'  \item{\code{year_label}}{string of year labels, defaults to \code{c("85-89", "90-94", "95-99", "00-04", "05-09", "10-14", "15-19")}}
#'  \item{\code{proj_year}}{projection year as numeric, defaults to \code{2015}}
#'  \item{\code{year_med}}{ median of year intervals, defaults to \code{c(1987, 1992, 1997, 2002, 2007, 2012, 2017)}}
#'  \item{\code{is.yearly}}{indicator for yearly model, defaults to \code{TRUE}}
#'  \item{\code{is.subnational}}{indicator for subnational model, defaults to \code{TRUE}}
#' }
#' @examples
#' \dontrun{
#' data(DemoData)
#' deta(DemoMap)
#' years <- levels(DemoData[[1]]$time)
#' 
#' data <- countrySummary_mult(births = DemoData, 
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
#' # combine data from multiple surveys
#' data_agg <- aggregateSurvey(data)
#' 
#' # Model fitting with INLA
#' years.all <- c(years, "15-19")
#'
#' 
#' fit <- fitINLA(data = data_agg, geo = NULL, Amat = NULL, 
#' year_names = years.all, year_range = c(1985, 2019),  
#' rw = 2, is.yearly=TRUE, 
#' m = 5, type.st = 4)
#' # Projection
#' out <- projINLA(fit, is.yearly = TRUE)
#' # National smoothed plot
#' plot(out, is.yearly=TRUE, is.subnational=FALSE) + ggplot2::ggtitle("National yearly model")
#' 
#' # National smoothed plot with the aggregated direct estimates
#' plot(out, is.yearly=TRUE, is.subnational=FALSE,  data.add = data_agg, 
#' option.add = list(point = "u5m", lower = "lower", upper = "upper"), 
#' color.add = "orange") + ggplot2::ggtitle("National yearly model") 
#' 
#' # National smoothed plot with the survey-specific direct estimates
#' plot(out, is.yearly=TRUE, is.subnational=FALSE,  data.add = data, 
#' option.add = list(point = "u5m", by = "surveyYears"), 
#' color.add = "darkblue") + ggplot2::ggtitle ("National yearly model") 
#'
#' 
#' 
#' fit <- fitINLA(data = data_agg, geo = geo, Amat = mat, 
#' year_names = years.all, year_range = c(1985, 2019),  
#' rw = 2, is.yearly=TRUE, 
#' m = 5, type.st = 4)
#' # Projection
#' out <- projINLA(fit, Amat = mat, is.yearly = TRUE)
#' 
#' # Subnational estimates
#' plot(out, is.yearly=TRUE, is.subnational=TRUE) + ggplot2::ggtitle("Subnational yearly model")
#'
#' 
#' # Subnational estimates with the aggregated direct estimates
#' plot(out, is.yearly=TRUE, is.subnational=TRUE,  data.add = data_agg, option.add = 
#' list(point = "u5m", lower = "lower", upper = "upper")) + 
#' ggplot2::ggtitle("Subnational yearly model") + facet_wrap(~region)
#'
#' 
#' # Subnational estimates with survey-specific direct estimates
#' plot(out, is.yearly=TRUE, is.subnational=TRUE,  data.add = data, option.add = 
#' list(point = "u5m", by = "surveyYears")) + 
#' ggplot2::ggtitle("Subnational yearly model") + facet_wrap(~region) 
#' 

#' 
#' }
#' @export
plot.projINLA <- function(x, year_label = c("85-89", "90-94", "95-99", "00-04", "05-09", "10-14", "15-19"), year_med = c(1987, 1992, 1997, 2002, 2007, 2012, 2017), is.yearly=TRUE, is.subnational = TRUE, proj_year = 2015, data.add = NULL, option.add = list(point = NULL, lower = NULL, upper = NULL, by = NULL), color.add = "black", dodge.width = 1, ...){

  if(is.null(proj_year)) {
    proj_year = 0
  }


  if(!is.null(data.add)){
    data.add <- data.add[, c("region", "years", as.character(option.add))]
    x <- merge(data.add, x, by = c("region", "years"), all.y = TRUE)
    x$add_x <- x[, option.add$point]
    if(!is.null(option.add$lower)){
      x$add_lower <- x[, option.add$lower]
      x$add_upper <- x[, option.add$upper]      
    }else{
      x$add_lower <- x$add_upper <- NULL
    }
    if(!is.null(option.add$by)){
        x$Comparisons <- x[, option.add$by]
        if(sum(is.na(x$Comparisons)) == dim(x)[1]) x$Comparisons <- "Direct"
    }else{
         x$Comparisons <- "Direct"
    }
  }

  is.periods <- x$years %in% year_label
  x$Year.num[is.periods] <- year_med[match(x$years[is.periods], year_label)]
  x$project <- "No"
  x$project[x$Year.num >= proj_year] <- "Yes"
  # fix for global variable issue
  Year.num <- NULL; region <- NULL; med <- NULL; q025 <- NULL; q975 <- NULL; project <- NULL; Comparisons <- NULL; add_x <- NULL; add_lower <- NULL; add_upper <- NULL
  
  if(is.subnational){
    g <- ggplot2::ggplot(ggplot2::aes(x = Year.num, y = med, ymin = q025, ymax = q975, color = region), data = x)
    my.dodge <- ggplot2::position_dodge(width = dodge.width)
  }else{
    g <- ggplot2::ggplot(ggplot2::aes(x = Year.num, y = med, ymin = q025, ymax = q975), data = x)
    my.dodge <- ggplot2::position_dodge(width = dodge.width/5)
  }
  if(!is.null(data.add)){
    my.dodgeadd <- ggplot2::position_dodge2(width = 0.15, padding = 0.1)
    g <- g + ggplot2::geom_point(data = subset(x, !is.na(Comparisons)), position = my.dodgeadd, ggplot2::aes(x = Year.num+0.5*dodge.width, y = add_x, shape = Comparisons), color = color.add)
    if(!is.null(option.add$lower)) g <- g + ggplot2::geom_errorbar(data = subset(x, !is.na(Comparisons)), position = my.dodgeadd, ggplot2::aes(x = Year.num+0.5*dodge.width, ymin = add_lower, ymax = add_upper), size = 0.5, width = .03, alpha = 0.35, color = color.add)

  }

  if(!is.yearly){
    g <- g + ggplot2::geom_point(position = my.dodge)
    g <- g + ggplot2::geom_line(position = my.dodge)
    g <- g + ggplot2::geom_errorbar(ggplot2::aes(linetype=project), size = .7, width = .05, position = my.dodge)
    g <- g + ggplot2::theme_bw() + ggplot2::xlab("Year") + ggplot2::ylab("U5MR")
    g <- g + ggplot2::scale_x_continuous(breaks=year_med, labels=year_label)
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
    g <- g + ggplot2::geom_point(shape = 17, size = 2.5, position = my.dodge, data=subset(x, is.periods==TRUE), alpha = 0.7)
    g <- g + ggplot2::geom_errorbar(ggplot2::aes(linetype=project), size = .7, width = .05, position = my.dodge, data=subset(x, is.periods==TRUE))
    g <- g + ggplot2::theme_bw() + ggplot2::xlab("Year") + ggplot2::ylab("U5MR")
  }
  g
}