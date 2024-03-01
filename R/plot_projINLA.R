#' Plot projection output.
#' @param x output from \code{\link{getSmoothed}}
#' @param year_label labels for the periods
#' @param year_med labels for the middle years in each period, only used when both yearly and period estimates are plotted. In that case, \code{year_med} specifies where each period estimates are aligned.
#' @param is.subnational logical indicator of whether the data contains subnational estimates
#' @param proj_year the first year where projections are made, i.e., where no data are available. 
#' @param data.add data frame for the Comparisons data points to add to the graph. This can be, for example, the raw direct estimates. This data frame is merged to the projections by column 'region' and 'years'. Except for these two columns, this dataset should not have Comparisons columns with names overlapping the getSmoothed output.
#' @param option.add list of options specifying the variable names for the points to plot, lower and upper bounds, and the grouping variable. This is intended to be used to add Comparisons estimates on the same plot as the smoothed estimates. See examples for details.  
#' @param color.add the color of the Comparisons data points to plot.
#' @param label.add the label of the Comparisons data points in the legend.
#' @param dodge.width the amount to add to data points at the same year to avoid overlap. Default to be 1.
#' @param plot.CI logical indicator of whether to plot the error bars.
#' @param per1000 logical indicator to plot mortality rates as rates per 1,000 live births. Note that the added comparison data should always be in the probability scale.
#' @param color.CI the color of the error bars of the credible interval.
#' @param alpha.CI the alpha (transparency) of the error bars of the credible interval.
#' @param ... optional arguments, see details
#' 
#' @method plot SUMMERproj
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_errorbar position_dodge position_dodge2 xlab ylab scale_x_continuous theme_bw scale_shape_discrete
#' @seealso \code{\link{getSmoothed}}
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
#' #  national model
#' years.all <- c(years, "15-19")
#' fit1 <- smoothDirect(data = data, geo = NULL, Amat = NULL, 
#'   year_label = years.all, year_range = c(1985, 2019), 
#'   rw = 2, is.yearly=FALSE, m = 5)
#' out1 <- getSmoothed(fit1)
#' plot(out1, is.subnational=FALSE)
#' 
#' #  subnational model
#' fit2 <- smoothDirect(data = data, geo = geo, Amat = mat, 
#'   year_label = years.all, year_range = c(1985, 2019), 
#'   rw = 2, is.yearly=TRUE, m = 5, type.st = 4)
#' out2 <- getSmoothed(fit2)
#' plot(out2, is.yearly=TRUE, is.subnational=TRUE)
#' 
#' 
#' }
#' 
#' @export
plot.SUMMERproj  <- function(x, year_label = c("85-89", "90-94", "95-99", "00-04", "05-09", "10-14", "15-19"), year_med = c(1987, 1992, 1997, 2002, 2007, 2012, 2017), is.subnational = TRUE, proj_year = 2015, data.add = NULL, option.add = list(point = NULL, lower = NULL, upper = NULL, by = NULL), color.add = "black", label.add = NULL, dodge.width = 1, plot.CI = NULL, per1000 = FALSE,  color.CI = NULL, alpha.CI = 0.5, ...){

  if(is.null(proj_year)) {
    proj_year = 0
  }
  if(is.null(plot.CI)){
    plot.CI <- !is.subnational
  }
  is.yearly = sum(x$is.yearly) > 0
  if(sum(!is.na(x$years)) == 0){
    is.temporal <- FALSE
  }else{
    is.temporal <- TRUE
  }

  if(!is.null(data.add)){
    if(!is.null(option.add$point)){
      for(ii in length(option.add$point)){
        if(option.add$point[ii] %in% colnames(x)){
         colnames(data.add)[which(colnames(data.add) == option.add$point[ii])] <- paste0(option.add$point[ii], ".addedtemp")
         option.add$point[ii] <- paste0(option.add$point[ii], ".addedtemp")
        }
      }
   }
    if("region" %in% colnames(data.add)){
      data.add <- data.add[, c("region", "years", as.character(option.add))]
      x <- merge(data.add, x, by = c("region", "years"), all.y = TRUE)      
    }else{
      data.add <- data.add[, c("years", as.character(option.add))]
      x <- merge(data.add, x, by = c("years"), all.y = TRUE)
    }
    x$add_x <- x[, option.add$point]
    if(!is.null(option.add$lower)){
      x$add_lower <- x[, option.add$lower]
      x$add_upper <- x[, option.add$upper]      
    }else{
      x$add_lower <- x$add_upper <- NULL
    }
    if(!is.null(option.add$by)){
      x$Comparisons <- x[, option.add$by]
        if(sum(is.na(x$Comparisons)) == dim(x)[1]) x$Comparisons <- ifelse(is.null(label.add), "Direct", label.add)
    }else{
        x$Comparisons <-  ifelse(is.null(label.add), "Direct", label.add)
    }
  }

  # deal with 1 year period
  period.1yr <- diff(sort(unique(x$years.num[x$is.yearly == FALSE])))[1] == 1
  if(sum(is.na(x$years.num)) == length(x$years.num)) period.1yr = FALSE

  is.periods <- x$years %in% year_label
  x$years.num[is.periods] <- year_med[match(x$years[is.periods], year_label)]
  x$project <- "No"
  x$project[x$years.num >= proj_year] <- "Yes"
  # fix for global variable issue
  years.num <- NULL; region <- NULL; median <- NULL; lower <- NULL; upper <- NULL; project <- NULL; Comparisons <- NULL; add_x <- NULL; add_lower <- NULL; add_upper <- NULL
  if(per1000){
    x$median <- x$median * 1000
    x$lower <- x$lower * 1000
    x$upper <- x$upper * 1000
    if(!is.null(data.add)){
      x$add_x <- x$add_x * 1000
      if(!is.null(option.add$lower)){
        x$add_lower <- x$add_lower * 1000
        x$add_upper <- x$add_upper * 1000   
        }   
    }
  }
  
 
  # spatial only model
  if(!is.temporal){
    g <- ggplot(data = x)
    g <- g + geom_point(data = subset(x, !is.na(Comparisons)), aes(x = region, y = add_x, shape = Comparisons), color = color.add)
    if(!is.null(option.add$lower)) g <- g + geom_errorbar(data = subset(x, !is.na(Comparisons)),  aes(x = region, ymin = add_lower, ymax = add_upper), size = 0.5, width = .03, alpha = 0.35, color = color.add)
    g <- g + geom_point(aes(x = region, y = median)) + xlab("Region") + ylab("")
    if(plot.CI){
            g <- g + geom_errorbar(aes(ymin = lower, ymax = upper), size = .5, width = .05, alpha = alpha.CI)
    } 

  }else{

  if(is.subnational){
    g <- ggplot(aes(x = years.num, color = region), data = x)
    my.dodge <- position_dodge(width = dodge.width * ifelse(plot.CI, 1, 0))
  }else{
    dodge.width <- dodge.width / 5
    g <- ggplot(aes(x = years.num), data = x)
    my.dodge <- position_dodge(width = dodge.width* ifelse(plot.CI, 1, 0))
  }
  if(!is.null(data.add)){
    my.dodgeadd <- position_dodge2(width = 0.15*dodge.width, padding = 0.1)
    g <- g + geom_point(data = subset(x, !is.na(Comparisons)), position = my.dodgeadd, aes(x = years.num+0.5*dodge.width, y = add_x, shape = Comparisons), color = color.add)
    if(!is.null(option.add$lower)) g <- g + geom_errorbar(data = subset(x, !is.na(Comparisons)), position = my.dodgeadd, aes(x = years.num+0.5*dodge.width, ymin = add_lower, ymax = add_upper), size = 0.5, width = .03, alpha = 0.35, color = color.add)
  }
  # period model
  if(!is.yearly){
    g <- g + geom_point(aes(y = median), position = my.dodge)
    g <- g + geom_line(aes(y = median), position = my.dodge)
    if(plot.CI & !is.null(color.CI)){
        g <- g + geom_errorbar(aes(ymin = lower, ymax = upper, linetype=project), size = .7, width = .05, position = my.dodge, color = color.CI, alpha = alpha.CI)
    }else if(plot.CI &length(unique(x$region))==1){
       g <- g + geom_errorbar(aes(ymin = lower, ymax = upper, linetype=project), size = .7, width = .05, position = my.dodge, color = "black", alpha = alpha.CI)
    }else if(plot.CI & is.null(color.CI)){
        g <- g + geom_errorbar(aes(ymin = lower, ymax = upper, linetype=project, color=region), size = .7, width = .05, position = my.dodge, alpha = alpha.CI)
    }
    g <- g + theme_bw() + xlab("Year") + ylab("")
    if(!period.1yr){
      g <- g + scale_x_continuous(breaks=year_med, labels=year_label)
    }else{
      g <- g + scale_x_continuous(breaks=scales::pretty_breaks())
    }
  
  # yearly model with only one color
  }else if(!is.subnational){
    g <- g + geom_point(aes(y = median), position = my.dodge, data=subset(x, is.periods==FALSE), alpha = 0.5, color = 1)
    g <- g + geom_line(aes(y = median), position = my.dodge, data=subset(x, is.periods==FALSE), alpha = 0.5, color = 1)
    
    g <- g + geom_point(aes(y = median), shape = 17, size = 2.5, position = my.dodge, data=subset(x, is.periods==TRUE), color = 2)
    if(plot.CI){
        g <- g + geom_errorbar(aes(ymin = lower, ymax = upper, linetype=project), size = .5, width = .05, position = my.dodge, data=subset(x, is.periods==FALSE), alpha = alpha.CI, color = "black")
        g <- g + geom_errorbar(aes(ymin = lower, ymax = upper, linetype=project), size = .7, width = .05, position = my.dodge, data=subset(x, is.periods==TRUE), color = "red")
    } 
    g <- g + theme_bw() + xlab("Year") + ylab("") + scale_x_continuous(breaks=scales::pretty_breaks())

  # yearly model with multiple colors
  }else if(is.subnational){
    g <- g + geom_point(aes(y = median), position = my.dodge, data=subset(x, is.periods==FALSE), alpha = 0.5)
    g <- g + geom_line(aes(y = median), position = my.dodge, data=subset(x, is.periods==FALSE), alpha = 0.5)
    # if(plot.CI) g <- g + geom_errorbar(aes(linetype=project), size = .5, width = .05, alpha = alpha.CI, position = my.dodge, color = color.CI)
    g <- g + geom_point(aes(y = median), shape = 17, size = 2.5, position = my.dodge, data=subset(x, is.periods==TRUE), alpha = 0.7)
    if(plot.CI && !is.null(color.CI)){
        g <- g + geom_errorbar(aes(ymin = lower, ymax = upper, linetype=project), size = .7, width = .05, position = my.dodge, data=subset(x, is.periods==TRUE), color = color.CI, alpha = alpha.CI)
    }else if(plot.CI){
        g <- g + geom_errorbar(aes(ymin = lower, ymax = upper, linetype=project, color=region), size = .7, width = .05, position = my.dodge, data=subset(x, is.periods==TRUE), alpha = alpha.CI)
    }
    g <- g + theme_bw() + xlab("Year") + ylab("") + scale_x_continuous(breaks=scales::pretty_breaks())
  }
  }
  if(per1000){
    g <- g + ylab("deaths per 1000 live births")
  }

  if(!is.null(label.add)) g <- g + scale_shape_discrete("") 
  g
}