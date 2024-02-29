#' Plot maps with uncertainty hatching.
#' 
#' This function visualizes the map with different variables. The input data frame can be either the long or wide format.
#'
#' @param data a data frame with variables to be plotted
#' @param variables vector of variables to be plotted. If long format of data is used, only one variable can be selected
#' @param values the column corresponding to the values to be plotted, only used when long format of data is used
#' @param labels vector of labels to use for each variable, only used when wide format of data is used
#' @param geo SpatialPolygonsDataFrame object for the map
#' @param by.data column name specifying region names in the data
#' @param by.geo variable name specifying region names in the data
#' @param is.long logical indicator of whether the data is in the long format, default to FALSE
#' @param lower column name of the lower bound of the CI
#' @param upper column name of the upper bound of the CI
#' @param lim fixed range of values for the variables to plot
#' @param lim.CI fixed range of the CI widths to plot
#' @param breaks.CI a vector of numerical values that decides the breaks in the CI widths to be shown
#' @param ncol number of columns for the output tabs
#' @param hatch color of the hatching lines.
#' @param border color of the polygon borders.
#' @param size line width of the polygon borders.
#' @param legend.label Label for the color legend.
#' @param per1000 logical indicator to plot mortality rates as rates per 1,000 live births. Note that the added comparison data should always be in the probability scale.
#' @param direction Direction of the color scheme. It can be either 1 (smaller values are darker) or -1 (higher values are darker). Default is set to 1.
#' @param ... unused.
#'
#' @author Zehang Richard Li, Katie Wilson
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
#' fit2 <- smoothDirect(data = data, geo = geo, Amat = mat, 
#'   year_label = years.all, year_range = c(1985, 2019), 
#'   rw = 2, is.yearly=TRUE, m = 5, type.st = 4)
#' out2 <- getSmoothed(fit2)
#' 
#' plot(out2, is.yearly=TRUE, is.subnational=TRUE)
#' 
#' hatchPlot(data = subset(out2, is.yearly==FALSE), geo = geo,
#' variables=c("years"), values = c("median"), 
#' by.data = "region", by.geo = "REGNAME", 
#' lower = "lower", upper = "upper", is.long=TRUE)
#' 
#' }
#' @importFrom viridis viridis_pal
#' @importFrom sp plot
#' @importFrom graphics text
#' @export 

hatchPlot <- function(data, variables, values = NULL, labels = NULL, geo, by.data, by.geo,  is.long = FALSE, lower, upper, lim = NULL, lim.CI = NULL, breaks.CI = NULL, ncol = 4, hatch = NULL, border = NULL, size = 1, legend.label = NULL,  per1000 = FALSE, direction = 1, ...){
    is.sf <- "sf" %in% class(geo)
    if(is.sf){
       geo =  as(geo, 'Spatial')
    }

   if (is.null(labels) & !is.long) {
        labels <- variables
    }
    if (is.null(labels) & is.long) {
        labels <- sort(unique(data[, variables]))
    }
    if (is.null(values) & is.long) {
        stop("values need to be specified for long format input.")
    }
    if (!is.long) {
        data2 <- data[, c(variables, by.data)]
        data2 <- reshape2::melt(data2)
        data2$variable <- factor(data2$variable, levels = variables)
        levels(data2$variable) <- labels
        data <- merge(data2, data, all.x=TRUE)
    }else {
        data$value <- data[, values]
        data$variable <- data[, variables]
        data$variable <- as.character(data$variable)
        data$variable <- factor(data$variable, levels = labels)
    }
    if(per1000){
      data$value <- data$value * 1000
      data[, upper] <- data[, upper] * 1000
      data[, lower] <- data[, lower] * 1000
    }
    npal = 100
    med.palette <- viridis::viridis_pal()(npal)
    if(direction == -1) med.palette <- med.palette[length(med.palette):1]

    if(is.null(lim)){
      zlim <- range(data$value, na.rm=TRUE)
    }else{
      zlim <- lim
    }
    data$value.col <- med.palette[floor((data$value - zlim[1])/(zlim[2] - zlim[1])*(npal-1))+1]

    if(is.null(lim.CI)){
      ylim <- range(c(data[, upper] - data[, lower]), na.rm = TRUE)
    }else{
      ylim <- lim.CI
    }
    if(is.null(breaks.CI)){
      breaks.CI <- seq(min(data[, upper] - data[, lower], na.rm = TRUE), 
                       max(data[, upper] - data[, lower], na.rm = TRUE), length.out = 11)
    }
    nbrks <- length(breaks.CI)
    brklabels <- paste(signif(breaks.CI[1:(nbrks-1)],2), signif(breaks.CI[2:nbrks],2), sep = " - ")
    dens <- (2:nbrks)*3

    if(is.null(hatch)){
      col.hatch <- viridis::viridis_pal(option="A")(11)[11]
      if(direction == -1) col.hatch <- viridis::viridis_pal(option="A")(11)[1]
    }else{
      col.hatch <- hatch
    }
    if(is.null(border)){
      col.border<- viridis::viridis_pal(option="E")(11)[8]
    }else{
      col.border <- border
    }
    
    nplot <- length(unique(data$variable))
    if(ncol == nplot) ncol <- ncol + 1
    m <- nplot / ncol
    m <- floor(m) + 1
    graphics::par(mfrow = c(m, ncol), mai = c(.25, 0.1,0.3,0.1), oma = c(0.5, 0.1, 0.1, 0.1))
    for(tt in levels(data$variable)){
        tmp <- data[data$variable == tt, ]
        tmp <- tmp[match(geo[[by.geo]], tmp[, by.data]),]
        sp::plot(geo, col = tmp$value.col, border = col.border,
           main  = tt, lwd = size)
      	geo@data$diff <- tmp[, upper] - tmp[, lower]
      	rrt1 <- geo@data$diff
      	sp::plot(geo,  density=dens[findInterval(rrt1, breaks.CI,
             	all.inside=TRUE)], add=T, col = col.hatch, border = FALSE)
    }
    nplot0 <- ncol - nplot %% ncol 

    legend.col <- function(col, lev, title = NULL){
      opar <- par
      n <- length(col)
      bx <- graphics::par("usr") 
      box.cx <- c(bx[1] + (bx[2] - bx[1]) / 1000,
      bx[1] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 30)
      box.cy <- c(bx[3], bx[3])
      box.sy <- (bx[4] - bx[3]) / n
       
      xx <- rep(box.cx, each = 2)
      graphics::par(xpd = TRUE)
      for(i in 1:n){ 
      yy <- c(box.cy[1] + (box.sy * (i - 1)),
      box.cy[1] + (box.sy * (i)),
      box.cy[1] + (box.sy * (i)),
      box.cy[1] + (box.sy * (i - 1)))
      graphics::polygon(xx, yy, col = col[i], border = col[i])
      }
      if(!is.null(title)) graphics::text(box.cx[1], box.cy[1] + box.sy * (n+2), title, pos=4, cex = 1.2)
      graphics::par(new = TRUE)
      graphics::plot(0, 0, type = "n",
      ylim = c(min(lev), max(lev)),
      xlim = c(bx[1], bx[2]),
      yaxt = "n", ylab = "",
      xaxt = "n", xlab = "",
      frame.plot = FALSE)
      z = graphics::axTicks(side = 2)
      graphics::text(box.cx[1], z, z, pos=4, cex = 1.2)
      par <- opar
    }
    for(tt in 1:nplot0){
        graphics::plot(0, type = "n", axes=FALSE, xlab="", ylab="")      
    }
    legend.col(col = med.palette, lev = zlim, title = legend.label)


    graphics::legend(x = 'center', inset = 0,
           legend = brklabels,
           col = rep('black',2), 
           cex = 1.4,  
           density = dens, bty = 'n', 
           title = "Uncertainty")

}
