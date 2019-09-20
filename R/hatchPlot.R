#' Makes map plot with uncertainty hatching.
#' 
#' This function visualizes the map with different variables. The input data frame can be either the long or wide format.
#'
#' @param data a data frame with variables to be plotted
#' @param variables vector of variables to be plotted. If long format of data is used, only one variable can be selected
#' @param values the column corresponding to the values to be plotted, only used when long format of data is used
#' @param labels vector of labels to use for each variable, only used when wide format of data is used
#' @param geo \code{geo} output from \code{\link{read_shape}}
#' @param by.data column name specifying region names in the data
#' @param by.geo variable name specifying region names in the data
#' @param is.long logical indicator of whether the data is in the long format, default to FALSE
#' @param lower column name of the lower bound of the CI
#' @param upper column name of the upper bound of the CI
#' @param lim fixed range of values for the variables to plot
#' @param lim.CI fixed range of the CI widths to plot
#' @param breaks.CI a vector of numerical values that decides the breaks in the CI widths to be shown
#' @param ncol number of columns for the output tabs
#' @param col.hatch color of the hatching lines.
#'
#' @examples
#' \dontrun{
#'  #TODO
#' }
#' @export 

hatchPlot <- function(data, variables, values = NULL, labels = NULL, geo, by.data, by.geo,  is.long = FALSE, lower, upper, lim = NULL, lim.CI = NULL, breaks.CI = NULL, ncol = 4, col.hatch = NULL){


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
        data <- data[, c(variables, by.data)]
        data <- reshape2::melt(data)
        data$variable <- factor(data$variable, levels = variables)
        levels(data$variable) <- labels
    }else {
        data$value <- data[, values]
        data$variable <- data[, variables]
        data$variable <- as.character(data$variable)
        data$variable <- factor(data$variable, levels = labels)
    }
    npal = 100
    med.palette <- viridis::viridis_pal()(npal)
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

    if(is.null(col.hatch)) col.hatch <- viridis::viridis_pal(option="A")(11)[11]
    col.border<- viridis::viridis_pal(option="E")(11)[8]
    
    m <- length(unique(data$variable)) / ncol
    m <- floor(m)
    if(length(unique(data$variable)) %% ncol != 0) m <- m + 1
    par(mfrow = c(m, ncol), mai = c(.25, 0.1,0.3,0.1), oma = c(0.5, 0.1, 0.1, 0.1))
    for(tt in levels(data$variable)){
        tmp <- data[data$variable == tt, ]
        tmp <- tmp[match(geo[[by.geo]], tmp[, by.data]),]
        plot(geo, col = tmp$value.col, border = col.border,
           main  = tt)
      	geo@data$diff <- tmp[, upper] - tmp[, lower]
      	rrt1 <- geo@data$diff
      	plot(geo,  density=dens[findInterval(rrt1, breaks.CI,
             	all.inside=TRUE)], add=T, col = col.hatch, border = FALSE)
    }





    legend.col <- function(col, lev, hadj=-2){
      opar <- par
      n <- length(col)
      bx <- par("usr")
      box.cx <- c(bx[1] + (bx[2] - bx[1]) / 1000,
      bx[1] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 30)
      box.cy <- c(bx[3], bx[3])
      box.sy <- (bx[4] - bx[3]) / n
       
      xx <- rep(box.cx, each = 2)
      par(xpd = TRUE)
      for(i in 1:n){ 
      yy <- c(box.cy[1] + (box.sy * (i - 1)),
      box.cy[1] + (box.sy * (i)),
      box.cy[1] + (box.sy * (i)),
      box.cy[1] + (box.sy * (i - 1)))
      polygon(xx, yy, col = col[i], border = col[i])
      }
      par(new = TRUE)
      plot(0, 0, type = "n",
      ylim = c(min(lev), max(lev)),
      yaxt = "n", ylab = "",
      xaxt = "n", xlab = "",
      frame.plot = FALSE)
      axis(side = 2, las = 2, tick = FALSE, line = .25, hadj=hadj)
      par <- opar
    }

    plot(1, type = "n", axes=FALSE, xlab="", ylab="")
    legend.col(col = med.palette, lev = zlim, hadj = -1.5)

    legend(x = 'center', inset = 0,
           legend = brklabels,
           col = rep('black',2), 
           cex = 1.25, 
           density = dens, bty = 'n', 
           title = "Uncertainty")

}
