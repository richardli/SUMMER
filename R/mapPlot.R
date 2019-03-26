#' Makes map plot.
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
#' @param size size of the border
#' @param removetab logical indicator to not show the tab label, only applicable when only one tab is present.
#' @param border color of the border
#' @param ncol number of columns for the output tabs
#' @param ylim range of the values to be plotted.
#' 
#' @examples
#' \dontrun{
#' data(DemoMap)
#' # Plotting data in the long format
#' dat <- data.frame(region = rep(c("central",  "eastern", "northern", "western"), 3),
#' year = rep(c(1980, 1990, 2000), each = 4),
#' values = stats::rnorm(12))
#' utils::head(dat)
#' mapPlot(dat, variables = "year", values = "values",
#' by.data = "region", geo = DemoMap$geo,
#' by.geo = "NAME_final", is.long = TRUE)

# Plotting data in the wide format
#' dat <- data.frame(region = c("central",  "eastern", "northern", "western"),
#' Year1 = stats::rnorm(4), Year2 = stats::rnorm(4),
#' Year3 = stats::rnorm(4))
#' utils::head(dat)
#' mapPlot(dat, variables = c("Year1", "Year2", "Year3"),
#'  labels = c(1980, 1990, 2000),
#' by.data = "region", geo = DemoMap$geo,
#' by.geo = "NAME_final", is.long = FALSE)
#' 
#' }
#' 
#' @export
mapPlot <- function(data, variables, values = NULL, labels = NULL, geo, by.data, by.geo, is.long = FALSE, size = 0.5, removetab = FALSE, border = "gray20", ncol = NULL, ylim = NULL){
    value <- group <- lat <- long <- NULL
    if (is.null(labels) & !is.long) {
        labels <- variables
    }
    if (is.null(labels) & is.long) {
        labels <- sort(unique(data[, variables]))
    }
    if (is.null(values) & is.long) {
        stop("values need to be specified for long format input.")
    }
    geo <- ggplot2::fortify(geo, region = by.geo)
    if (!is.long) {
        data <- data[, c(variables, by.data)]
        data <- reshape2::melt(data)
        data$variable <- factor(data$variable, levels = variables)
        levels(data$variable) <- labels
    }
    else {
        data$value <- data[, values]
        data$variable <- data[, variables]
        data$variable <- as.character(data$variable)
        data$variable <- factor(data$variable, levels = labels)
    }
    geo2 <- merge(geo, data, by = "id", by.y = by.data)
    g <- ggplot2::ggplot(geo2)
    g <- g + ggplot2::geom_polygon(ggplot2::aes(x = long, y = lat, 
        group = group, fill = value), color = border, size = size)
    if(length(unique(data$variable)) > 1 || removetab == FALSE){
        if(is.null(ncol)){
            g <- g + ggplot2::facet_wrap(~variable)
        }else{
            g <- g + ggplot2::facet_wrap(~variable, ncol = ncol)
        }
    }
    if(!is.null(ylim)){
        g <- g + ggplot2::scale_fill_viridis_c(lim = ylim)
    }else{
        g <- g + ggplot2::scale_fill_viridis_c()
    } 

    return(g)
}


## Visualize covariates (Space-invariant)