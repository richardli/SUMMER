#' Plot region-level variables on a map
#' 
#' This function visualizes the map with different variables. The input data frame can be either the long or wide format.
#'
#' @param data a data frame with variables to be plotted. When it is null, a map is produced.
#' @param variables vector of variables to be plotted. If long format of data is used, only one variable can be selected
#' @param values the column corresponding to the values to be plotted, only used when long format of data is used
#' @param labels vector of labels to use for each variable, only used when wide format of data is used
#' @param geo SpatialPolygonsDataFrame object for the map
#' @param by.data column name specifying region names in the data
#' @param by.geo variable name specifying region names in the data
#' @param is.long logical indicator of whether the data is in the long format, default to FALSE
#' @param size size of the border
#' @param removetab logical indicator to not show the tab label, only applicable when only one tab is present.
#' @param border color of the border
#' @param ncol number of columns for the output tabs
#' @param ylim range of the values to be plotted.
#' @param legend.label Label for the color legend.
#' @param per1000 logical indicator to plot mortality rates as rates per 1,000 live births. Note that the added comparison data should always be in the probability scale.
#' @param clean remove all coordinates for a cleaner layout, default to TRUE.
#' @param size.label size of the label of the regions.
#' @param add.adj logical indicator to add edges between connected regions.
#' @param color.adj color of the adjacency matrix edges.
#' @param alpha.adj alpha level (transparency) of the adjacency matrix edges.
#' @param direction Direction of the color scheme. It can be either 1 (smaller values are darker) or -1 (higher values are darker). Default is set to 1.
#' @param cut a vector of values to cut the continuous scale color to discrete intervals. 
#' @importFrom sp proj4string
#' @importFrom shadowtext geom_shadowtext
#' @importFrom sp Polygon
#' @importFrom stats setNames
#' @author Zehang Richard Li
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
mapPlot <- function(data = NULL, variables, values = NULL, labels = NULL, geo, by.data, by.geo, is.long = FALSE, size = 0.5, removetab = FALSE, border = "gray20", ncol = NULL, ylim = NULL, legend.label = NULL,  per1000 = FALSE, clean = TRUE, size.label = 2, add.adj = FALSE, color.adj = "red", alpha.adj = 0.85, direction = 1, cut = NULL){
    value <- group <- lat <- long <- x0 <- x1 <- y0 <- y1 <- id <- name <- variable <- NULL

    # Simple Map Plot
    if(is.null(data)){
        geo1 <- ggplot2::fortify(geo, region = by.geo)
        geo2 <- by(geo1, geo1$id, function(x) {sp::Polygon(x[c('long', 'lat')])@labpt})
        centroids <- stats::setNames(do.call("rbind.data.frame", geo2), c('long', 'lat'))
        centroids$name <- names(geo2) 

        g <- ggplot2::ggplot() + ggplot2::geom_polygon(data=geo1, ggplot2::aes(x = long, y = lat, group = group, fill = id), color = border, alpha = .3) + ggplot2::coord_map()  

        if(add.adj){
            mat <- getAmat(geo, names = geo[[by.geo]])
            centroids2 <- centroids[match(colnames(mat), centroids$name), ]
            edges <- data.frame(x0 = rep(NA, sum(mat)/2), x1 = rep(NA, sum(mat)/2), y0 = rep(NA, sum(mat)/2), y1 = rep(NA, sum(mat)/2))
            index <- 1
            for(i in 1:dim(mat)[1]){
                for(j in 1:i){
                    if(mat[i, j] == 1){
                        edges[index, ] <- c(centroids2[i, 1], centroids2[j, 1], centroids2[i, 2], centroids2[j, 2])
                        index <- index + 1
                }
            }
        }
        g <- g + ggplot2::geom_segment(data = edges, ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1), color = color.adj, alpha = alpha.adj)
        }

        g <- g + ggplot2::scale_fill_discrete(guide = FALSE) + shadowtext::geom_shadowtext(data = centroids, ggplot2::aes(x = long, y = lat, label = name), check.overlap = TRUE, size = size.label, color = border, bg.colour='white') + ggplot2::theme_bw()

        if(clean){
            g <- g + ggplot2::theme_bw() + ggplot2::theme(legend.title=ggplot2::element_text(size=ggplot2::rel(0.7)), axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(), axis.ticks.x=ggplot2::element_blank(), axis.title.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.ticks.y=ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())
        }
        return(g)
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
    has.coord <- !is.na(sp::proj4string(geo))
    regions <- as.character(unique(geo[[by.geo]]))
    geo <- ggplot2::fortify(geo, region = by.geo)
    if (!is.long) {
        data <- data[, c(variables, by.data)]
        nonexist <- regions[regions %in% as.character(data[, by.data]) == FALSE]
        if(length(nonexist) > 0){
            data.add <- data[rep(1, length(nonexist)), ]
            data.add[, 1:ncol(data.add)] <- NA
            data.add[, by.data] <- nonexist
            data <- rbind(data, data.add)
            warning(paste0("The following areas not in the dataset to plot: ", paste(nonexist, collapse = ", ")))
        }
        data <- reshape2::melt(data)
        data$variable <- factor(data$variable, levels = variables)
        levels(data$variable) <- labels
    }
    else {
        data$value <- data[, values]
        data$variable <- data[, variables]
        nonexist.any <- NULL
        alllevels <- unique(as.character(data$variable))
        for(jj in 1:length(alllevels)){
            sub <- subset(data, variable == alllevels[jj])
            nonexist <- regions[regions %in% as.character(sub[, by.data]) == FALSE]
            if(length(nonexist) > 0){
                    data.add <- data[rep(1, length(nonexist)), ]
                    data.add[, 1:ncol(data.add)] <- NA
                    data.add[, by.data] <- nonexist
                    data.add[, "variable"] <- alllevels[jj]
                    data <- rbind(data, data.add)
            }
            nonexist.any <- c(nonexist.any, nonexist)
        }
        if(length(nonexist.any) > 0) warning(paste0("The following areas contain missing values: ", paste(unique(nonexist.any), collapse = ", ")))            
        data$variable <- as.character(data$variable)
        data$variable <- factor(data$variable, levels = labels)
        data <- data[, c(by.data, "variable", "value")]
    }
    if(per1000){
        data$value <- data$value * 1000
    }
    geo2 <- merge(geo, data, by = "id", by.y = by.data)
    g <- ggplot2::ggplot(geo2)
    if(!is.null(cut)){
        g <- g + ggplot2::geom_polygon(ggplot2::aes(x = long, y = lat, 
            group = group, fill = cut(value, cut)), color = border, size = size)
    }else{
        g <- g + ggplot2::geom_polygon(ggplot2::aes(x = long, y = lat, 
            group = group, fill = value), color = border, size = size)
    }
    if(length(unique(data$variable)) > 1 || removetab == FALSE){
        if(is.null(ncol)){
            g <- g + ggplot2::facet_wrap(~variable)
        }else{
            g <- g + ggplot2::facet_wrap(~variable, ncol = ncol)
        }
    }
    if(is.null(legend.label)){
        legend.label <- "Value"
    }
    if(!is.null(ylim) && is.null(cut)){
        g <- g + ggplot2::scale_fill_viridis_c(legend.label, lim = ylim, direction = direction)
    }else if(is.null(cut)){
        g <- g + ggplot2::scale_fill_viridis_c(legend.label, direction = direction)
    }else if(!is.null(cut)){
        g <- g + ggplot2::scale_fill_viridis_d(legend.label, direction = direction)
    }
    if(has.coord) g <- g + ggplot2::coord_map()

    if(clean){
        g <- g + ggplot2::theme_bw() + ggplot2::theme(legend.title=ggplot2::element_text(size=ggplot2::rel(0.7)), axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(), axis.ticks.x=ggplot2::element_blank(), axis.title.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.ticks.y=ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())
    }

    return(g)
}