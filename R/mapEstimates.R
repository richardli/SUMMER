#' Mapping estimates for svysae object
#'
#' 
#' @param x syvsae object
#' @param geo.data sf object containing polygon data for the small areas. One of the columns should be named domain and contain the domain labels.
#' @param variable The posterior summary variable to plot. May be one of "median", "mean", or "var".
#' @param viridis.option viridis color scheme
#'
#' @return ggplot containing map of small area posterior summary statistics
#' 
#' @export
#'
#' @examples 
#' \dontrun{
#' data(DemoData2)
#' data(DemoMap2)
#' library(survey)
#' des0 <- svydesign(ids = ~clustid+id, strata = ~strata,
#'                   weights = ~weights, data = DemoData2, nest = TRUE)
#' Xmat <- aggregate(age~region, data = DemoData2, FUN = mean)
#' geo.data <- sf::st_as_sf(DemoMap2$geo)
#' geo.data$domain <- geo.data$REGNAME
#' cts.res <- smoothArea(tobacco.use ~ 1,
#'                       domain = ~region,
#'                       design = des0,
#'                       adj.mat = DemoMap2$Amat, 
#'                       pc.u = 1,
#'                       pc.alpha = 0.01,
#'                       pc.u.phi = 0.5,
#'                       pc.alpha.phi = 2/3,
#'                       return.samples = TRUE)
#' mapEstimates(cts.res, geo.data = geo.data, variable = "median")
#' mapEstimates(cts.res, geo.data = geo.data, variable = "var")
#' }
mapEstimates <- function(x, geo.data, variable, viridis.option = "viridis") {
  combined_est <- do.call(rbind, x[attr(x, "method.names")])
  
  # join estimates and geo, by domain column
  plot_dat <- merge(geo.data, combined_est, all.x = TRUE)
  ggplot2::ggplot(plot_dat, ggplot2::aes(fill = .data[[variable]])) + 
    ggplot2::geom_sf() + 
    ggplot2::facet_wrap(~method) + 
    ggplot2::theme_bw() +
    ggplot2::scale_fill_viridis_c(name = variable, option = viridis.option) 
}
