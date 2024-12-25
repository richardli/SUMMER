#' Plot heatmap comparing pairwise posterior exceedence probabilities for svysae object
#' 
#'
#' @param x  an object in the S3 class of svysae, fhModel, or clusterModel. Plots are created for all models in this object.
#' @param posterior.sample Matrix of posteriors samples of area level quantities with one row for each area and one column for each sample. This argument may be specified to only provide a heatmap for the desired samples.
#' @param title Optional parameter changing the title of the plot
#' @param return.plot Logical indicator for whether the ggplot object is returned
#'
#' @return ggplot containing heat map of pairwise comparisons
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
#' 
#' cts.res <- smoothArea(tobacco.use ~ 1,
#'                       domain = ~region,
#'                       design = des0,
#'                       adj.mat = DemoMap2$Amat, 
#'                       pc.u = 1,
#'                       pc.alpha = 0.01,
#'                       pc.u.phi = 0.5,
#'                       pc.alpha.phi = 2/3,
#'                       return.samples = TRUE)
#' compareEstimates(cts.res)
#' }
compareEstimates <- function(x,
                             posterior.sample = NULL,
                             title = NULL, 
                             return.plot = FALSE) {
  
  
  x_att <- attributes(x)

  if (is.null(posterior.sample)) {
      # USING SURVEYPREV CLASSES
      if(x_att$class %in% c("fhModel", "clusterModel")){
          if ("admin2_post" %in% x_att$names){
              sample_list=list(t(x$admin2_post))
          }else{
              sample_list=list(t(x$admin1_post))}
      # USING SUMMER OBJECTS
      }else{
          sample_list <- x[paste0(x_att$inla.fitted, ".sample")]}
   # USING SAMPLES       
   }else{
      sample_list <- list(posterior.sample)
   }

  for (i in seq_along(sample_list)) {
    current_samples <- sample_list[[i]]
    domain_medians <- apply(current_samples, 1, median)
    median_order <- order(domain_medians)
    
    current_samples <- current_samples[median_order, ]
    # get indices for all possible pairs of admin 1 areas
    domain_comb <- combn(length(x_att$domain.names), 2)
    
    plot_dat <- data.frame(
      Domain1 = x_att$domain.names[median_order][domain_comb[1,]],
      Domain2 = x_att$domain.names[median_order][domain_comb[2,]]
    )
    plot_dat$Domain1 <- factor(plot_dat$Domain1,
                               x_att$domain.names[median_order])
    plot_dat$Domain2 <- factor(plot_dat$Domain2,
                               x_att$domain.names[median_order])
    plot_dat$Prob = apply(
      domain_comb, 2, 
      function(x) mean(current_samples[x[1],] > current_samples[x[2],])
    )
    plot_dat$est = NA
    
    median_dat <- data.frame(
      Domain1 = x_att$domain.names[median_order], 
      Domain2 = "Median", 
      Prob = NA,
      est = domain_medians[median_order]
      
    )
    extra_cols <- c("", "Median", "Interval")
    # combine into single tibble for plotting
    plot_dat <- rbind(plot_dat, median_dat)
    
    plot_dat$Domain2 <- 
      factor(plot_dat$Domain2, levels = c(levels(plot_dat$Domain1), extra_cols))
    
    g_heat <- ggplot2::ggplot(data = plot_dat,
                              ggplot2::aes(x = Domain2, y = Domain1, fill = Prob)) + 
      ggplot2::theme_minimal() + 
      ggplot2::theme(legend.position="bottom", 
                     axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1, size = 10), 
                     panel.grid.major = ggplot2::element_blank(), 
                     panel.grid.minor = ggplot2::element_blank(),
                     text = ggplot2::element_text(size = 10)) +
      ggplot2::geom_tile() + 
      ggplot2::geom_text(
        ggplot2::aes(
          label = ifelse(is.na(Prob), "", sprintf("%0.2f", round(Prob, digits = 2))),
            color = Prob > .5
          ), 
        size = 5
      ) + 
      ggplot2::geom_text(
        data = plot_dat,
        ggplot2::aes(
          label = ifelse(is.na(est), "", sprintf("%0.2f", round(est, digits = 2)))
          ),
        size = 5
      ) +
      ggplot2::coord_equal() + 
      ggplot2::scale_fill_viridis_c(name = "Probability of Domain 1 > Domain 2",
                                    na.value = "white") + 
      ggplot2::scale_color_manual(values = c("white", "black"),  guide = "none") +
      ggplot2::labs(x = "Domain 2",
                    y = "Domain 1")
    if (!is.null(title)) {
      g_heat <- g_heat + ggplot2::labs(title = title)
    } else if (is.null(posterior.sample)) {
      g_heat <-  g_heat + ggplot2::labs(title = x_att$inla.fitted[i])
    }
    if(return.plot) return(g_heat)
    suppressWarnings(print(g_heat))
  }
}
