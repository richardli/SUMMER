#' Function to calculate and plot densities of the projected estimates from INLA for each time and region.
#' 
#' 
#'
#' @param inla_mod output from \code{\link{smoothDirect}}
#' @param results output from \code{\link{getSmoothedDensity}}. This argument can be specified to avoid calculating densities again when only the visualization changes.
#' @param year_range range corresponding to year label
#' @param year_label vector of year string vector
#' @param Amat adjacency matrix
#' @param nsim number of simulations
#' @param weight.strata a data frame with three columns, years, region, and proportion of each strata for the corresponding time period and region. 
#' @param verbose logical indicator whether to print progress messages from inla.posterior.sample.
#' @param mc number of monte carlo draws to approximate the marginal prevalence/hazards for binomial model. If mc = 0, analytical approximation is used. The analytical approximation is invalid for hazard modeling with more than one age groups.
#' @param include_time_unstruct logical indicator whether to include the temporal unstructured effects (i.e., shocks) in the smoothed estimates.
#' @param year_plot similar to year_label, a vector indicate which years to plot
#' @param byyear logical indicator for whether the output uses years as facets. 
#' @param ncol number of columns in the output figure.
#' @param per1000 logical indicator to multiply results by 1000.
#' @param order order of regions when byyear is set to TRUE. Negative values indicate regions are ordered from high to low posterior medians from top to bottom. Positive values indicate from low to high. 0 indicate alphabetic orders.
#' @param ... additional configurations passed to inla.posterior.sample.
#' 
#' @return a data frame of the calculated densities and a ggplot figure.
#' @seealso \code{\link{plot.SUMMERproj}}
#' @examples
#' \dontrun{
#' years <- levels(DemoData[[1]]$time)
#' 
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
#' ## Plot marginal posterior densities over time
#' density <- getSmoothedDensity(fit1,  year_label = years.all, year_plot = years.all, 
#' ncol = 4, byyear = FALSE)
#' density$g
#' 
#' #  subnational model
#' fit2 <- smoothDirect(data = data, geo = geo, Amat = mat, 
#'   year_label = years.all, year_range = c(1985, 2019), 
#'   rw = 2, is.yearly=TRUE, m = 5, type.st = 1)
#' 
##' # Plot marginal posterior densities over time (regions are ordered alphabetically)
#' density <- getSmoothedDensity(fit2, Amat = mat, year_label = years.all, 
#' year_plot = years.all, ncol = 4)
#' density$g
#' 
##' # Re-order the regions 
#' density <- getSmoothedDensity(fit2, Amat = mat,  year_plot = years.all,
#'  ncol = 4, per1000 = TRUE, order = -1)
#' density$g
#' 
#' # Show each region (instead of each year) in a panel 
#' ## Instead of recalculate the posteriors, we can use previously calculated densities as input 
#' density <- getSmoothedDensity(results = density, year_plot = years.all, 
#' ncol = 4, byyear=FALSE, per1000 = TRUE)
#' density$g
#' 
#' # Show more years
#' density <- getSmoothedDensity(results = density, year_plot = c(1990:2019), 
#' ncol = 4, byyear=FALSE, per1000 = TRUE)
#' density$g
#'
#' 
#' 
#' }
#' 

#' @export
getSmoothedDensity <- function(inla_mod=NULL, results = NULL, year_range = c(1985, 2019), year_label = c("85-89", "90-94", "95-99", "00-04", "05-09", "10-14", "15-19"), Amat = NULL, nsim = 1000, weight.strata = NULL, verbose = FALSE, mc = 0, include_time_unstruct = FALSE, year_plot = NULL, byyear = TRUE, ncol = 4, per1000 = FALSE, order = 0, ...){

      years <- x <- y <- `..x..` <- region <- NA

      ########################
      ## Binomial methods
      ########################
    if(!is.null(inla_mod$family)){
            ## Experimental function, finish this later
      stop("This function is experimental and has not been implemented for cluster-level model yet.")

      ########################
      ## Mercer et al. methods
      ########################
    }else{

      if (!isTRUE(requireNamespace("INLA", quietly = TRUE))) {
        stop("You need to install the packages 'INLA'. Please run in your R terminal:\n install.packages('INLA', repos='https://www.math.ntnu.no/inla/R/stable')")
      }
      # If INLA is installed, then attach the Namespace (so that all the relevant functions are available)
      if (isTRUE(requireNamespace("INLA", quietly = TRUE))) {
        if (!is.element("INLA", (.packages()))) {
          attachNamespace("INLA")
        }

      if(is.null(results)){ 
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
        
        names <- expand.grid(area = region_nums, time = timelabel.yearly)
        mod <- inla_mod$fit
        lincombs.info <- inla_mod$lincombs.info
   
          results <- NULL
          for(i in 1:length(timelabel.yearly)){
            for(j in 1:length(region_names)){
                index <- lincombs.info$Index[lincombs.info$District == region_nums[j] & lincombs.info$Year == i]
                tmp <- data.frame(INLA::inla.tmarginal(expit, mod$marginals.lincomb.derived[[index]]))
                tmp$region <- region_nums[j]
                tmp$years <- timelabel.yearly[i]
                results <- rbind(results, tmp)
            }
          }
          results$is.yearly <- !(results$years %in% year_label)
          results$years.num <- suppressWarnings(as.numeric(as.character(results$years)))
          if(region_names[1] != "All"){
            results$region <- region_names[results$region]
          }else{
            results$region <- "All"
          }
          results$years <- factor(results$years, levels = timelabel.yearly)

          # reorder areas
          if(order != 0 && byyear && length(region_names) > 1){
            tmp <- data.frame(region = region_names, median = NA)
            for(j in 1:length(region_names)){
                index <- lincombs.info$Index[lincombs.info$District == region_nums[j] & lincombs.info$Year == length(timelabel.yearly)]
                tmp$median[j]<- INLA::inla.qmarginal(0.5, mod$marginals.lincomb.derived[[index]])
            }
            tmp <- tmp[order(tmp$median, decreasing = (order > 0)), ]
            results$region <- factor(results$region, levels = tmp$region)
          }else{
            results$region <- factor(results$region, rev(sort(as.character(unique(results$region)))))
          }      
      }else{
        results <- results$data
        region_names <- unique(results$region)
        timelabel.yearly <- levels(results$years)
        if(order != 0) warning("Plotting pre-calculated densities, order argument is ignored.")
      }

      results.plot <- results
      if(per1000) results.plot$x <- 1000 * results.plot$x
      if(is.null(year_plot)){
        year_plot <- year_label 
      }
      if(byyear){
          g <- ggplot2::ggplot(subset(results.plot, years %in% year_plot), ggplot2::aes(x = x, y = region, height = y, fill = ..x..)) 
      }else{
         results.plot$years <- factor(results.plot$years, levels = rev(timelabel.yearly))
         g <- ggplot2::ggplot(subset(results.plot, years %in% year_plot), ggplot2::aes(x = x, y = years, height = y, fill = ..x..)) 
      }
      g <- g + ggridges::geom_density_ridges_gradient(stat="identity", alpha = 0.5) 
      g <- g + ggplot2::scale_fill_viridis_c(option = "B") + ggplot2::theme_bw()  + ggplot2::ylab("") + ggplot2::theme(legend.position = 'none') + ggplot2::xlab("")
      if(byyear){
        if(length(year_plot) > 1) g <- g + ggplot2::facet_wrap(~years, ncol = ncol)
      }else{
        if(length(region_names) > 1) g <- g + ggplot2::facet_wrap(~region, ncol = ncol)
      }


      return(list(data = results, g = g))
      }
  }

}


