#' Calculate and plot posterior densities of the projected estimates
#' 
#' The function \code{ridgePlot} replaces the previous function name \code{getSmoothedDensity} (before version 1.0.0).
#'
#' @param x output from \code{\link{smoothDirect}} for the smoothed direct estimates, or \code{\link{smoothCluster}} for the cluster-level estimates.
#' @param nsim number of posterior draws to take. Only used for cluster-level models when \code{draws} is NULL. Otherwise the posterior draws in \code{draws} will be used instead without resampling.
#' @param draws Output of \code{\link{getSmoothed}} with \code{save.draws} set to TRUE. This argument allows the previously sampled draws (by setting \code{save.draws} to be TRUE) be used in new aggregation tasks. This argument is only used for cluster-level models.   
#' @param year_plot A vector indicate which years to plot
#' @param strata_plot Name of the strata to plot. If not specified, the overall is plotted.
#' @param by.year logical indicator for whether the output uses years as facets. 
#' @param ncol number of columns in the output figure.
#' @param scale numerical value controlling the height of the density plots.
#' @param per1000 logical indicator to multiply results by 1000.
#' @param order order of regions when by.year is set to TRUE. Negative values indicate regions are ordered from high to low posterior medians from top to bottom. Positive values indicate from low to high. 0 indicate alphabetic orders.
#' @param direction Direction of the color scheme. It can be either 1 (smaller values are darker) or -1 (higher values are darker). Default is set to 1.
#' @param results output from \code{\link{ridgePlot}} returned object with \code{save.density = TRUE}. This argument can be specified to avoid calculating densities again when only the visualization changes.
#' @param save.density Logical indicator of whether the densities will be returned with the ggplot object. If set to TRUE, the output will be a list consisting of (1) a data frame of computed densities and (2) a ggplot object of the plot. 
#' @param ... additional configurations passed to inla.posterior.sample.
#' 
#' @return ridge plot of the density, and  if \code{save.density = TRUE}, also a data frame of the calculated densities 
#' @seealso \code{\link{plot.SUMMERproj}}
#' @author Zehang Richard Li
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
#'   rw = 2, m = 5)
#' ## Plot marginal posterior densities over time
#' ridgePlot(fit1, year_plot = years.all, 
#'           ncol = 4, by.year = FALSE)
#' 
#' #  subnational model
#' fit2 <- smoothDirect(data = data, geo = DemoMap$geo, Amat = DemoMap$Amat, 
#'   year_label = years.all, year_range = c(1985, 2019), 
#'   rw = 2, m = 5, type.st = 1)
#' 
##' # Plot marginal posterior densities over time (regions are ordered alphabetically)
#' ridgePlot(fit2, year_plot = years.all, ncol = 4)
#' 
##' # Re-order the regions and save the density to avoid re-compute later
#' density <- ridgePlot(fit2, year_plot = years.all,
#'  ncol = 4, per1000 = TRUE, order = -1, save.density = TRUE)
#' density$g
#' 
#' # Show each region (instead of each year) in a panel 
#' ## Instead of recalculate the posteriors, we can use previously calculated densities as input 
#' ridgePlot(results = density, year_plot = years.all, 
#' ncol = 4, by.year=FALSE, per1000 = TRUE)
#' 
#' # Show more years
#' ridgePlot(results = density, year_plot = c(1990:2019), 
#' ncol = 4, by.year=FALSE, per1000 = TRUE)
#'
#' 
#' 
#' }
#' 

#' @export
ridgePlot <- function(x=NULL, nsim = 1000, draws = NULL, year_plot = NULL, strata_plot = NULL, by.year = TRUE, ncol = 4, scale = 2, per1000 = FALSE, order = 0, direction = 1, results = NULL, save.density = FALSE, ...){

      years <-  y <- `..x..` <- region <- NA

    
      if (!isTRUE(requireNamespace("INLA", quietly = TRUE))) {
        stop("You need to install the packages 'INLA'. Please run in your R terminal:\n  install.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/stable'), dep=TRUE)")
      }
      # If INLA is installed, then attach the Namespace (so that all the relevant functions are available)
      if (isTRUE(requireNamespace("INLA", quietly = TRUE))) {
        if (!is.element("INLA", (.packages()))) {
          attachNamespace("INLA")
        }

     is.density <- FALSE
     if(!is.null(results)){
        results <- results$data
        region_names <- unique(results$region)
        timelabel.yearly <- levels(results$years)
        if(order != 0) warning("Plotting pre-calculated densities, order argument is ignored.")
        is.density <- "y" %in% colnames(results) 

     }else{ 
        Amat <- x$Amat
        if(is.null(Amat) && is.null(draws)){
          region_names <- "All"
          region_nums <- 0
        
        }else if(!is.null(draws)){
          region_names <- unique(draws$overall$region)
          region_nums <- 1:length(region_names)
        
        }else{
          region_names <- colnames(Amat)
          region_nums <- 1:length(region_names)
        }      
        ########################
        ## Binomial methods
        ########################
        # when posterior draws exist, make that into x 
        if(!is.null(draws)){
          if("draws.est" %in% draws == FALSE && is(x, "list")) stop("draws argument is not correctly specified. It should be the whole output of getSmoothed() function.")
          message("Use posterior draws from input.")
          x <- draws
        # when draws does not exist
        }else if(is.null(x)){
          # only take draws for cluster-level models
          if(!is.null(x$family)){
            message("Draws not specified. Use getSmoothed() to calculate posterior draws. Please be aware this does not take into account strata weighting.")
            draws <- getSmoothed(x, Amat = Amat, nsim = nsim, save.draws = TRUE)
            # overwrite x...
            x <- draws
          }
        }else if(is.null(x$fit)){
          stop("Neither fitted object nor posterior draws are provided.")
        }

        # at this point x is an object from getSmoothed or smoothDirect
        ##
        ##  If the input is from getSmoothed
        ##
        if((is(x, "list") || is(x, "SUMMERprojlist")) && is.null(x$fit)){
          is.density <- FALSE
          if(is.null(x$draws.est.overall)){
            stop("Posterior draws not found. Please rerun getSmoothed() with save.draws = TRUE.")
          }
          tmp <- NULL
          if(is.null(strata_plot)){
            draws.plot <- x$draws.est.overall
          }else{
             draws.plot <- NULL
             counter <- 1
             for(i in 1:length(x$draws.est)){
              if(x$draws.est[[i]]$strata == draws.plot){
                draws.plot[[counter]] <- x$draws.est[[i]]
                counter <- counter + 1
              }
            }
          }

          # draw.est is ordered
          for(i in 1:length(draws.plot)) {
            if(draws.plot[[i]]$years %in% tmp) next
            tmp <- c(tmp, draws.plot[[i]]$years)
          }
          year_label <- tmp 
        
          timelabel.yearly <- year_label
          results <- NULL

          for(i in 1:length(timelabel.yearly)){
            for(j in 1:length(region_names)){
              draw_est_j <- draws.plot[lapply(draws.plot, '[[',"region")==region_names[j]]
              draw_est_ij <- draw_est_j[lapply(draw_est_j, '[[', "years")==timelabel.yearly[i]]
              tmp <- data.frame(draws = draw_est_ij[[1]]$draws)
              tmp$region <- region_nums[j]
              tmp$years <- timelabel.yearly[i]
              results <- rbind(results, tmp)
            }
          }
          results$years.num <- suppressWarnings(as.numeric(as.character(results$years)))
          results$x <- results$draws
          if(region_names[1] != "All"){
            results$region <- region_names[results$region]
          }else{
            results$region <- "All"
          }
          results$years <- factor(results$years, levels = timelabel.yearly)

          # reorder areas
          if(order != 0 && by.year && length(region_names) > 1){
            tmp <- data.frame(region = region_names, median = NA)
            for(j in 1:length(region_names)){
                tmp$median[j]<- median(results[results$region == region_names[j] & results$years == timelabel.yearly[length(timelabel.yearly)], ]$draws, na.rm = T)
            }
            tmp <- tmp[order(tmp$median, decreasing = (order > 0)), ]
            results$region <- factor(results$region, levels = tmp$region)
          }else{
            results$region <- factor(results$region, rev(sort(as.character(unique(results$region)))))
          }      
        ########################
        ## Mercer et al. methods
        ########################
        }else{
          is.density <- TRUE
          is.yearly = x$is.yearly
          year_label <- x$year_label
          if(is.yearly){
            timelabel.yearly <- c(x$year_range[1] : x$year_range[2], year_label)
          }else{
            timelabel.yearly <- year_label
          }

          names <- expand.grid(area = region_nums, time = timelabel.yearly)
          mod <- x$fit
          lincombs.info <- x$lincombs.info

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
          if(order != 0 && by.year && length(region_names) > 1){
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
        }      
    }

      results.plot <- results
      if(per1000) results.plot$x <- 1000 * results.plot$x
      if(is.null(year_plot)){
        year_plot <- year_label 
      }
      # plot calculated density
      if(by.year && is.density){
         g <- ggplot2::ggplot(subset(results.plot, years %in% year_plot), ggplot2::aes(x = x, y = region, height = y, fill = ..x..)) 
      }else if(is.density){
         results.plot$years <- factor(results.plot$years, levels = rev(timelabel.yearly))
         g <- ggplot2::ggplot(subset(results.plot, years %in% year_plot), ggplot2::aes(x = x, y = years, height = y, fill = ..x..)) 
      }
      if(is.density) g <- g + ggridges::geom_density_ridges_gradient(stat="identity", alpha = 0.5, size = 0.3) 

      # plot draws
      if(by.year && !is.density){
        g <- ggplot2::ggplot(subset(results.plot, years %in% year_plot), ggplot2::aes(x = x, y = region)) 
      }else if(!is.density){
        results.plot$years <- factor(results.plot$years, levels = rev(timelabel.yearly))
        g <- ggplot2::ggplot(subset(results.plot, years %in% year_plot), ggplot2::aes(x = x, y = years)) 
      }
      if(!is.density) g <- g + ggridges::geom_density_ridges_gradient(ggplot2::aes(fill = ..x..), scale = scale, size = 0.3, alpha = 0.5)

      g <- g + ggplot2::scale_fill_viridis_c(option = "D", direction = direction) + ggplot2::theme_bw()  + ggplot2::ylab("") + ggplot2::theme(legend.position = 'none') + ggplot2::xlab("")
      if(by.year){
        if(length(year_plot) > 1) g <- g + ggplot2::facet_wrap(~years, ncol = ncol)
      }else{
        if(length(region_names) > 1) g <- g + ggplot2::facet_wrap(~region, ncol = ncol)
      }

      if(save.density){
        return(list(data = results, g = g))        
      }else{
        return(g)
      }

  }

}



#' @export
#' @rdname ridgePlot
getSmoothedDensity <- ridgePlot