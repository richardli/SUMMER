#' Calculate and plot posterior densities of the projected estimates
#' 
#' The function \code{ridgePlot} replaces the previous function name \code{getSmoothedDensity} (before version 1.0.0).
#'
#' @param x output from \code{\link{smoothDirect}} for the smoothed direct estimates, or \code{\link{smoothCluster}} for the cluster-level estimates.
#' @param nsim number of posterior draws to take. Only used for cluster-level models when \code{draws} is NULL. Otherwise the posterior draws in \code{draws} will be used instead without resampling.
#' @param draws Output of \code{\link{getSmoothed}} with \code{save.draws} set to TRUE. This argument allows the previously sampled draws (by setting \code{save.draws} to be TRUE) be used in new aggregation tasks. This argument is only used for cluster-level models.   
#' @param year.plot A vector indicate which years to plot
#' @param year_plot `r lifecycle::badge("deprecated")` replaced by year.plot
#' @param strata.plot Name of the strata to plot. If not specified, the overall is plotted.
#' @param strata_plot `r lifecycle::badge("deprecated")` replaced by strata.plot
#' @param by.year logical indicator for whether the output uses years as facets. 
#' @param ncol number of columns in the output figure.
#' @param scale numerical value controlling the height of the density plots.
#' @param per1000 logical indicator to multiply results by 1000.
#' @param order order of regions when by.year is set to TRUE. Negative values indicate regions are ordered from high to low posterior medians from top to bottom. Positive values indicate from low to high. 0 indicate alphabetic orders.
#' @param direction Direction of the color scheme. It can be either 1 (smaller values are darker) or -1 (higher values are darker). Default is set to 1.
#' @param linewidth width of the ridgeline.
#' @param results output from \code{\link{ridgePlot}} returned object with \code{save.density = TRUE}. This argument can be specified to avoid calculating densities again when only the visualization changes.
#' @param save.density Logical indicator of whether the densities will be returned with the ggplot object. If set to TRUE, the output will be a list consisting of (1) a data frame of computed densities and (2) a ggplot object of the plot. 
#' @param ... additional configurations passed to inla.posterior.sample.
#' 
#' @return ridge plot of the density, and  if \code{save.density = TRUE}, also a data frame of the calculated densities 
#' @seealso \code{\link{plot.SUMMERproj}}
#' @importFrom methods as
#' @importFrom stats rnorm
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
#'   year.label = years.all, year.range = c(1985, 2019), 
#'   rw = 2, m = 5)
#' ## Plot marginal posterior densities over time
#' ridgePlot(fit1, year.plot = years.all, 
#'           ncol = 4, by.year = FALSE)
#' 
#' #  subnational model
#' fit2 <- smoothDirect(data = data, geo = DemoMap$geo, Amat = DemoMap$Amat, 
#'   year.label = years.all, year.range = c(1985, 2019), 
#'   rw = 2, m = 5, type.st = 1)
#' 
##' # Plot marginal posterior densities over time (regions are ordered alphabetically)
#' ridgePlot(fit2, year.plot = years.all, ncol = 4)
#' 
##' # Re-order the regions and save the density to avoid re-compute later
#' density <- ridgePlot(fit2, year.plot = years.all,
#'  ncol = 4, per1000 = TRUE, order = -1, save.density = TRUE)
#' density$g
#' 
#' # Show each region (instead of each year) in a panel 
#' ## Instead of recalculate the posteriors, we can use previously calculated densities as input 
#' ridgePlot(results = density, year.plot = years.all, 
#' ncol = 4, by.year=FALSE, per1000 = TRUE)
#' 
#' # Show more years
#' ridgePlot(results = density, year.plot = c(1990:2019), 
#' ncol = 4, by.year=FALSE, per1000 = TRUE)
#'
#' 
#' # Example using surveyPrev package output
#' 
#' library(surveyPrev)
#' dhsData <- getDHSdata(country = "Rwanda", indicator = "nmr", year = 2019)
#' data <- getDHSindicator(dhsData, indicator = "nmr")
#' geo <- getDHSgeo(country = "Rwanda", year = 2019)
#' poly.adm1 <- geodata::gadm(country="RWA", level=1, path=tempdir())
#' poly.adm1 <- sf::st_as_sf(poly.adm1)
#' poly.adm2 <- geodata::gadm(country="RWA", level=2, path=tempdir())
#' poly.adm2 <- sf::st_as_sf(poly.adm2)
#' cluster.info <- clusterInfo(geo = geo, 
#'               poly.adm1 = poly.adm1, 
#'               poly.adm2 = poly.adm2,
#'                             by.adm1 = "NAME_1", 
#'                             by.adm2 = "NAME_2")
#' 
#' fit1 <- directEST(data = data, cluster.info = cluster.info,  admin = 1)
#' fit2 <- directEST(data = data, cluster.info = cluster.info,  admin = 2) 
#' ridgePlot(fit1, direction = -1)
#' ridgePlot(fit2, direction = -1)
#' 
#' }
#' 

#' @export
ridgePlot <- function(x=NULL, nsim = 1000, draws = NULL, year.plot = NULL, year_plot = deprecated(), strata.plot = NULL, strata_plot = deprecated(), by.year = TRUE, ncol = 4, scale = 2, per1000 = FALSE, order = 0, direction = 1, linewidth = 0.5, results = NULL, save.density = FALSE, ...){


      if (lifecycle::is_present(year_plot)) {
          lifecycle::deprecate_soft("2.0.0", "ridgePlot(year_plot)", "ridgePlot(year.plot)")
          year.plot <- year_plot
      }
      if (lifecycle::is_present(strata_plot)) {
          lifecycle::deprecate_soft("2.0.0", "ridgePlot(strata_plot)", "ridgePlot(strata.plot)")
          strata.plot <- strata_plot
      }
      years <-  y <- `..x..` <- region <- value <- region.name <- admin2.name.short <- NA

      # FOR SURVEYPREV INPUT

      if(class(x) %in% c("fhModel", "clusterModel", "directEST")){
           x_att <- attributes(x)
           domain.names <- x_att$domain.names
        # USING SURVEYPREV CLASSES: smoothed model
        if(x_att$class %in% c("fhModel", "clusterModel")){
          if ("admin2_post" %in% x_att$names){
              samples = x$admin2_post
          }else{
              samples = x$admin1_post
          }
        # USING SURVEYPREV CLASSES: direct est
        }else{
          if("res.admin1" %in% names(x)){
            domain.names <- x$res.admin1$admin1.name
            samples <- matrix(NA, nsim, length(domain.names))
            for(i in 1:dim(x$res.admin1)[1]){
              samples[, i] <- expit(rnorm(nsim, mean = x$res.admin1$direct.logit.est[i], 
                                          sd = (x$res.admin1$direct.logit.prec[i])^(-1/2)))
            }
          }else{
            domain.names <- x$res.admin2$admin2.name.full
            samples <- matrix(NA, nsim, length(domain.names))
            for(i in 1:dim(x$res.admin2)[1]){
              samples[, i] <- expit(rnorm(nsim, mean = x$res.admin2$direct.logit.est[i], 
                                          sd = (x$res.admin2$direct.logit.prec[i])^(-1/2)))
            }
          }
        }
        # samples is now a nsamp * nregion matrix
        samples.long <- data.frame(region.name = rep(domain.names, each = nrow(samples)), 
                                   value = as.numeric(samples))
        if("res.admin2" %in% names(x)){
          upper <- x$res.admin2[, c("admin2.name.full", "admin1.name")]
          upper$admin2.name.short <- NA
            for(i in 1:dim(upper)[1]){
              k <- nchar(as.character(upper$admin1.name[i]))
              upper$admin2.name.short[i] <- substr(upper$admin2.name.full[i], 
                                 start = k+2, 
                                 stop = nchar(upper$admin2.name.full[i]))
            }
          colnames(upper) <- c("region.name", "group.name", "admin2.name.short")
          samples.long <- dplyr::left_join(samples.long, upper)
        }

        n.levels <- dim(samples)[2]
        ridge.max <- max(samples.long$value)*1.03
        if(ridge.max>0.95){ridge.max=1}
        ridge.min <- min(samples.long$value)*0.97
        if(ridge.min<0.05){ridge.min=0}

        g <- ggplot2::ggplot(samples.long) +
            aes(x = value, y = region.name) + 
            ggridges::geom_density_ridges_gradient(aes(fill = ggplot2::after_stat(x))) +
            ggplot2::scale_fill_viridis_c(lim = c(ridge.min,ridge.max), direction = direction) + ggplot2::theme(legend.position = "none") + xlab("") + ylab("")
        if("group.name" %in% colnames(samples.long)){
          g <- g + aes(y = admin2.name.short) + ggplot2::facet_wrap(~group.name, scale = "free_y")
        }
        return(g)
       }

    
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
          if(is.null(strata.plot)){
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
          year.label <- tmp 
        
          timelabel.yearly <- year.label
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
          year.label <- x$year.label
          if(is.yearly){
            timelabel.yearly <- c(x$year.range[1] : x$year.range[2], year.label)
          }else{
            timelabel.yearly <- year.label
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
          results$is.yearly <- !(results$years %in% year.label)
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
      if(is.null(year.plot)){
        year.plot <- year.label 
      }
      # plot calculated density
      if(by.year && is.density){
         g <- ggplot2::ggplot(subset(results.plot, years %in% year.plot), ggplot2::aes(x = x, y = region, height = y, fill = ..x..)) 
      }else if(is.density){
         results.plot$years <- factor(results.plot$years, levels = rev(timelabel.yearly))
         g <- ggplot2::ggplot(subset(results.plot, years %in% year.plot), ggplot2::aes(x = x, y = years, height = y, fill = ..x..)) 
      }
      if(is.density) g <- g + ggridges::geom_density_ridges_gradient(stat="identity", alpha = 0.5, linewidth = linewidth) 

      # plot draws
      if(by.year && !is.density){
        g <- ggplot2::ggplot(subset(results.plot, years %in% year.plot), ggplot2::aes(x = x, y = region)) 
      }else if(!is.density){
        results.plot$years <- factor(results.plot$years, levels = rev(timelabel.yearly))
        g <- ggplot2::ggplot(subset(results.plot, years %in% year.plot), ggplot2::aes(x = x, y = years)) 
      }
      if(!is.density) g <- g + ggridges::geom_density_ridges_gradient(ggplot2::aes(fill = ..x..), scale = scale, alpha = 0.5, linewidth = linewidth)

      g <- g + ggplot2::scale_fill_viridis_c(option = "D", direction = direction) + ggplot2::theme_bw()  + ggplot2::ylab("") + ggplot2::theme(legend.position = 'none') + ggplot2::xlab("")
      if(by.year){
        if(length(year.plot) > 1) g <- g + ggplot2::facet_wrap(~years, ncol = ncol)
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
