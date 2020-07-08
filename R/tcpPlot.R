#' Discrete-color maps based on the True Classification Probabilities
#' 
#' 
#'
#' @param draws a posterior draw object from \code{\link{getSmoothed}}
#' @param geo SpatialPolygonsDataFrame object for the map
#' @param by.geo variable name specifying region names in geo
#' @param year_plot vector of year string vector to be plotted.
#' @param ncol number of columns in the output figure.
#' @param per1000 logical indicator to multiply results by 1000.
#' @param thresholds a vector of thresholds (on the mortality scale) defining the discrete color scale of the maps. 
#' @param intervals number of quantile intervals defining the discrete color scale of the maps. Required when thresholds are not specified.
#' @param size.title  a numerical value giving the amount by which the plot title should be magnified relative to the default.
#' @param legend.label Label for the color legend.
#' @param border color of the border
#' @param size size of the border
#' @return a list of True Classification Probability (TCP) tables, a list of individual spplot maps, and a gridded array of all maps.
#' @import data.table
#' @import RColorBrewer
#' @importFrom grDevices colorRampPalette
#' @importFrom methods is
#'
#' @author Tracy Qi Dong, Zehang Richard Li
#' @references Tracy Qi Dong, and Jon Wakefield. (2020) \emph{Modeling and presentation of vaccination coverage estimates using data from household surveys.} arXiv preprint arXiv:2004.03127.
#' @examples
#' \dontrun{
#' library(dplyr)
#' data(DemoData)
#' # Create dataset of counts, unstratified
#' counts.all <- NULL
#' for(i in 1:length(DemoData)){
#'   counts <- getCounts(DemoData[[i]][, c("clustid", "time", "age", "died",
#'                                         "region")],
#'             variables = 'died', by = c("age", "clustid", "region", 
#'                                          "time"))
#'   counts <- counts %>% mutate(cluster = clustid, years = time, Y=died)
#'   counts$strata <- NA
#'   counts$survey <- names(DemoData)[i] 
#'   counts.all <- rbind(counts.all, counts)
#' }
#' 
#' # fit cluster-level model on the periods
#' periods <- levels(DemoData[[1]]$time)
#' fit <- smoothCluster(data = counts.all, 
#'       Amat = DemoMap$Amat, 
#'       time.model = "rw2", 
#'       st.time.model = "rw1",
#'       strata.time.effect =  TRUE, 
#'       survey.effect = TRUE,
#'       family = "betabinomial",
#'       year_label = c(periods, "15-19"))
#' est <- getSmoothed(fit, nsim = 1000, save.draws=TRUE)
#' 
#' tcp <- tcpPlot(est, DemoMap$geo, by.geo = "REGNAME", interval = 3, year_plot = periods) 
#' tcp$g
#' }
#' 

#' @export
tcpPlot <- function(draws, geo, by.geo = NULL, year_plot = NULL, ncol = 4, per1000 = FALSE, thresholds = NULL, intervals = 3, size.title = 0.7, legend.label = NULL, border = "gray20", size = 0.5){

  grp_val <- long <- lat <- group <- NA

  if(is(draws, "data.frame") || (is(draws, "list") && !is.null(draws$fit))){
    stop("TCP plot has not been implemented with smoothed direct estimates for now...")

  }else if(is.null(draws$draws.est.overall)){
     stop("Posterior draws not found. Please rerun getSmoothed() with save.draws = TRUE.")
 
  }else{
      if (!isTRUE(requireNamespace("INLA", quietly = TRUE))) {
        stop("You need to install the packages 'INLA'. Please run in your R terminal:\n  install.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/stable'), dep=TRUE)")
      }
      if (isTRUE(requireNamespace("INLA", quietly = TRUE))) {
        if (!is.element("INLA", (.packages()))) {
          attachNamespace("INLA")
        }
        
        geo@data$region_names <- as.character(geo@data[, by.geo])
        region_names <- geo@data$region_names
        region_nums <- 1:length(region_names)
        n_admin <- length(region_names)
        draws.plot <- draws$draws.est.overall

        tmp <- NULL
        # draw.est is ordered
        for(i in 1:length(draws.plot)) {
        if(draws.plot[[i]]$years %in% tmp) next
            tmp <- c(tmp, draws.plot[[i]]$years)
        }
        timelabel.yearly <- tmp 
        if(!is.null(year_plot)) timelabel.yearly <- year_plot
        
        n_years <- length(timelabel.yearly)
        
        n_samp <- length(draws.plot[[1]]$draws)
        
        ##################
        
        postsamp_mt_list <- vector(mode = "list", length = n_years)
        postsamp_vt <- NULL
        
        for (i in 1:n_years){
          # i <- 1
          postsamp_mt_list[[i]]$years <- timelabel.yearly[i]
          
          postsamp_mt <- matrix(0, nrow = n_admin, ncol = n_samp)
          for(j in 1:n_admin){
            # j <- 1
            draw_est_j <- draws.plot[lapply(draws.plot, '[[',"region")==region_names[j]]
            draw_est_ij <- draw_est_j[lapply(draw_est_j, '[[', "years")==timelabel.yearly[i]]
            postsamp_mt[j, ] <- draw_est_ij[[1]]$draws
            postsamp_vt <- append(postsamp_vt, draw_est_ij[[1]]$draws)
          }
          postsamp_mt_list[[i]]$postsamp_mt <- postsamp_mt
        }
        
        ##############
        
        get_measure_dt <- function(postsamp_mt, pred_dt, grp_thresh){
          lower <- upper <- colid <- rowid <- value <- J <- grp <- TCP <- NA
          n_grp <- length(grp_thresh)-1
          grp_dt <- data.table::data.table(grp = 1:n_grp,
                                           lower = grp_thresh[1:(n_grp)],
                                           upper = grp_thresh[2:(n_grp+1)])
          n_area <- nrow(postsamp_mt)
          n_postsamp <- ncol(postsamp_mt)
          grp_cnt_mt <- matrix(0, nrow = n_area, ncol = n_grp)
          for (i in 1:n_grp){
            grp_cnt_mt[, i] <- apply(postsamp_mt, 1, 
                                     function(x){sum(x > grp_dt[i, lower] & x <= grp_dt[i, upper])})
          }
          DF <- data.frame(grp_cnt_mt)
          DT <- data.table::data.table(value = unlist(DF, use.names=FALSE), 
                                       colid = 1:nrow(DF), 
                                       rowid = rep(names(DF), each=nrow(DF)))
          data.table::setkey(DT, colid, value)
          grp_cnt_dt_max <- data.table::as.data.table(DF)
          grp_cnt_dt_max[, "grp"] <- DT[J(unique(colid), DT[J(unique(colid)), value, mult="last"]), rowid, mult="first"]
          grp_cnt_dt_max[, "TCP" := 0]
          for (i in 1:n_grp){
            idx_grp <- which(grp_cnt_dt_max[, grp] == paste0("X", i))
            grp_cnt_dt_max[idx_grp, "TCP"] <- grp_cnt_dt_max[idx_grp, paste0("X", i), with = F]/n_postsamp
          }
          pred_dt[, "grp"] <- grp_cnt_dt_max[, grp]
          pred_dt[, "TCP"] <- grp_cnt_dt_max[, TCP]
          return(pred_dt)
        }
        
        if (is.null(thresholds)){
          K <- intervals
          intv <- 1/K
          quant_vt <- seq(0, 1, intv)
          quant_val_vt <- (quant_vt[1:K]+quant_vt[2:(K+1)])/2
          L_vt <- quantile(postsamp_vt, quant_vt)
          L_vt[1] <- quantile(postsamp_vt, 0.01)
          L_vt[K+1] <- quantile(postsamp_vt, 0.99)
          L_val_vt <- quantile(postsamp_vt, quant_val_vt)
        }else{
          # thresholds <- L_vt
          K <- length(thresholds)-1
          L_vt <- thresholds
          L_val_vt <- (thresholds[1:K]+thresholds[2:(K+1)])/2
        }
        
        L_dt <- data.table::data.table(grp = paste0("X", 1:K),
                                       grp_low = L_vt[1:K],
                                       grp_up = L_vt[2:(K+1)],
                                       grp_val = L_val_vt)
                           
        
        manual.col <- colorRampPalette(RColorBrewer::brewer.pal(8, "YlGnBu"))(length(L_dt$grp_val))
        color.match <- manual.col[1:length(L_dt$grp_val)]
        lookup_dt <- data.table::data.table(grp_val = sort(L_dt$grp_val),
                                            col = color.match)
        labelat <- sort(unique(c(L_dt$grp_low, 
                                 L_dt$grp_up)))
        if (per1000){
          labeltext <- format(round(labelat*1000, 1), nsmall = 1)
        }else{
          labeltext <- format(round(labelat, 3), nsmall = 3)
        }
        
        which.plot <- match(year_plot, timelabel.yearly)
        n_plot <- length(which.plot)
        map_list <- TCP_list <- vector(mode = "list", length = n_plot)
        toplot <- NULL
        for (i in 1:n_plot){
          # i <- 1
          TCP_list[[i]]$years <- timelabel.yearly[which.plot[i]]
          measure_dt <- get_measure_dt(postsamp_mt = postsamp_mt_list[[which.plot[i]]]$postsamp_mt,
                                       pred_dt = data.table::data.table(region_names = region_names),
                                       grp_thresh = L_vt)
          measure_dt <- merge(measure_dt, L_dt, by = "grp", all.x = T)
          
          TCP_list[[i]]$TCPtable <- data.frame(region_names = measure_dt$region_names,
                                               TCP = measure_dt$TCP,
                                               interval = measure_dt$grp,
                                               interval_lower = measure_dt$grp_low,
                                               interval_upper = measure_dt$grp_up)
          
          ###########################################################
          ## Using SPPLOT
          ###########################################################
          shp_plot <- sp::merge(geo, measure_dt, by = "region_names")
          shp_plot <- sp::merge(shp_plot, lookup_dt, by = "grp_val")
          shp_plot$grp_val <- as.factor(shp_plot$grp_val)
          col_regions <- as.vector(lookup_dt[grp_val %in% shp_plot$grp_val, col])
          
          DiscreteMap <- sp::spplot(shp_plot, zcol = "grp_val",
                                    main = list(label = paste0("Year ", timelabel.yearly[which.plot[i]], ", ATCP = ", format(round(mean(measure_dt$TCP), 2), nsmall = 2)),
                                                cex = size.title),
                                    xlab = "", ylab = "",
                                    col.regions = col_regions,
                                    colorkey = list(col = color.match,
                                                    at = labelat,
                                                    labels = list(at = labelat, labels = labeltext)))
                                
          map_list[[i]] <- DiscreteMap

          ###########################################################
          ## Using ggplot? 
          ###########################################################
          subplot <- data.frame(measure_dt)
          subplot$label <- paste0("Year ", timelabel.yearly[which.plot[i]], 
                                ", ATCP = ", format(round(mean(measure_dt$TCP), 2)))
          toplot <- rbind(toplot, subplot)
        }

        ###########################################################
        ## Using ggplot? 
        ###########################################################
        has.coord <- !is.na(sp::proj4string(geo))
        geo <- ggplot2::fortify(geo, region = by.geo)
        geo2 <- merge(geo, toplot, by = "id", by.y = "region_names")
        g <- ggplot2::ggplot(geo2)
        g <- g + ggplot2::geom_polygon(ggplot2::aes(x = long, y = lat, 
                group = group, fill = grp_val), color = border, size = size)
        g <- g + ggplot2::facet_wrap(~label, ncol = ncol)

        # Create palette
        color_cont <- NULL 
        for(i in 1:dim(lookup_dt)[1]) color_cont <- c(color_cont, rep(as.character(lookup_dt[i, 2]), round((labelat[i+1] - labelat[i]) * 1e5)))
        g <- g + ggplot2::scale_fill_gradientn(legend.label, colours = color_cont, limits=range(labelat),  breaks = labelat, labels = labeltext)

        # Make legends taller
        # panel_height = ggplot2::unit(1,"npc") - sum(ggplot2::ggplotGrob(g)[["heights"]][-3]) - ggplot2::unit(1,"line")
        if(is.null(ncol)){
          nrow.panel <- ggplot2::wrap_dims(length(unique(geo2$label)))[1]
        }else{
          nrow.panel <- ceiling(length(unique(geo2$label)) / ncol)
        } 
        # panel_height <- ggplot2::unit(0.9/as.numeric(nrow.panel),"npc") 
        g <- g + ggplot2::guides(fill= ggplot2::guide_colorbar(title.theme = ggplot2::element_blank(), barheight=ggplot2::rel(14))) 

        if(has.coord) g <- g + ggplot2::coord_map()

        # Clean themes
        g <- g + ggplot2::theme_bw() + ggplot2::theme(legend.title=ggplot2::element_text(size=ggplot2::rel(size.title)), axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(), axis.ticks.x=ggplot2::element_blank(), axis.title.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.ticks.y=ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())
        
        return(list(TCP_list = TCP_list, map_list = map_list, g = g, labelat = labelat, labeltext = labeltext, lookup_dt = lookup_dt))
        
      }
    
  }
}