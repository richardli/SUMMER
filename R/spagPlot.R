#' Makes spaghetti plot.
#'
#' @param countryname Country name as a string
#' @param results output from \code{\link{projINLA}}
#' @param geo \code{geo} output from \code{\link{read_shape}}
#' @param countrysum output from \code{\link{countrySummary_mult}}
#' @param inlamod output from \code{\link{fitINLA}}
#' @examples
#' \dontrun{
#' data(Uganda)
#' data(UgandaMap)
#' geo <- UgandaMap$geo
#' mat <- UgandaMap$Amat
#' years <- c("85-89", "90-94", "95-99", "00-04", "05-09", "10-14")
#' 
#' # Get direct estimates
#' u5m <- countrySummary_mult(births = Uganda, years = years, idVar = "id", 
#' regionVar = "region", timeVar = "time", clusterVar = "~clustid+id", 
#' ageVar = "age", weightsVar = "weights", geo.recode = NULL)
#' 
#' # Get hyper priors
#' priors <- simhyper(R = 2, nsamp = 1e+05, nsamp.check = 5000, Amat = mat)
#' 
#' # Fit INLA models
#' data <- data[data$region %in% c("central","eastern","northern","western"),]
#' inla_model <- fitINLA(data = data, geo = geo, Amat = mat, year_names = years, priors = priors)
#' 
#' # Projection
#' surveylabel <- paste0("DHS ", unique(data$surveyYears)) 
#' results_rw2 <- projINLA(data = data, inla_mod = inla_model, years = years, geo = geo, 
#'                      newyear = "15-19", quantiles = c(0.025,0.5,0.975))
#' 
#' # Plot results
#' spagPlot(countryname = "Uganda", results = results_rw2, geo = geo, 
#'            countrysum = data, inlamod = inla_model)
#' }
#' @export
spagPlot <- function(countryname, results, geo, countrysum, inlamod) {
    out <- list(countryname = countryname, results.rw2 = results, geo = geo, 
                data.HT = countrysum, model.rw2 = inlamod)
    surveylabel <- paste0("DHS ", unique(countrysum$surveyYears))
    out$data.HT$survey.label <- surveylabel[out$data.HT$survey]
  
    years <- names(out$results.rw2)
    n.years <- length(years)
    n.area <- dim(out$geo)[1]
    areasCap <- out$geo@data$DHSREGEN
    if (is.null(areasCap)) {
        areasCap <- out$geo@data$NAME_final
    }
    areas.smooth <- rownames(out$results.rw2[[1]])

    
    plot.res <- expand.grid(District = areas.smooth, Year = years)
    plot.res$sd <- plot.res$mean <- plot.res$q90 <- plot.res$q10 <- plot.res$q975 <- plot.res$q025 <- plot.res$med <- NA
    
    
    for (i in 1:(n.years)) {
        plot.res$med[plot.res$Year == years[i]] <- apply(out$results.rw2[[i]], 1, function(x) {
            stats::median(expit(x))
        })
        plot.res$q025[plot.res$Year == years[i]] <- apply(out$results.rw2[[i]], 1, function(x) {
            stats::quantile(expit(x), 0.025)
        })
        plot.res$q975[plot.res$Year == years[i]] <- apply(out$results.rw2[[i]], 1, function(x) {
            stats::quantile(expit(x), 0.975)
        })
        plot.res$q10[plot.res$Year == years[i]] <- apply(out$results.rw2[[i]], 1, function(x) {
            stats::quantile(expit(x), 0.1)
        })
        plot.res$q90[plot.res$Year == years[i]] <- apply(out$results.rw2[[i]], 1, function(x) {
            stats::quantile(expit(x), 0.9)
        })
        plot.res$mean[plot.res$Year == years[i]] <- apply(out$results.rw2[[i]], 1, function(x) {
            mean(expit(x))
        })
        plot.res$sd[plot.res$Year == years[i]] <- apply(out$results.rw2[[i]], 1, function(x) {
            stats::sd(expit(x))
        })
    }
    
    
    plot.res$ratio <- plot.res$q975/plot.res$q025
    plot.res$unc_ratio <- plot.res$ratio > 3
    plot.res$unc_ratio0 <- plot.res$ratio >= 2 & plot.res$ratio < 2.25
    plot.res$unc_ratio1 <- plot.res$ratio >= 2.25 & plot.res$ratio < 2.5
    plot.res$unc_ratio2 <- plot.res$ratio >= 2.5
    
    med.palette <- RColorBrewer::brewer.pal(n = 9, name = "Purples")
    med.int <- classInt::classIntervals(round(plot.res$med, 3), n = 9, style = "jenks")
    med.col <- classInt::findColours(med.int, med.palette)
    plot.res$med.col <- med.col
    

    time <- utils::head(years,-1)
    
    cols <- grDevices::rainbow(n.area)
    ymin <- 0
    ymax <- 0.3
    ymax2 <- max(plot.res$q975)
    if (ymax2 > ymax) 
        ymax <- ymax2
    # graphics::par(mfrow = c(1, 1))
    # graphics::plot(NA, xlim=c(1,length(time)), xaxt="n",ylim = c(ymin, ymax),xlab="Year",ylab = "Smooth Estimate")
    # graphics::axis(side = 1, at = c(1:length(time)),labels = time)
    # for (i in 1:n.area) {
    #     jit.time <- jitter(1:length(time),.25)
    #     tmp <- plot.res[plot.res$District == areas.smooth[i] & plot.res$Year %in% head(years,-1), ]
    #     graphics::lines(jit.time, tmp$med, lwd = 2, col = cols[i])
    #     
    #     for (k in 1:(n.years - 1)) {
    #         graphics::segments(jit.time[k], tmp$q025[k], jit.time[k], tmp$q975[k], col = cols[i])
    #     }
    # }
    # graphics::legend("topright", cex = 0.75, bty = "n", lwd = 2, col = cols, legend = areasCap)
    ### Consider for later to make nicer ggplots
    df <- data.frame(Time= rep(time,each=n.area), plot.res[plot.res$Year %in% utils::head(years,-1),])
    df$Time <- factor(df$Time, levels(df$Time)[order(utils::head(levels(plot.res$Year),-1))])
    levels(df$District) <- levels(areasCap)
    pd <- ggplot2::position_dodge(width=0.3)
    # fix for global variable issue
    Time <- NULL; District <- NULL; med <- NULL; q025 <- NULL; q975 <- NULL
    ggplot2::ggplot(df, ggplot2::aes(x=Time,y=med,group=District,colour=District)) +
      ggplot2::geom_line(ggplot2::aes(lty=District),lwd=1) +
      ggplot2::ylim(ymin,ymax) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=q025, ymax=q975), width=.3,position=pd, show.legend=FALSE) +
      #geom_point() +
      ggplot2::labs(x = "Time", y = "Smooth Estimate", title = countryname, colour = "Region", lty="Region") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) 
}
