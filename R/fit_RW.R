#' Function to get RW2 output
#' 
#' See package vignette for usage details.
#' 
#'
#' @param data country summary data from \code{\link{countrySummary_mult}}
#' @param inla_mod output from \code{\link{fitINLA}}
#' @param years years string vector
#' @param geo geographic polygon object
#' @param newyear string of years for projection, defaults to \code{'15-19'}
#' @param quantiles quantiles desired, defaults to \code{c(0.025,0.5,0.975)}
#' 
#' @return Results from RW2 model fit, including projection.
#' 
#' @export
fit_RW <- function(data, inla_mod, years, geo, newyear = "15-19", quantiles = c(0.025, 0.5, 0.975)) {
    # surveylabel <- paste0('DHS ', unique(data$surveyYears)) timelabel <- years countrylabel <- countryname
    
    pdata <- data.frame(data)
    pdata$years <- factor(pdata$years, levels = years)
    # model.smoothed <- NULL
    
    newyears <- c(years, newyear)
    proj.index <- length(years) + 1
    
    # inla_proj <- NULL
    results.rw2 <- vector("list", length(newyears))
    
    for (i in 1:length(newyears)) {
        results.rw2[[i]] <- matrix(NA, length(geo$NAME_final), 1000)
        rownames(results.rw2[[i]]) <- geo$NAME_final
    }
    
    for (j in 1:length(geo$NAME_final)) {
        # target.num <- j
        target <- geo$NAME_final[j]
        proj.rw2 <- projINLA_multi(fitted = inla_mod, proj.time = newyears[proj.index], ntime = length(years), which.area = target, quantiles = quantiles, 
            return_raw = TRUE)
        
        for (i in 1:length(newyears)) {
            results.rw2[[i]][j, ] <- proj.rw2[i, ]
        }
    }
    
    names(results.rw2) <- newyears
    
    return(results.rw2)
}
