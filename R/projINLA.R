#' Function to obtain projected estimates from INLA for each time and region.
#' 
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
#' }
#' 
#' 
#' @export
projINLA <- function(data, inla_mod, years, geo, newyear = "15-19", quantiles = c(0.025, 0.5, 0.975)) {
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
