#' Obtain the Horvitz-Thompson direct estimates and standard errors using delta method for a single survey.
#'
#' 
#' @param births A matrix child-month data from \code{\link{getBirths}}
#' @param years String vector of the year intervals used
#' @param idVar Variable name for ID, typically 'v002'
#' @param regionVar Variable name for region, typically 'v024', for older surveys might be 'v101'
#' @param timeVar Variable name for time, typically 'per5'
#' @param ageVar Variable name for age group, default assumes the variable is called 'ageGrpD'
#' @param weightsVar Variable name for sampling weights, typically 'v005'
#' @param clusterVar Variable name for cluster, typically '~v001 + v002'
#' @param geo.recode The recode matrix to be used if region name is not consistent across different surveys. See \code{\link{ChangeRegion}}.
#'
#' @return a matrix of period-region summary of the Horvitz-Thompson direct estimates, the standard errors using delta method for a single survey, the 95\% confidence interval, and the logit of the estimates.
#' @seealso \code{\link{countrySummary_mult}}
#' @examples
#' \dontrun{
#' data(DemoData)
#' years <- c("85-89", "90-94", "95-99", "00-04", "05-09", "10-14")
#' u5m <- countrySummary(births = DemoData[[1]],  years = years, idVar = "id", 
#' regionVar = "region", timeVar = "time", clusterVar = "~clustid+id", 
#' ageVar = "age", weightsVar = "weights", geo.recode = NULL)
#' }
#' @export
countrySummary <- function(births, years, idVar = "v002", regionVar = "region", timeVar = "per5", clusterVar = "~v001+v002",
                           ageVar = "ageGrpD", weightsVar = "v005", geo.recode = NULL) {
    # check all elements are provided
    if (is.null(births)) {
        stop("No births file specified!")
    }
    if (is.null(years)) {
        stop("Names of the 5-year intervals not provided!")
    }
    
    if (!regionVar %in% colnames(births)) {
        if ("v101" %in% colnames(births)) {
            colnames(births)[which(colnames(births) == "v101")] <- regionVar
            warning("region variable not defined: using v101", immediate. = TRUE)
        } else if ("v024" %in% colnames(births)) {
            colnames(births)[which(colnames(births) == "v024")] <- regionVar
            warning("region variable not defined: using v024", immediate. = TRUE)
        } else {
            stop("region variable not defined, and no v101 or v024!")
        }
    }
    
    # recode geo information
    if (!is.null(geo.recode)) {
        births <- ChangeRegion(births, Bmat = geo.recode, regionVar = regionVar)
    }
    
    # create new variables
    births$region0 <- births[, regionVar]
    births$id0 <- births[, idVar]
    births$weights0 <- births[, weightsVar]
    births$time0 <- births[, timeVar]
    births$age0 <- births[, ageVar]
    
    # fix for input data issue
    if (sum(c(0, 1, 12, 24, 36, 48) %in% births$age0) == 6) {
        births$age0[births$age0 == 0] <- "0"
        births$age0[births$age0 == 1] <- "1-11"
        births$age0[births$age0 == 12] <- "12-23"
        births$age0[births$age0 == 24] <- "24-35"
        births$age0[births$age0 == 36] <- "36-47"
        births$age0[births$age0 == 48] <- "48-59"
    }
    
   
    time_inconsistent <- which(!(births$time0 %in% years))
    if (length(time_inconsistent) > 0) {
        warning(paste("Name for time periods are inconsistent. Found the following levels in data:", unique(births$time0[time_inconsistent])), 
            immediate. = TRUE)
    }
    
    
    if (is.null(births$strata)) {
        stop("Strata not defined.")
    }
    if (is.null(clusterVar)){
        clusterVar <- paste0("~", idVar)
        warning(paste("Cluster not specified, use", clusterVar, "instead"), immediate. = TRUE)
    }

    options(survey.lonely.psu = "adjust")
    my.svydesign <- survey::svydesign(ids = stats::formula(clusterVar), strata = ~strata, nest = T, weights = ~weights0, 
        data = births)
    
    # get region list, sorted alphabetically
    regions_list <- as.character(sort(names(table(births$region0))[as.vector(table(births$region0) != 0)]))
    regions_num <- 1:length(regions_list)
    
    # add 'All' to all regions
    regions_list <- c("All", regions_list)
    regions_num <- c(0, regions_num)
    
    # create result data frame
    results <- data.frame(region = rep(regions_list, each = length(years)))
    results$region_num <- rep(regions_num, each = length(years))
    results$years <- rep(years, length(regions_list))
    # add empty variables
    results$var.est <- results$logit.est <- results$upper <- results$lower <- results$u5m <- NA
    
    # updated helper function: region.time.HT notes: the original function only works when the selected area-time combination has
    # data. Add a new line of codes so that when there is no data for selected combination, return a line of NA values updated
    # Aug 17, 2015: Enable area = 'All'
    
    region.time.HT.withNA <- function(which.area, which.time) {
        #cat(".")
        # Address visible binding note
        time0 <- NULL
        region0 <- NULL
        if (which.area == "All") {
            tmp <- subset(my.svydesign, (time0 == which.time))
        } else {
            tmp <- subset(my.svydesign, (time0 == which.time & region0 == as.character(which.area)))
        }
        if (dim(tmp)[1] == 0) {
            return(rep(NA, 5))
        } else if (sum(tmp$variables$died) == 0) {
            warning(paste0(which.area, " ", which.time, " has no death, set to NA\n"),immediate. = TRUE)
            return(rep(NA, 5))
        } else {
            glm.ob <- survey::svyglm(died ~ (-1) + factor(age0), design = tmp, family = stats::quasibinomial, maxit = 50)
            return(get.est.withNA(glm.ob))
        }
    }
    
    # updated helper function: get.est notes: the original function only works when all age group exist, if some age groups are
    # missing, the dimension from covariance matrix will not match the ns vector note: not vary satisfying modification yet since
    # the 'factor(ageGrpD)'
    
    get.est.withNA <- function(glm.ob) {
        ## initialize with the full covariance matrix
        V <- matrix(0, 6, 6)
        betas <- rep(NA, 6)
        labels <- c("0", "1-11", "12-23", "24-35", "36-47", "48-59")
        labels <- paste("factor(age0)", labels, sep = "")
        colnames(V) <- rownames(V) <- labels
        names(betas) <- labels
        # now get the regression covariance matrix
        V2 <- stats::vcov(glm.ob)
        # check if the columns match
        if (length(which(colnames(V2) %in% colnames(V) == FALSE)) > 0) {
            stop("Error for input age group names!")
        }
        # fill in the full V
        V[rownames(V2), colnames(V2)] <- V2
        
        ## same fill for betas under 5 child mortality for males in 00-02 ##
        betas2 <- summary(glm.ob)$coef[, 1]
        betas[names(betas2)] <- betas2
        
        # (removed) avoid NA by replacing with 0 This is wrong way to handle NA!  messed up the mean estimates
        # betas[which(is.na(betas))] <- 0
        
        ns <- c(1, 11, 12, 12, 12, 12)
        probs <- expit(betas)
        
        u5m.est <- (1 - prod((1 - probs)^ns, na.rm = TRUE))  #*1000  
        
        ## partial derivatives ##
        gamma <- prod((1 + exp(betas))^ns, na.rm = TRUE)
        derivatives <- (gamma)/(gamma - 1) * ns * expit(betas)
        
        ## now handle the NA in derivatives, it won't affect mean estimation
        derivatives[which(is.na(derivatives))] <- 0
        ## Items to return ##
        var.est <- t(derivatives) %*% V %*% derivatives
        lims <- logit(u5m.est) + stats::qnorm(c(0.025, 0.975)) * sqrt(c(var.est))
        return(c(u5m.est, expit(lims), logit(u5m.est), var.est))
    }
    
    # run for all combinations
    x <- mapply(region.time.HT.withNA, which.area = results$region, which.time = results$years)
    results[, 4:8] <- t(x)
     
    return(results)
}


