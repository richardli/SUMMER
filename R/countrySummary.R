#' Obtain the Horvitz-Thompson direct estimates and standard errors using delta method for a single survey.
#'
#' 
#' @param births A matrix child-month data from \code{\link{getBirths}}
#' @param years String vector of the year intervals used
#' @param regionVar Variable name for region in the input births data.
#' @param timeVar Variable name for the time period indicator in the input births data.
#' @param ageVar Variable name for age group. This variable need to be in the form of "a-b" where a and b are both ages in months. For example, "1-11" means age between 1 and 11 months, including both end points. An exception is age less than one month can be represented by "0" or "0-0".
#' @param weightsVar Variable name for sampling weights, typically 'v005'
#' @param clusterVar Variable name for cluster, typically '~v001 + v002'
#' @param geo.recode The recode matrix to be used if region name is not consistent across different surveys. See \code{\link{ChangeRegion}}.
#' @param national.only Logical indicator to obtain only the national estimates
#'
#' @return a matrix of period-region summary of the Horvitz-Thompson direct estimates by region and time period specified in the argument, the standard errors using delta method for a single survey, the 95\% confidence interval, and the logit of the estimates.
#' @seealso \code{\link{countrySummary_mult}}
#' @examples
#' \dontrun{
#' data(DemoData)
#' years <- c("85-89", "90-94", "95-99", "00-04", "05-09", "10-14")
#' u5m <- countrySummary(births = DemoData[[1]],  years = years, 
#' regionVar = "region", timeVar = "time", clusterVar = "~clustid+id", 
#' ageVar = "age", weightsVar = "weights", geo.recode = NULL)
#' }
#' @export
countrySummary <- function(births, years, regionVar = "region", timeVar = "time", clusterVar = "~v001+v002",  ageVar = "age", weightsVar = "v005", geo.recode = NULL, national.only = FALSE) {
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
        ns <- c(1, 11, 12, 12, 12, 12)
    }
    labels <- as.character(unique(births$age0))
    ns <- rep(1, length(labels))
    for(i in 1:length(labels)){
        if(labels[i] == "0"){
            ns[i] <- 1
            next
        }
        tmp <- as.numeric(strsplit(as.character(labels[i]), "-")[[1]])
        ns[i] <- tmp[2] - tmp[1] + 1
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
        stop("Cluster not defined")
    }

    options(survey.lonely.psu = "adjust")
    my.svydesign <- survey::svydesign(ids = stats::formula(clusterVar), strata = ~strata, nest = T, weights = ~weights0, 
        data = births)
    
    # get region list, sorted alphabetically
    regions_list <- as.character(sort(names(table(births$region0))[as.vector(table(births$region0) != 0)]))
    regions_num <- 1:length(regions_list)
    
    # add 'All' to all regions
    if(national.only){
        regions_list <- c("All")
        regions_num <- c(0)
    }else{
        regions_list <- c("All", regions_list)
        regions_num <- c(0, regions_num)        
    }
    
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
        } else if(length(unique(tmp$variables$age0)) > 1){
            glm.ob <- survey::svyglm(died ~ (-1) + factor(age0), design = tmp, family = stats::quasibinomial, maxit = 50)
            if(dim(summary(glm.ob)$coef)[1] < length(labels)){
                bins.nodata <- length(labels) - dim(summary(glm.ob)$coef)[1] 
                if(bins.nodata >= length(ns)/2){
                    warning(paste0(which.area, " ", which.time, " has no observation in more than half of the age bins, set to NA\n"),immediate. = TRUE)
                    return(rep(NA, 5))
                }
                # This can only happen for the last one or several bins, since person-month are cumulative
                ns.comb <- ns
                ns.comb[length(ns) - bins.nodata] <- sum(ns[c(length(ns) - bins.nodata) : length(ns)])
                warning(paste0(which.area, " ", which.time, " has no observations in ", bins.nodata, " age bins. They are collapsed with previous bins\n"),immediate. = TRUE)
                return(get.est.withNA(glm.ob, labels, ns.comb))
            }
            return(get.est.withNA(glm.ob, labels, ns))
        } else {
            glm.ob <- survey::svyglm(died ~ 1, design = tmp, family = stats::quasibinomial, maxit = 50)
            var.est <- stats::vcov(glm.ob)
            mean <- expit(summary(glm.ob)$coefficient[1])
            lims <- expit(logit(mean) + stats::qnorm(c(0.025, 0.975)) * sqrt(c(var.est)))

            ## Alternative calculation
            # p.i <- survey::svymean(~died, design = tmp)
            # var.i <- as.numeric(vcov(p.i)) 
            # p.i <- as.numeric(p.i)   
            # ht <- log(p.i/(1-p.i))
            # ht.v <- var.i/(p.i^2*(1-p.i)^2)
            # ht.prec <- 1/ht.v
           return(c(mean, lims, logit(mean), var.est))

        }
    }
    
    # updated helper function: get.est notes: the original function only works when all age group exist, if some age groups are
    # missing, the dimension from covariance matrix will not match the ns vector note: not vary satisfying modification yet since
    # the 'factor(ageGrpD)'
    
    get.est.withNA <- function(glm.ob, labels, ns) {
        ## initialize with the full covariance matrix
        # K <- dim(summary(glm.ob)$coef)[1]
        K <- length(labels)
        V <- matrix(0, K, K)
        betas <- rep(NA, K)
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
        
        # ns <- c(1, 11, 12, 12, 12, 12)
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
    results$survey <- NA
    results$logit.prec <- 1/results$var.est
    results <- results[, c("region", "years", "u5m", "lower", "upper", "logit.est", "var.est", "region_num", "survey", "logit.prec")]
     
    return(results)
}


