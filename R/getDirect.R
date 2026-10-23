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
#' @param Ntrials Variable for the total number of person-months if the input data (births) is in the compact form.
#' @param geo.recode The recode matrix to be used if region name is not consistent across different surveys. See \code{changeRegion}.
#' @param national.only Logical indicator to obtain only the national estimates
#' @param CI the desired confidence interval to calculate
#' @param var.fix Logical; if \code{TRUE}, we attempt a target-specific synthetic-PSU correction when Rule 2 detects no between-PSU hazard contrast. The default preserves the usual estimates and attributes.
#' @param refStrataVar Optional name of a reference-group column. With \code{NULL}, we pool all real observations across areas and sampling strata separately within each period and age band. An explicit column must have one nonmissing value per design stratum within each period. We never infer this column from \code{births$strata}.
#' @return A data frame of period-region direct estimates, confidence limits, logit estimates, logit variances and precisions, and a survey placeholder. When \code{var.fix = TRUE}, we also attach \code{fixed_areas}, a character vector of triggered regions (including \code{"All"} when triggered), and \code{variance.fix}, a data frame with one row per triggered region-period. These attributes include unresolved attempts; \code{fixed_areas} does not imply successful correction.
#' @details
#' With \code{var.fix = TRUE}, we apply variance fix if there is a single PSU or 
#' if every exposed PSU's weighted hazard  equals the area's weighted hazard in every exposed age band. 
#'
#' We construct references before converting compact death counts to proportions.
#' Within each reference group and period, \eqn{K} counts all real PSUs with
#' retained exposure, including PSUs without exposure in a particular age band.
#' We average each PSU's exposure-weighted mean sampling weight to obtain
#' \eqn{w^*}. For each age band, synthetic exposure is
#' \eqn{E^*_a = \sum_i w_i E_{ia}/(K w^*)}, and synthetic deaths are
#' \eqn{D^*_a = \sum_i w_i D_{ia}/(K w^*)}. 
#' Observations with missing area labels still contribute to references and
#' national estimates. Zero-exposure rows do not contribute to the correction.
#'
#' We add one synthetic PSU in each design stratum present in a triggered target,
#' using that stratum's reference group. 
#' We retain the zero-death guard after augmentation. Targets without data remain
#' missing
#' @seealso \code{\link{getDirectList}}
#' @author Zehang Richard Li, Qianyu Dong, Bryan Martin, Laina Mercer
#' @references Li, Z., Hsiao, Y., Godwin, J., Martin, B. D., Wakefield, J., Clark, S. J., & with support from the United Nations Inter-agency Group for Child Mortality Estimation and its technical advisory group. (2019). \emph{Changes in the spatial distribution of the under-five mortality rate: Small-area analysis of 122 DHS surveys in 262 subregions of 35 countries in Africa.} PloS one, 14(1), e0210645.
#' @references Mercer, L. D., Wakefield, J., Pantazis, A., Lutambi, A. M., Masanja, H., & Clark, S. (2015). \emph{Space-time smoothing of complex survey data: small area estimation for child mortality.} The annals of applied statistics, 9(4), 1889.
#' @examples
#' \dontrun{
#' data(DemoData)
#' years <- c("85-89", "90-94", "95-99", "00-04", "05-09", "10-14")
#' mean <- getDirect(births = DemoData[[1]],  years = years, 
#' regionVar = "region", timeVar = "time", clusterVar = "~clustid+id", 
#' ageVar = "age", weightsVar = "weights", geo.recode = NULL)
#' }
#' @export
getDirect <- function(births, years, regionVar = "region", timeVar = "time", clusterVar = "~v001+v002",  ageVar = "age", weightsVar = "v005",  Ntrials = NULL, geo.recode = NULL, national.only = FALSE, CI = 0.95, var.fix = FALSE, refStrataVar = NULL) {
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
            message("region variable not defined: using v101")
        } else if ("v024" %in% colnames(births)) {
            colnames(births)[which(colnames(births) == "v024")] <- regionVar
            message("region variable not defined: using v024")
        } else {
            stop("region variable not defined, and no v101 or v024!")
        }
    }
    
    # recode geo information
    if (!is.null(geo.recode)) {
        births <- changeRegion(births, Bmat = geo.recode, regionVar = regionVar)
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
    }else{
        births$strata <- factor(births$strata)
    }
    if (is.null(clusterVar)){
        stop("Cluster not defined")
    }
    fix <- NULL
    audit <- .u5mr_fix_audit()
    if (isTRUE(var.fix)) {
        # Input: births with death counts, cluster formula, and optional exposure/reference column names (Ntrials/refStrataVar).
        # Output: fix stores compact counts/references, cluster design, retained-row mask, regions, target row indexes and phantom ID prefix/strata.
        fix <- .prepare_u5mr_fix(births, stats::formula(clusterVar), Ntrials, refStrataVar)
        # We omit zero-exposure rows before dividing death counts by exposure.
        births <- births[fix$retain, , drop = FALSE]
    }
    if(!is.null(Ntrials)){
        births$n <- births[,Ntrials]
        births$died <- births$died / births$n
        message("Compact input is used.")
    }
    cluster.formula <- stats::formula(clusterVar)
    cluster.vars <- all.vars(cluster.formula)
    compact.by <- unique(c(cluster.vars, "region0", "time0", "age0", "weights0", "strata", "died"))
    # Compact identical person-months using frequency weights, preserving age order.
    if (is.null(Ntrials) && length(cluster.vars) > 0 &&
        all(attr(stats::terms(cluster.formula), "term.labels") %in% cluster.vars) &&
        all(compact.by %in% names(births)) &&
        !any(cluster.vars %in% c("died", "n")) &&
        all(births$died %in% c(0, 1))) {
        compact <- data.table::as.data.table(births)[,
            list(n = .N), by = compact.by]
        # Retain the original boundary fit if an age-domain has no weighted survivors.
        if (nrow(compact) < nrow(births) &&
            all(compact[, any(get("died") == 0 & get("weights0") > 0, na.rm = TRUE),
                        by = c("region0", "time0", "age0")]$V1)) {
            births <- as.data.frame(compact)
            Ntrials <- "n"
        }
    }
    options(survey.lonely.psu = "adjust")
    # We retain the original design for untriggered targets. If no exposure
    # remains, the correction path returns missing estimates without fitting.
    my.svydesign <- if (isTRUE(var.fix) && !nrow(births)) NULL else survey::svydesign(ids = cluster.formula, strata = ~strata, nest = T, weights = ~weights0,
        data = births)
    
    # get region list, sorted alphabetically
    regions_list <- as.character(sort(names(table(births$region0))[as.vector(table(births$region0) != 0)]))
    # We retain areas whose rows had zero exposure. seq_along also permits a
    # national-only output when all real observations have missing area labels.
    if (isTRUE(var.fix)) regions_list <- fix$regions
    regions_num <- if (isTRUE(var.fix)) seq_along(regions_list) else 1:length(regions_list)
    
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
    results$var.est <- results$logit.est <- results$upper <- results$lower <- results$mean <- NA

    if (isTRUE(var.fix)) {
        targets <- vector("list", nrow(results))
        phantoms <- list()
        augmented.n <- nrow(fix$data)
        for (i in seq_len(nrow(results))) {
            # Row numbers index compact counts in their original order. Rule 2
            # uses only the real observations in this area-period.
            rows <- .u5mr_target(fix, results$region[i], results$years[i], rows.only = TRUE)
            target <- fix$data[rows, , drop = FALSE]
            targets[[i]] <- list(rows = rows)
            if (.u5mr_target_rule2(target)) {
                # Input: real target counts and area-period labels. Output: phantom death/exposure
                # counts, an audit row, and whether the augmented target still satisfies Rule 2.
                augmented <- .u5mr_augment(fix, target, results$region[i], results$years[i])
                phantoms[[length(phantoms) + 1L]] <- augmented$data
                audit <- rbind(audit, augmented$audit)
                targets[[i]] <- list(rows = c(rows, augmented.n + seq_len(nrow(augmented$data))),
                    audit.row = nrow(audit), degenerate = augmented$degenerate,
                    strata = unique(augmented$data$strata))
                augmented.n <- augmented.n + nrow(augmented$data)
            }
        }
        if (length(phantoms)) {
            augmented <- do.call(rbind, c(list(fix$data), phantoms))
            augmented$died <- augmented$died / augmented$n
            augmented.design <- survey::svydesign(ids = fix$formula, strata = ~strata,
                nest = TRUE, weights = ~weights0, data = augmented)
            # Cache original first-stage PSU counts for strata without usable
            # phantom exposure in a particular corrected target.
            original.psus <- my.svydesign$fpc$sampsize[
                match(fix$strata, as.character(my.svydesign$variables$strata)), 1L]
        }
    }
    
    # updated helper function: region.time.HT notes: the original function only works when the selected area-time combination has
    # data. Add a new line of codes so that when there is no data for selected combination, return a line of NA values updated
    # Aug 17, 2015: Enable area = 'All'
    
    region.time.HT.withNA <- function(which.area, which.time, target.index) {
        #cat(".")
        # Address visible binding note
        time0 <- NULL
        region0 <- NULL
        # We start each target with the unaugmented design and its input format.
        design <- my.svydesign
        trials <- !is.null(Ntrials)
        audit.row <- NULL
        augmented <- NULL
        # We record the outcome on every return from a triggered fit. A target
        # still satisfying Rule 2 is unresolved even if roundoff gives var > 0.
        finish <- function(value) {
            if (!is.null(audit.row)) {
                reason <- if (augmented$degenerate) "degenerate_reference" else
                    if (!is.finite(value[4])) "nonfinite_estimate" else
                    if (!is.finite(value[5]) || value[5] <= 0) "nonpositive_variance" else "augmented"
                audit$status[audit.row] <<- if (reason == "augmented") "fixed" else "unresolved"
                audit$reason[audit.row] <<- reason
            }
            value
        }
        if (isTRUE(var.fix)) {
            augmented <- targets[[target.index]]
            if (!length(augmented$rows)) return(finish(rep(NA, 5)))
            audit.row <- augmented$audit.row
            if (!is.null(audit.row)) {
                # Select this target's real rows and its own phantom rows,
                # including for national targets, from the shared design.
                tmp <- augmented.design[augmented$rows, ]
                # Another target may have added a phantom in a stratum that
                # has no usable reference exposure here. Restore its real count.
                missing.phantom <- !tmp$variables$strata %in% augmented$strata
                tmp$fpc$sampsize[missing.phantom, 1L] <- original.psus[
                    match(tmp$variables$strata[missing.phantom], fix$strata)]
                # We use n as exposure weights, including fractional phantom counts.
                trials <- TRUE
            }
        }
        if (is.null(audit.row)) {
            if (which.area == "All") {
                tmp <- subset(design, (time0 == which.time))
            } else {
                tmp <- subset(design, (time0 == which.time & region0 == as.character(which.area)))
            }
        }
        # We keep the existing empty/zero-death guards after augmentation. If the
        # reference also has no deaths, the corrected target can still return NA.
        if (dim(tmp)[1] == 0) {
            return(finish(rep(NA, 5)))
        } else if (sum(tmp$variables$died) == 0) {
            message(paste0(which.area, " ", which.time, " has no death, set to NA"))
            return(finish(rep(NA, 5)))
        } else if(length(unique(tmp$variables$age0)) > 1){
            if(!trials){
                glm.ob <- survey::svyglm(died ~ (-1) + factor(age0), design = tmp, family = stats::quasibinomial, maxit = 50)
            }else{
                glm.ob <- survey::svyglm(died ~ (-1) + factor(age0), design = tmp, family = stats::quasibinomial, maxit = 50, weights = tmp$variables$n)

            }
            if(length(stats::coef(glm.ob)) < length(labels)){
                bins.nodata <- length(labels) - length(stats::coef(glm.ob))
                if(bins.nodata >= length(ns)/2){
                    message(paste0(which.area, " ", which.time, " has no observation in more than half of the age bins, set to NA"))
                    return(finish(rep(NA, 5)))
                }
                # This can only happen for the last one or several bins, since person-month are cumulative
                ns.comb <- ns
                ns.comb[length(ns) - bins.nodata] <- sum(ns[c(length(ns) - bins.nodata) : length(ns)])
                message(paste0(which.area, " ", which.time, " has no observations in ", bins.nodata, " age bins. They are collapsed with previous bins"))
                return(finish(get.est.withNA(glm.ob, labels, ns.comb)))
            }
            return(finish(get.est.withNA(glm.ob, labels, ns)))
        } else {
              if(!trials){
                 glm.ob <- survey::svyglm(died ~ 1, design = tmp, family = stats::quasibinomial, maxit = 50)
              }else{
                  glm.ob <- survey::svyglm(died ~ 1, design = tmp, family = stats::quasibinomial, maxit = 50, weights = tmp$variables$n)
              }
            var.est <- stats::vcov(glm.ob)
            mean <- expit(stats::coef(glm.ob)[1])
            lims <- expit(logit(mean) + stats::qnorm(c((1 - CI) / 2, 1 - (1 - CI) / 2)) * sqrt(c(var.est)))

            ## Alternative calculation
            # p.i <- survey::svymean(~died, design = tmp)
            # var.i <- as.numeric(vcov(p.i)) 
            # p.i <- as.numeric(p.i)   
            # ht <- log(p.i/(1-p.i))
            # ht.v <- var.i/(p.i^2*(1-p.i)^2)
            # ht.prec <- 1/ht.v
           return(finish(c(mean, lims, logit(mean), var.est)))

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
        betas2 <- stats::coef(glm.ob)
        betas[names(betas2)] <- betas2
        
        # (removed) avoid NA by replacing with 0 This is wrong way to handle NA!  messed up the mean estimates
        # betas[which(is.na(betas))] <- 0
        
        # ns <- c(1, 11, 12, 12, 12, 12)
        probs <- expit(betas)
        
        mean.est <- (1 - prod((1 - probs)^ns, na.rm = TRUE))  #*1000  
        
        ## partial derivatives ##
        gamma <- prod((1 + exp(betas))^ns, na.rm = TRUE)
        derivatives <- (gamma)/(gamma - 1) * ns * expit(betas)
        
        ## now handle the NA in derivatives, it won't affect mean estimation
        derivatives[which(is.na(derivatives))] <- 0
        ## Items to return ##
        var.est <- t(derivatives) %*% V %*% derivatives
        lims <- logit(mean.est) + stats::qnorm(c((1 - CI) / 2, 1 - (1 - CI) / 2)) * sqrt(c(var.est))
        return(c(mean.est, expit(lims), logit(mean.est), var.est))
    }
    
    # run for all combinations
    x <- mapply(region.time.HT.withNA, which.area = results$region, which.time = results$years,
        target.index = seq_len(nrow(results)))
    results[, 4:8] <- t(x)
    results$survey <- NA
    results$logit.prec <- 1/results$var.est
    results <- results[, c("region", "years", "mean", "lower", "upper", "logit.est", "var.est", "survey", "logit.prec")]
     
    # We preserve the usual columns. fixed_areas lists all triggered regions;
    # variance.fix distinguishes successful and unresolved area-period attempts.
    if (isTRUE(var.fix)) {
        attr(results, "fixed_areas") <- unique(as.character(audit$region))
        attr(results, "variance.fix") <- audit
    }
    return(results)
}
