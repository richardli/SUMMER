#' Function to obtain projected estimates from INLA for each time and region.
#' 
#' 
#'
#' @param inla_mod output from \code{\link{fitINLA}}
#' @param year_range range corresponding to year label
#' @param year_label vector of year string vector
#' @param Amat adjacency matrix
#' @param nsim number of simulations
#' @param weight.strata a data frame with three columns, years, region, and proportion of each strata for the corresponding time period and region. 
#' @param verbose logical indicator whether to print progress messages from inla.posterior.sample.
#' @param ... additional configurations passed to inla.posterior.sample.
#' 
#' @return Results from RW2 model fit, including projection.
#' 
#' @examples
#' \dontrun{
#' years <- levels(DemoData[[1]]$time)
#' 
#' # obtain direct estimates
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
#' fit1 <- fitINLA(data = data, geo = NULL, Amat = NULL, 
#'   year_label = years.all, year_range = c(1985, 2019), 
#'   rw = 2, is.yearly=FALSE, m = 5)
#' out1 <- getSmoothed(fit1)
#' plot(out1, is.subnational=FALSE)
#' 
#' #  subnational model
#' fit2 <- fitINLA(data = data, geo = geo, Amat = mat, 
#'   year_label = years.all, year_range = c(1985, 2019), 
#'   rw = 2, is.yearly=TRUE, m = 5, type.st = 4)
#' out2 <- getSmoothed(fit2, Amat = mat)
#' plot(out2, is.yearly=TRUE, is.subnational=TRUE)
#' 
#' 
#' }
#' 
#' 
#' @export
getSmoothed <- function(inla_mod, year_range = c(1985, 2019), year_label = c("85-89", "90-94", "95-99", "00-04", "05-09", "10-14", "15-19"), 
                            Amat = NULL, nsim = 1000, weight.strata = NULL, verbose = FALSE, ...){

      years <- NA

      ########################
      ## Binomial methods
      ########################
    if(!is.null(inla_mod$family)){


        # check strata weights are properly specified
        stratalabels <- inla_mod$strata.base
        other <- rownames(inla_mod$fit$summary.fixed)
        other <- other[grep("strata", other)]
        other <- gsub("strata", "", other)
        stratalabels <- c(stratalabels, other)
        if(inla_mod$is.yearly) year_label <- c(year_range[1]:year_range[2], year_label)
        if(!inla_mod$is.yearly) year_label <- year_label

        err <- NULL
        weight.strata.by <- NULL

        if(!is.null(weight.strata)){

          if(sum(stratalabels %in% colnames(weight.strata)) != length(stratalabels)) stop(paste0("weight.strata argument not specified correctly. It requires the following columns: ", paste(stratalabels, collapse = ", ")))

           if(sum(c("region", "years") %in% colnames(weight.strata)) == 2){
              for(i in year_label){
                tmp <- colnames(Amat)[which(colnames(Amat) %in% subset(weight.strata, years == i)$region == FALSE)]
                if(length(tmp) > 0) err <- c(err, paste(tmp, i))
              }
              if(!is.null(err)){
                stop(paste0("The following region-year combinations are not present in the strata weights: ", paste(err, collapse = ", ")))
              }
              weight.strata.by <- c("region", "years")
            }else if("region" %in% colnames(weight.strata)){
              warning("Time is not specified, assuming the strata weights are static.")
                tmp <- colnames(Amat)[which(colnames(Amat) %in% weight.strata$region == FALSE)]
                if(length(tmp) > 0){
                    stop(paste0("The following regions are not present in the strata weights: ", paste(tmp, collapse = ", ")))
                 }
              weight.strata.by <- c("region")
            }else if("years" %in% colnames(weight.strata)){
                tmp <- year_label[which(year_label %in% weight.strata$years == FALSE)]
                if(length(tmp) > 0){
                    stop(paste0("The following time periods are not present in the strata weights: ", paste(tmp, collapse = ", ")))
                 }
              weight.strata.by <- c("years")
            }else{
               stop("weight.strata argument not specified correctly. It requires one of the following columns: region, years.")
            }  
        }

       # get a subset of fields only
       cs <- inla_mod$fit$misc$configs$contents$tag
       cs <- cs[cs != "Predictor"]
       cs <- cs[cs != "nugget.id"]
       select <- list()
       for(i in 1:length(cs)){
          select[[i]] <- 0
          names(select)[i] <- cs[i]
       }

        sampAll <- INLA::inla.posterior.sample(n = nsim, result = inla_mod$fit, intern = TRUE, selection = select, verbose = verbose, ...)

        fields <- rownames(sampAll[[1]]$latent)        
        pred <- grep("Predictor", fields)
        time.struct <- grep("time.struct", fields)
        time.unstruct <- grep("time.unstruct", fields)
        region.struct <- grep("region.struct", fields)
        region.unstruct <- grep("region.unstruct", fields)
        if(length(region.unstruct)==0) region.struct <- region.struct[1:(length(region.struct)/2)]
        region.int <- grep("region.int", fields)
        time.area <- grep("time.area", fields)
        age <- grep("age", fields)
        strata <- grep("strata", fields)
        T <- length(time.struct)
        N <- dim(Amat)[1]
          

        # organize output
        out1 <- expand.grid(strata = stratalabels, time = 1:T, area = 1:N)
        out2 <- expand.grid(time = 1:T, area = 1:N)
        out1$lower <- out1$upper <- out1$mean <- out1$median <- NA
        out2$lower <- out2$upper <- out2$mean <- out2$median <- NA
        out1$years <- year_label[out1$time]
        out2$years <- year_label[out2$time]
        out1$region <- colnames(Amat)[out1$area]
        out2$region <- colnames(Amat)[out2$area]

        if(is.null(weight.strata)){
          warning("No strata weights has been supplied. Equal weights are used to calculate the overall estimates. Please interpret results with caution or use only the stratified estimates.")
          for(tt in stratalabels){
            out2[, tt] <- 1/length(stratalabels)
          }

        }else{
          out2 <- merge(out2, weight.strata, by = weight.strata.by)
          strata.index <- match(stratalabels, colnames(out2))
          out2 <- out2[with(out2, order(area, time)), ]
        }
      

        ## T blocks, each with region 1 to N
        theta <- matrix(NA, nsim, T * N)
        tau <- rep(NA, nsim)
        ## K blocks, each with age 1 to G
        beta <- matrix(NA, nsim, (length(strata)+1) * length(age))

        time.area.order <- inla_mod$time.area[with(inla_mod$time.area, order( time.unstruct, region_number)), "time.area"]

        for(i in 1:nsim){
          if(length(time.area)> 0){
            # TODO: double check this ordering
            theta[i, time.area.order] <- sampAll[[i]]$latent[time.area]
          }else{
            theta[i, ] <- sampAll[[i]]$latent[region.int]
          }
          theta[i, ] <- theta[i, ] + rep(sampAll[[i]]$latent[time.struct], each = N)
          theta[i, ] <- theta[i, ] + rep(sampAll[[i]]$latent[time.unstruct], each = N)
          theta[i, ] <- theta[i, ] + rep(sampAll[[i]]$latent[region.struct], T)
          if(length(region.unstruct)> 0){
            theta[i, ] <- theta[i, ] + rep(sampAll[[i]]$latent[region.unstruct], T)
          }

          beta[i, ] <- rep(sampAll[[i]]$latent[age], length(strata) + 1)
          beta[i, ] <- beta[i, ] + rep(c(0, sampAll[[i]]$latent[strata]), each = length(age))
          if(inla_mod$family == "binomial"){
            tau[i] <-exp(sampAll[[i]]$hyperpar[["Log precision for nugget.id"]])
          }
        }


        ########################
        ## Beta-Binomial methods
        ########################        
        if(inla_mod$family == "betabinomial"){
        # Does nothing, this is already marginal
        
        ########################
        ## Logistic-Normal-Binomial methods
        ########################    
        }else{
          # again, column operation here
          k <- 16 * sqrt(3) / 15 / base::pi
          theta <- theta / sqrt(1 + k^2 / tau)
        }

        # Put hazards together
        index1 <- 1
        index2 <- 1
        for(j in 1:N){
          for(i in 1:T){
            draws <- matrix(NA, nsim, length(strata) + 1)
            # For each strata
            for(k in 1:(length(strata)+1)){
              # column add
              draws.hazards <- expit(theta[, j + (i-1)*N] + beta[, ((k-1)*length(age)+1) :(k*length(age))])
              draws.mort <- (1 - draws.hazards[, 1])^(inla_mod$age.n[1])
              for(tt in 2:dim(draws.hazards)[2]){
                  draws.mort <- draws.mort * (1 - draws.hazards[, tt])^(inla_mod$age.n[tt])
              }
              draws.mort <- 1 - draws.mort
              draws[, k] <- draws.mort
              out1[index1, c("lower", "median", "upper")] <- quantile(draws.mort, c(0.025, 0.5, 0.975))
              out1[index1, "mean"] <- mean(draws.mort)
              index1 <- index1 + 1
            }
            # aggregate to time-area estimates
            prop <- out2[index2, strata.index] 
            draws <- apply(draws, 1, function(x, p){sum(x * p)}, prop)
            out2[index2, c("lower", "median", "upper")] <- quantile(draws, c(0.025, 0.5, 0.975))
            out2[index2, "mean"] <- mean(draws)
            index2 <- index2 + 1
          }
        }

      out1$is.yearly <- !(out1$years %in% year_label)
      out1$years.num <- suppressWarnings(as.numeric(as.character(out1$years)))
      out2$is.yearly <- !(out2$years %in% year_label)
      out2$years.num <- suppressWarnings(as.numeric(as.character(out2$years)))
   
      # ## deal with 1 year case
      # tmp <- as.numeric(strsplit(out1$years[out1$time == 1][1], "-")[[1]])
      # if(tmp[1] == tmp[2]){
      #     out1$years.num <- 1900 + tmp[1] + out1$time - 1
      #     out2$years.num <- 1900 + tmp[1] + out2$time - 1
      # }
      out1$years <- factor(out1$years, year_label)
      out2$years <- factor(out2$years, year_label)
      class(out1) <- c("SUMMERproj", "data.frame")
      class(out2) <- c("SUMMERproj", "data.frame")
      return(list(overall = out2, stratified = out1)) 


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
      results <- expand.grid(District = region_nums, Year = timelabel.yearly)
      results$median <- results$lower <- results$upper <- results$logit.median <- results$logit.lower <- results$logit.upper <- NA
      mod <- inla_mod$fit
      lincombs.info <- inla_mod$lincombs.info

      for(i in 1:length(timelabel.yearly)){
        for(j in 1:length(region_names)){
            index <- lincombs.info$Index[lincombs.info$District == region_nums[j] & lincombs.info$Year == i]
            tmp.logit <- INLA::inla.rmarginal(nsim, mod$marginals.lincomb.derived[[index]])
            marg <- INLA::inla.tmarginal(expit, mod$marginals.lincomb.derived[[index]])
            tmp <- INLA::inla.rmarginal(nsim, marg)

            results$median[results$District == region_nums[j] & results$Year == timelabel.yearly[i]] <- stats::median(tmp)
            results$upper[results$District == region_nums[j] & results$Year == timelabel.yearly[i]] <- stats::quantile(tmp, .975)
            results$lower[results$District == region_nums[j] & results$Year == timelabel.yearly[i]] <- stats::quantile(tmp, .025)
            results$logit.median[results$District == region_nums[j] & results$Year == timelabel.yearly[i]] <- stats::median(tmp.logit)
            results$logit.upper[results$District == region_nums[j] & results$Year == timelabel.yearly[i]] <- stats::quantile(tmp.logit, .975)
            results$logit.lower[results$District == region_nums[j] & results$Year == timelabel.yearly[i]] <- stats::quantile(tmp.logit, .025)

        }
      }
      results$is.yearly <- !(results$Year %in% year_label)
      results$years.num <- suppressWarnings(as.numeric(as.character(results$Year)))
      #  ## deal with 1 year case
      # results$years <- as.character(results$years)
      # tmp <- as.numeric(strsplit(results$years[results$time == 1][1], "-")[[1]])
      # if(tmp[1] == tmp[2]){
      #     results$years.num <- 1900 + tmp[1] + results$time - 1
      # }
      if(region_names[1] != "All"){
        results$District <- region_names[results$District]
      }else{
        results$District <- "All"
      }
      colnames(results)[which(colnames(results) == "District")] <- "region"
      colnames(results)[which(colnames(results) == "Year")] <- "years"
      # Add S3 method
      class(results) <- c("SUMMERproj", "data.frame")
      return(results)
      }
  }


}


