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
#' @param mc number of monte carlo draws to approximate the marginal prevalence/hazards for binomial model. If mc = 0, analytical approximation is used. The analytical approximation is invalid for hazard modeling with more than one age groups.
#' @param include_time_unstruct logical indicator whether to include the temporal unstructured effects (i.e., shocks) in the smoothed estimates.
#' @param ... additional configurations passed to inla.posterior.sample.
#' 
#' @return Results from RW2 model fit, including projection.
#' @seealso \code{\link{plot.SUMMERproj}}
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

#' @export
getSmoothed <- function(inla_mod, year_range = c(1985, 2019), year_label = c("85-89", "90-94", "95-99", "00-04", "05-09", "10-14", "15-19"), 
                            Amat = NULL, nsim = 1000, weight.strata = NULL, verbose = FALSE, mc = 0, include_time_unstruct = FALSE, ...){

      years <- NA

      ########################
      ## Binomial methods
      ########################
    if(!is.null(inla_mod$family)){
       if("region.struct" %in% names(inla_mod$fit$summary.random) == FALSE){
        warning("No spatial random effects in the model. Set Amat to NULL", immediate. = TRUE)
        Amat <- NULL
       }

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
            }else if("region" %in% colnames(weight.strata) && dim(weight.strata)[1] > 1){
              warning("Time is not specified, assuming the strata weights are static.", immediate.=TRUE)
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
               warning("no region or years in the strata proportion. Treat proportion as constant.", immediate.=TRUE)
               weight.strata.by = "Constant"
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


       message("Starting posterior sampling...")
        sampAll <- INLA::inla.posterior.sample(n = nsim, result = inla_mod$fit, intern = TRUE, selection = select, verbose = verbose, ...)
       message("Finished posterior sampling, cleaning up results now...")

        fields <- rownames(sampAll[[1]]$latent)        
        pred <- grep("Predictor", fields)
        time.unstruct <- grep("time.unstruct", fields)
        region.struct <- grep("region.struct", fields)
        region.unstruct <- grep("region.unstruct", fields)
        if(length(region.unstruct)==0) region.struct <- region.struct[1:(length(region.struct)/2)]
        region.int <- grep("region.int", fields)
        time.area <- grep("time.area", fields)
        age <- grep("age", fields)
        age.nn <- inla_mod$age.n
        strata <- grep("strata", fields)
        T <- length(time.unstruct) 
        if(!is.null(Amat)){
            N <- dim(Amat)[1]          
        }else{
            N <- 1
        }

        # Handle replicated RW
        time.struct <- grep("time.struct", fields)
        if(length(age) > 0){
          newindex <- rep(NA, T * length(age))
          for(i in 1:length(age)){
            where <- ((inla_mod$age.rw.group[i] - 1) * T + 1) : (inla_mod$age.rw.group[i] * T)
            newindex[((i-1)*T+1):(i*T)] <- where
          }          
          time.struct <- time.struct[newindex]
        }

        

        if(length(age) == 0){
          age.length  <- 1
        }else{
          age.length <- length(age)
        }
          

        # organize output
        out1 <- expand.grid(strata = stratalabels, time = 1:T, area = 1:N)
        out2 <- expand.grid(time = 1:T, area = 1:N)
        out1$lower <- out1$upper <- out1$mean <- out1$median <- out1$variance <- NA
        out2$lower <- out2$upper <- out2$mean <- out2$median <- out2$variance <- NA
        out1$years <- year_label[out1$time]
        out2$years <- year_label[out2$time]
        if(N > 1){              
          out1$region <- colnames(Amat)[out1$area]
          out2$region <- colnames(Amat)[out2$area]
        }else{
          out1$region <- out2$region <- "All"
        }

        if(is.null(weight.strata)){
          warning("No strata weights has been supplied. Equal weights are used to calculate the overall estimates. Please interpret results with caution or use only the stratified estimates.", immediate.=TRUE)
          for(tt in stratalabels){
            out2[, tt] <- 1/length(stratalabels)
          }

        }else if(weight.strata.by == "Constant"){
          for(tt in stratalabels){
            out2[, tt] <- weight.strata[1, match(tt, colnames(weight.strata))]
          }
          strata.index <- match(stratalabels, colnames(out2))
        }else{
          out2 <- merge(out2, weight.strata, by = weight.strata.by)
          strata.index <- match(stratalabels, colnames(out2))
          out2 <- out2[with(out2, order(area, time)), ]
        }
      

        ## T blocks, each with region 1 to N
        theta <- matrix(0, nsim, T * N)
        ## Age blocks, each with age 1 to T
        theta.rw <- matrix(NA, nsim, age.length * T)
        tau <- rep(NA, nsim)
        ## K blocks, each with age 1 to G
        beta <- matrix(NA, nsim, length(stratalabels) * age.length)

        time.area.order <- inla_mod$time.area[with(inla_mod$time.area, order( time.unstruct, region_number)), "time.area"]

        for(i in 1:nsim){
          if(length(time.area)> 0){
            theta[i, ] <- sampAll[[i]]$latent[time.area[time.area.order]]
          }else if(length(region.int) > 0){
            theta[i, ] <- sampAll[[i]]$latent[region.int]
          }else{
            theta[i, ] <- 0
          }
          if(include_time_unstruct) theta[i, ] <- theta[i, ] + rep(sampAll[[i]]$latent[time.unstruct], each = N)
          if(N > 1) theta[i, ] <- theta[i, ] + rep(sampAll[[i]]$latent[region.struct], T)
          theta.rw[i, ] <- sampAll[[i]]$latent[time.struct]

          if(length(region.unstruct)> 0){
            theta[i, ] <- theta[i, ] + rep(sampAll[[i]]$latent[region.unstruct], T)
          }

          if(length(age) > 0){
              beta[i, ] <- rep(sampAll[[i]]$latent[age], length(stratalabels))
          }else{
            beta[i, ] <- 0
          }

          if(length(strata) == length(stratalabels)){
              beta[i, ] <- beta[i, ] + rep(sampAll[[i]]$latent[strata], each = age.length)
          }else{
              beta[i, ] <- beta[i, ] + rep(c(0, sampAll[[i]]$latent[strata]), each = age.length)
          }
          if(inla_mod$family == "binomial"){
            tau[i] <-exp(sampAll[[i]]$hyperpar[["Log precision for nugget.id"]])
          }
        }


        ########################
        ## Beta-Binomial methods
        ########################        
        if(inla_mod$family == "binomial" && (mc == 0)){
          # again, column operation here
          k <- 16 * sqrt(3) / 15 / base::pi
          theta <- theta / sqrt(1 + k^2 / tau)
          theta.rw <- theta.rw / sqrt(1 + k^2 / tau)
          beta <- beta / sqrt(1 + k^2 / tau)
        }

        # Put hazards together
        draw.temp <- draws.hazards <- NA
        index1 <- 1
        index2 <- 1
        for(j in 1:N){
          # Monte Carlo approximation of the marginal effects that are shared across time
          if(inla_mod$family == "binomial" && mc > 0){
            sd.temp <- matrix(1/sqrt(tau), nsim, mc)
            err.temp <- matrix(stats::rnorm(nsim*mc, mean = matrix(0, nsim, mc), sd = sd.temp), nsim, mc)
          }

          for(i in 1:T){
            draws <- matrix(NA, nsim, length(stratalabels))
            # For each strata
            for(k in 1:(length(stratalabels))){
              # column add
              # This is a matix of dimension nsim * age.length
              age.rw <- theta.rw[, (1:(T*age.length)) %% T == ifelse(i == T, 0, i)]  
              draws.hazards <- theta[, j + (i-1)*N] + beta[, ((k-1)*age.length+1) :(k*age.length)] + age.rw            
              # Monte Carlo approximation of the marginal effects
              if(inla_mod$family == "binomial" && mc > 0){
                # sd.temp <- matrix(1/sqrt(tau), nsim, mc)
                # err.temp <- matrix(stats::rnorm(nsim*mc, mean = matrix(0, nsim, mc), sd = sd.temp), nsim, mc)
                for(tt in 1:age.length){
                    draws.temp <- matrix(draws.hazards[, tt], nsim, mc)
                    draws.temp <- expit(draws.temp + err.temp)
                    draws.hazards[, tt] <- apply(draws.temp, 1, mean)
                }
              }else{                 
                draws.hazards <- expit(draws.hazards)
              }


              if(!is.null(age.nn) && age.length > 1){
                draws.mort <- (1 - draws.hazards[, 1])^(age.nn[1])
                for(tt in 2:dim(draws.hazards)[2]){
                    draws.mort <- draws.mort * (1 - draws.hazards[, tt])^(age.nn[tt])
                }
                draws.mort <- 1 - draws.mort
              }else{
                draws.mort <- draws.hazards
              }
              draws[, k] <- draws.mort
              out1[index1, c("lower", "median", "upper")] <- quantile(draws.mort, c(0.025, 0.5, 0.975))
              out1[index1, "mean"] <- mean(draws.mort)
              out1[index1, "variance"] <- var(draws.mort)
              index1 <- index1 + 1
            }
            # aggregate to time-area estimates
            prop <- out2[index2, strata.index] 
            draws <- apply(draws, 1, function(x, p){sum(x * p)}, prop)
            out2[index2, c("lower", "median", "upper")] <- quantile(draws, c(0.025, 0.5, 0.975))
            out2[index2, "mean"] <- mean(draws)
            out2[index2, "variance"] <- var(draws)
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


