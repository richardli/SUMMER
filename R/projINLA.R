#' Extract smoothed estimates.
#' 
#' 
#'
#' @param inla_mod output from \code{\link{smoothDirect}} or \code{\link{smoothCluster}}
#' @param nsim number of simulations, only applicable for the cluster-level model. The smooth direct model always draws 1e5 samples from the marginal distribution since the computation is faster.
#' @param weight.strata a data frame with two columns specifying time and region, followed by columns specifying proportion of each strata for each region. This argument specifies the weights for strata-specific estimates on the probability scale.
#' @param weight.frame a data frame with three columns, years, region, and the weight of each frame for the corresponding time period and region. This argument specifies the weights for frame-specific estimates on the logit scale. Notice this is different from weight.strata argument. 
#' @param verbose logical indicator whether to print progress messages from inla.posterior.sample.
#' @param mc number of monte carlo draws to approximate the marginal prevalence/hazards for binomial model. If mc = 0, analytical approximation is used. The analytical approximation is invalid for hazard modeling with more than one age groups.
#' @param include_time_unstruct  Indicator whether to include the temporal unstructured effects (i.e., shocks) in the smoothed estimates from cluster-level model. The argument only applies to the cluster-level models (from \code{\link{smoothCluster}}). Default is FALSE which excludes all unstructured temporal components. If set to TRUE all  the unstructured temporal random effects will be included. Alternatively, if this is specified as a vector of   subset of year labels (as in the year_label argument), only the unstructured terms in the corresponding time periods will be added to the prediction.
#' @param include_subnational logical indicator whether to include the spatial and space-time interaction components in the smoothed estimates. If set to FALSE, only the main temporal trends are returned.
#' @param CI Desired level of credible intervals
#' @param draws Posterior samples drawn from the fitted model. This argument allows the previously sampled draws (by setting save.draws to be TRUE) be used in new aggregation tasks.  
#' @param save.draws Logical indicator whether the raw posterior draws will be saved. Saved draws can be used to accelerate aggregations with different weights.
#' @param ... Unused arguments, for users with fitted object from the package before v1.0.0, arguments including Amat, year_label, and year_range can still be specified manually.
#' 
#' @return A data frame or a list of data frames of S3 class SUMMERproj, which contains the smoothed estimates. 
#' 
#' @seealso \code{\link{plot.SUMMERproj}}
#' @author Zehang Richard Li
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
#' fit1 <- smoothDirect(data = data, Amat = NULL, 
#'   year_label = years.all, year_range = c(1985, 2019), 
#'   rw = 2, is.yearly=FALSE, m = 5)
#' out1 <- getSmoothed(fit1)
#' plot(out1, is.subnational=FALSE)
#' 
#' #  subnational model
#' fit2 <- smoothDirect(data = data, Amat = mat, 
#'   year_label = years.all, year_range = c(1985, 2019), 
#'   rw = 2, is.yearly=TRUE, m = 5, type.st = 4)
#' out2 <- getSmoothed(fit2)
#' plot(out2, is.yearly=TRUE, is.subnational=TRUE)
#' 
#' 
#' }
#' 

#' @export
getSmoothed <- function(inla_mod, nsim = 1000, weight.strata = NULL, weight.frame = NULL, verbose = FALSE, mc = 0, include_time_unstruct = FALSE, CI = 0.95, draws = NULL, save.draws = FALSE, include_subnational = TRUE, ...){

      years <- region <- age.diff <- NA
      lowerCI <- (1 - CI) / 2
      upperCI <- 1 - lowerCI
      save.draws.est <- save.draws
      msg <- NULL

      if(!is.null(inla_mod$year_range)){
        year_range <- inla_mod$year_range
      }else{
        warning("The fitted object was from an old version of SUMMER, please specify 'year_range' argument when calling getSmoothed()")
      }
      if(!is.null(inla_mod$year_label)){
        year_label <- inla_mod$year_label
      }else if(inla_mod$is.temporal){
        warning("The fitted object was from an old version of SUMMER, please specify 'year_label' argument when calling getSmoothed()")
      }
      if(!is.null(inla_mod$has.Amat)){
        Amat <- inla_mod$Amat
      }else{
        warning("The fitted object was from an old version of SUMMER, please specify 'Amat' argument when calling getSmoothed()")
      }
      if(!inla_mod$is.temporal){
        year_label <- NA
      }
      
      ########################
      ## Cluster level methods
      ########################
    if(!is.null(inla_mod$family)){
       if("region.struct" %in% names(inla_mod$fit$summary.random) == FALSE && !is.null(Amat)){
        warning("No spatial random effects in the model. Amat not used", immediate. = TRUE)
        Amat <- matrix(1,1,1)
        colnames(Amat) <- rownames(Amat) <- inla_mod$fit$.args$data$region[1]
       }
        is.dynamic <- as.logical(inla_mod$strata.time.effect)
        if(length(is.dynamic) == 0) is.dynamic = FALSE

        if(!is.dynamic){
          # check strata weights are properly specified: static strata effect case
          stratalabels <- stratalabels.orig <- inla_mod$strata.base                  
          if(length(stratalabels) == 1 && stratalabels[1] == ""){
            stratalabels <-  "strata_all"
            stratalabels.orig <- "strata_all"
          }
          other <- rownames(inla_mod$fit$summary.fixed)
          other <- other[grep("strata", other)]
          other <- gsub("strata", "", other)
          stratalabels <- c(stratalabels, other)
          stratalabels.orig <- c(stratalabels.orig, other)
          frame.strata <- (inla_mod$fit$.args$data)[, c("age", "age.idx", "age.rep.idx")]
          frame.strata <- unique(frame.strata)        
          frame.strata <- frame.strata[order(frame.strata$age.idx), ]
          framelabels <- "frame_all"
          multi.frame <- FALSE
        }else{
          # check strata weights are properly specified: dynamic strata effect case
          # here we have (age group) x (frame-strata) with no base level 
          if("frame" %in% colnames(inla_mod$fit$.args$data)){
            frame.strata <- (inla_mod$fit$.args$data)[, c("age", "age.idx", "age.rep.idx", "strata", "frame")]
            frame.strata <- unique(frame.strata)
            frame.strata <- frame.strata[!is.na(frame.strata$frame), ]            
            frame.strata$strata.orig <- NA
            for(i in 1:dim(frame.strata)[1]){
              frame.strata$strata.orig[i] <- gsub(paste0(frame.strata$frame[i], "-"), "", frame.strata$strata[i])
            }
          }else{
            frame.strata <- (inla_mod$fit$.args$data)[, c("age", "age.idx", "strata", "age.rep.idx")]
            frame.strata <- unique(frame.strata)
            frame.strata$strata.orig <- frame.strata$strata
          }
          frame.strata <- frame.strata[order(frame.strata$age.idx), ]
          stratalabels <- as.character(unique(frame.strata$strata))
          stratalabels.orig <- as.character(unique(frame.strata$strata.orig))
          framelabels <- as.character(unique(frame.strata$frame))
         
          if(length(framelabels)==0){
            framelabels <- "frame_all"
            multi.frame <- FALSE
          }else if(length(framelabels)==1){
            multi.frame <- FALSE
          }else if(sum(frame.strata$frame != frame.strata$strata) == 0){
            multi.frame <- FALSE
          }else{
            multi.frame <- TRUE
          }
        }
       
        if(inla_mod$is.yearly) year_label <- c(year_range[1]:year_range[2], year_label)
        if(!inla_mod$is.yearly) year_label <- year_label
        err <- NULL
        weight.strata.by <- NULL

        if(!is.null(weight.strata)){
          if((!"frame" %in% colnames(inla_mod$fit$.args$data)) && "frame" %in% colnames(weight.strata)){
            stop("frame variable is not in the fitted model, but exists in the weight.strata data frame. Please remove the column from weight.strata and use a single set of weights (that are not specific to sampling frames).")
          }

          if(!is.dynamic && sum(stratalabels %in% colnames(weight.strata)) != length(stratalabels)){
            stop(paste0("weight.strata argument not specified correctly. It requires the following columns: ", paste(stratalabels, collapse = ", ")))
          }
          if(is.dynamic && "frame" %in% colnames(frame.strata) && !("frame" %in% colnames(weight.strata))){
            stop("frame column was provided in the fitted object but not the weights.")
          }


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

       if(is.null(draws)){
          message("Starting posterior sampling...")
            sampAll <- INLA::inla.posterior.sample(n = nsim, result = inla_mod$fit, intern = TRUE, selection = select, verbose = verbose)
           message("Cleaning up results...")
     
       }else{
          message("Use posterior draws from input.")
          sampAll <- draws
          nsim <- length(draws)
       }

        fields <- rownames(sampAll[[1]]$latent)        
        T <- max(inla_mod$fit$.args$data$time.struct) 
        if(is.na(T)) T <- 1

        # Construct the AA matrix for the  output data frame.
        if(!"region.struct:1" %in% fields) inla_mod$fit$.args$data$region.struct = 1

        rep.time <- length(unique(inla_mod$age.rw.group)) > 0
        cols <- c("time.struct", "time.unstruct", "region.struct", "time.area", "strata", "age", "age.idx", "age.intercept", "age.diff")
        if(rep.time){
          cols <- c(cols, "age.rep.idx")
        } 
        if(!is.null(inla_mod$covariate.names)){
          cols <- c(cols, inla_mod$covariate.names)
        }
        A <- unique(inla_mod$fit$.args$data[, cols])
        if(sum(stratalabels == "strata_all") == length(stratalabels)){
          A$strata <- "strata_all"
        }
        # A might contain  missing  levels. Make an expanded version
        # strata may be nested within age
        if(inla_mod$strata.time.effect){
          AA  <-  expand.grid(time.area = unique(A$time.area), 
                              age = unique(A$age))
          AA <- merge(AA, unique(A[, c("age", "strata")]), by = "age")
        }else{
          AA  <-  expand.grid(time.area = unique(A$time.area), 
                            strata =  unique(A$strata),  
                            age = unique(A$age))
          if(sum(AA$strata != "") == 0) AA$strata <- "strata_all"
        }
        AA <- merge(AA, unique(A[, c("time.struct", "time.unstruct", "region.struct",  "time.area")]), by = "time.area")
        if("strata" %in% colnames(A)){
          if(rep.time){ 
            tmp <- unique(A[, c("strata", "age", "age.idx", "age.rep.idx", "age.intercept", "age.diff")])
          }else{            
            tmp <- unique(A[, c("strata", "age", "age.idx", "age.intercept", "age.diff")])
          }        
        }else{
          if(rep.time){ 
            tmp <- unique(A[, c("age", "age.idx", "age.rep.idx", "age.intercept", "age.diff")])
          }else{            
            tmp <- unique(A[, c("age", "age.idx", "age.intercept", "age.diff")])
          }    
        }
        # check the filler data contains any NA, which will make the merge call go wrong.
        if(sum(!is.na(tmp$age.diff)) > 0){
          tmp <- subset(tmp, !is.na(age.diff))
        }else{
          tmp <- tmp[, colnames(tmp) != "age.diff"]
        }
        if("strata" %in% colnames(tmp)){
          AA <- merge(AA, tmp, by = c("age", "strata"))
        }else{
          AA <- merge(AA, tmp, by = c("age"))
        }

        if(!is.null(inla_mod$covariate.names)){
          # adding covariates
          Asub <- data.frame(A[, inla_mod$covariate.names])
          if(sum(!is.na(A$region.struct)) > 0) Asub$region.struct <- A$region.struct
          if(sum(!is.na(A$time.unstruct)) > 0) Asub$time.unstruct <- A$time.unstruct
          Asub <- unique(Asub)
          AA <- merge(AA, Asub)
        }
        if(rep.time){
          AA$time.struct <- AA$time.struct + (AA$age.rep.idx - 1)  * T
        }
        if(!inla_mod$is.temporal){
          AA$time.struct <- AA$time.unstruct <- 1
        }        
        AA$age.intercept <- paste0("age.intercept", AA$age.intercept, ":1")
        AA$age.diff <- paste0("age.diff", AA$age.diff, ":1")

        # # When there's only one age group, smoothCluster uses the generic intercept
        # # if(length(unique(AA$age.intercept)) == 1) AA$age.intercept <- "(Intercept):1"
        
        # Checking age intercept is tricky, directly check if intercept is in the posterior draws...
        if("(Intercept):1" %in% fields){
          AA$intercept <- "(Intercept):1"
        }else{
          AA$intercept <- NA
        }

        #  AA.loc is the same  format as AA, but with location index
        AA.loc <- AA
        AA.loc$age.intercept  <- match(AA.loc$age.intercept, fields)
        AA.loc$age.diff  <- match(AA.loc$age.diff, fields)
        AA.loc$age <- NA
        # For constant case, base strata will end up being NA here, which is fine.
        if(!is.dynamic) AA.loc$strata <- paste0("strata", AA.loc$strata, ":1")
        AA.loc$strata  <- match(AA.loc$strata, fields)
        if(!is.null(inla_mod$covariate.names)){
            AA.loc[, inla_mod$covariate.names] <- NA
        }
        AA.loc$intercept <- match(AA.loc$intercept, fields)

        AA.loc$time.area  <- match(paste0("time.area:", AA.loc$time.area), fields)
        # Update time.area as the row index of the correct samples
        # when region.int and time.int is used
        if("region.int:1" %in% rownames(sampAll[[1]]$latent)){
          AA.loc$time.area <- (AA.loc$time.unstruct - 1) * dim(Amat)[1] + AA.loc$region.struct
          AA.loc$time.area  <- match(paste0("region.int:", AA.loc$time.area), fields)
        }
        AA.loc$time.struct  <- match(paste0("time.struct:", AA.loc$time.struct), fields)
        AA.loc$region.struct  <- match(paste0("region.struct:", AA.loc$region.struct), fields)
        AA.loc$time.unstruct  <- match(paste0("time.unstruct:", AA.loc$time.unstruct), fields)


        # if include_time_unstruct is a logical indicator
        if(is.logical(include_time_unstruct)){
          if(!include_time_unstruct) AA.loc$time.unstruct <- NA
        }else{
          # if it is a vector $TODO
          which.include <- which(year_label %in% as.character(include_time_unstruct))
          included <- which(AA$time.unstruct %in% which.include) 
          not_included <- which(AA$time.unstruct %in% which.include == FALSE ) 
          AA.loc$time.unstruct[not_included] <- NA
          text <- paste0("The IID temporal components are included in the following time periods: ", paste(year_label[which.include], collapse = ", "))
          message(text)
          msg <- paste0(msg, text)
        }
        slope <- grep("time.slope.group", fields)
        slope0 <- grep("time.slope:1", fields)
        if(length(slope) > 0){
           AA$tstar <- (AA$time.unstruct - (T + 1)/2) / (T + 1)
           AA$slope  <- match(paste0("time.slope.group", AA$age.rep.idx, ":1"), fields)
        }else if(length(slope0) > 0){
           AA$tstar <- (AA$time.unstruct - (T + 1)/2) / (T + 1)
           AA$slope  <- match(paste0("time.slope:1"), fields)
        }else{
          AA$tstar <- NA
          AA$slope <- "time.slope:NA"
        }
        st.slope <- grep("st.slope.id", fields)
        if(length(st.slope)>0){
           AA$ststar <- (AA$time.unstruct - (T + 1)/2) / (T + 1)
           AA$st.slope  <- match(paste0("st.slope.id:", AA$region.struct), fields)
        }else{
           AA$ststar <-  NA
           AA$st.slope <- "st.slope.id:NA"
        }
        AA.loc$age.idx <- AA.loc$age.rep.idx <- NA


        if(!include_subnational){
          AA.loc$time.area <- NA
          AA.loc$region.struct <- NA
        }

        # time.unstruct <- grep("time.unstruct", fields)
        # region.struct <- grep("region.struct", fields)
        # region.unstruct <- grep("region.unstruct", fields)
        # if(length(region.unstruct)==0) region.struct <- region.struct[1:(length(region.struct)/2)]
        # region.int <- grep("region.int", fields)
        # time.area <- grep("time.area", fields)
        age <- match(paste0("age", frame.strata$age, ":1"), fields)
        age.nn <- inla_mod$age.n#[frame.strata$age.idx]
        # if(!is.dynamic){
        #   strata <- grep("strata", fields)
        # }else{
        #   strata <- NULL
        # }
        # T <- length(time.unstruct) 
        if(!is.null(Amat)){
            N <- dim(Amat)[1]          
        }else{
            N <- 1
        }

        # Handle replicated RW
        # # The static case: time.struct is of length age * T, here age is not always 6, can be age-frame-strata comb 
        # time.struct <- grep("time.struct", fields)
        # if(length(age) > 0){
        #   newindex <- rep(NA, T * length(age))
        #   for(i in 1:length(age)){
        #     where <- ((inla_mod$age.rw.group[i] - 1) * T + 1) : (inla_mod$age.rw.group[i] * T)
        #     newindex[((i-1)*T+1):(i*T)] <- where
        #   }          
        #   time.struct <- time.struct[newindex]
        # }
      
        if(length(age) == 0){
          age.length  <- 1
        }else{
          age.length <- length(age)
        }
          

        # organize output
        out1 <- expand.grid(strata = stratalabels, time = 1:T, area = 1:N)
        out2 <- expand.grid(frame = framelabels, time = 1:T, area = 1:N)
        out3 <- expand.grid(time = 1:T, area = 1:N)
        out1$lower <- out1$upper <- out1$mean <- out1$median <- out1$variance <- NA
        out2$lower <- out2$upper <- out2$mean <- out2$median <- out2$variance <- NA
        out3$lower <- out3$upper <- out3$mean <- out3$median <- out3$variance <- NA
        out1$years <- year_label[out1$time]
        out2$years <- year_label[out2$time]
        out3$years <- year_label[out3$time]
        if(N > 1){              
          out1$region <- colnames(Amat)[out1$area]
          out2$region <- colnames(Amat)[out2$area]
          out3$region <- colnames(Amat)[out3$area]
        }else{
          out1$region <- out2$region <- out3$region <- "All"
        }

        if(is.null(weight.strata)){
          if(length(stratalabels) > 1){
              text <- "No strata weights has been supplied. Overall estimates are not calculated."
              message(text)
              msg <- paste0(msg, "\n", text)
          }else{
              text <- "No stratification in the model. Overall and stratified estimates are the same"
              message(text)
              msg <- paste0(msg,"\n",  text)
          }
          if(!is.null(Amat)){
            weight.strata <- expand.grid(region = colnames(Amat), frame = framelabels)
            weight.strata.by <- "region"
          }else{
            weight.strata <- data.frame(frame = framelabels)
            weight.strata.by <- "Constant"
          }
          
          for(tt in stratalabels.orig){
            weight.strata[, tt] <- ifelse(length(stratalabels) > 1, 0, 1)
          }
        }

        # if(weight.strata.by[1] == "Constant" && ){
        #   ## TODO: is this outdated?
        #   for(tt in stratalabels){
        #     out2[, tt] <- weight.strata[1, match(tt, colnames(weight.strata))]
        #   }
        #   strata.index <- match(stratalabels, colnames(out2))
        # } 

        if(weight.strata.by[1] == "Constant"){
          if("frame" %in% colnames(weight.strata)){
            out2 <- merge(out2, weight.strata, by = 'frame')
          }else{
            out2 <- cbind(out2, data.frame(weight.strata))
          }
          strata.index <- match(stratalabels.orig, colnames(out2))
        }else{          
          if(length(framelabels)  == 1){
            out2 <- merge(out2[, colnames(out2)!="frame"], weight.strata, by = weight.strata.by)
          }else{
            weight.strata.by <- c(weight.strata.by, "frame")
            out2 <- merge(out2, weight.strata, by = weight.strata.by)
          } 
          strata.index <- match(stratalabels.orig, colnames(out2))
          out2 <- out2[with(out2, order(area, time)), ]
        }
      

        # ## T blocks, each with region 1 to N
        # theta <- matrix(0, nsim, T * N)
        # ## Age blocks, each with time 1 to T
        # theta.rw <- matrix(NA, nsim, age.length * T)
        # tau <- rep(NA, nsim)
        # ## K blocks, each with age 1 to G
        # if(is.dynamic){
        #   # dynamic case age.length includes strata interactions already
        #   beta <- matrix(NA, nsim, age.length)
        # }else{
        #   beta <- matrix(NA, nsim, length(stratalabels) * age.length)
        # }

        # time.area.order <- inla_mod$time.area[with(inla_mod$time.area, order( time.unstruct, region_number)), "time.area"]

        # for(i in 1:nsim){
        #   if(length(time.area)> 0){
        #     theta[i, ] <- sampAll[[i]]$latent[time.area[time.area.order]]
        #   }else if(length(region.int) > 0){
        #     theta[i, ] <- sampAll[[i]]$latent[region.int]
        #   }else{
        #     theta[i, ] <- 0
        #   }
        #   if(include_time_unstruct) theta[i, ] <- theta[i, ] + rep(sampAll[[i]]$latent[time.unstruct], each = N)
        #   if(N > 1) theta[i, ] <- theta[i, ] + rep(sampAll[[i]]$latent[region.struct], T)
         
        #   theta.rw[i, ] <- sampAll[[i]]$latent[time.struct]
        #   if(!is.null(slope)){
        #      tstar <- rep(1:T, length(inla_mod$age.rw.group))
        #      tstar <- (tstar - T/2) / sd(1:T)
        #      theta.rw[i, ] <- theta.rw[i, ] + sampAll[[i]]$latent[slope][rep(inla_mod$age.rw.group, each=T)] * tstar
        #   }

        #   if(length(region.unstruct)> 0){
        #     theta[i, ] <- theta[i, ] + rep(sampAll[[i]]$latent[region.unstruct], T)
        #   }

        #   if(length(age) > 0 && !is.dynamic){
        #       beta[i, ] <- rep(sampAll[[i]]$latent[age], length(stratalabels))
        #   }else if(length(age) > 0){
        #       beta[i, ] <- sampAll[[i]]$latent[age]
        #   }else{
        #     beta[i, ] <- 0
        #   }

      
        #   if(length(strata) == length(stratalabels)){
        #       beta[i, ] <- beta[i, ] + rep(sampAll[[i]]$latent[strata], each = age.length)
        #   }else{
        #       beta[i, ] <- beta[i, ] + rep(c(0, sampAll[[i]]$latent[strata]), each = age.length)
        #   }
          # if(inla_mod$family == "binomial"){
          #   tau[i] <-exp(sampAll[[i]]$hyperpar[["Log precision for nugget.id"]])
          # }
        # }
        tau <- rep(NA, nsim)
        theta <- matrix(0, nsim, dim(AA)[1])
        for(i in 1:nsim){
          draw <- sampAll[[i]]$latent
          theta[i,  ] <-  apply(AA.loc, 1, function(x, ff){sum(ff[x], na.rm=TRUE)}, draw)
          add.slope <- draw[AA$slope] * AA$tstar
          add.slope[is.na(add.slope)] <- 0
          add.slope.st <- draw[AA$st.slope] * AA$ststar
          add.slope.st[is.na(add.slope.st)] <- 0
          theta[i,  ] <-  theta[i,  ] + add.slope + add.slope.st
          
          if(!is.null(inla_mod$covariate.names)){
            for(xx in inla_mod$covariate.names){
                covariate <- AA[, xx]
                slope <- draw[paste0(xx, ":1"), ]
                add.cov.effect <- covariate * slope
                add.cov.effect[is.na(add.cov.effect)] <- 0
                theta[i,  ] <-  theta[i,  ] + add.cov.effect
            }
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
          # theta.rw <- theta.rw / sqrt(1 + k^2 / tau)
          # beta <- beta / sqrt(1 + k^2 / tau)
        }

        # Put hazards together
        draw.temp <- draws.hazards <- NA
        draws.est <- NULL
        index.draws.est <- 1
        draws.est.overall <- NULL
        index.draws.est.overall <- 1
        index1 <- 1
        index2 <- 1
        if(N == 1 && AA$region.struct[1] == 0) AA$region.struct <- 1
        for(j in 1:N){
          # Monte Carlo approximation of the marginal effects that are shared across time
          if(inla_mod$family == "binomial" && mc > 0){
            sd.temp <- matrix(1/sqrt(tau), nsim, mc)
            err.temp <- matrix(stats::rnorm(nsim*mc, mean = matrix(0, nsim, mc), sd = sd.temp), nsim, mc)
          }

          for(i in 1:T){
            sub <-  which(AA$time.unstruct==i & AA$region.struct == j)
            AA.sub <- AA[sub, ] 
            draws.sub <- theta[, sub, drop = FALSE]
            draws.sub.agg <- matrix(NA, nsim, length(stratalabels))
            # For each strata
            for(k in 1:length(stratalabels)){
              strata.sub <- which(AA.sub$strata == stratalabels[k])
              draws.hazards <- draws.sub[, strata.sub, drop=FALSE]
              # Monte Carlo approximation of the marginal effects
              if(inla_mod$family == "binomial" && mc > 0){
                for(tt in 1:dim(draws.hazards)[2]){
                    draws.temp <- matrix(draws.hazards[, tt], nsim, mc)
                    draws.temp <- expit(draws.temp + err.temp)
                    draws.hazards[, tt] <- apply(draws.temp, 1, mean)
                }
              }else{                 
                draws.hazards <- expit(draws.hazards)
              }
              draws.mort <- rep(1, dim(draws.hazards)[1])
              for(tt in 1:dim(draws.hazards)[2]){
                draws.mort <- draws.mort * (1 - draws.hazards[, tt])^age.nn[AA.sub[strata.sub, "age.idx"][tt]]                
              }
              draws.mort <- 1 - draws.mort
              draws.sub.agg[, k] <- draws.mort
              ## Strata specific output
              index1 <- which(out1$time == i & out1$area == j & out1$strata == stratalabels[k])
              
              # save stratified mortality draws
              if(save.draws.est){
                  draws.est[[index.draws.est]] <- list(years = year_label[i], region = colnames(Amat)[j], strata = stratalabels[k], draws = draws.mort)
                  index.draws.est <- index.draws.est + 1
              }
              out1[index1, c("lower", "median", "upper")] <- stats::quantile(draws.mort, c(lowerCI, 0.5, upperCI))
              out1[index1, "mean"] <- mean(draws.mort)
              out1[index1, "variance"] <- var(draws.mort)
            }

            # aggregate across strata
            index2 <- which(out2$area == j & out2$time == i)
            draws.sub.agg.sum <- matrix(NA, nsim, length(index2))
            for(k in 1:length(index2)){
              prop <- out2[index2[k], strata.index]
              if(!multi.frame){
                  cols <- match(stratalabels.orig, stratalabels)
              }else{
                  cols <- match(paste(out2$frame[index2[k]], stratalabels.orig, sep = "-"), stratalabels)
              }

              draws.sub.agg.sum[, k] <- apply(draws.sub.agg[, cols, drop=FALSE], 1, function(x, p){sum(x * p)}, prop)

              # save overall mortality draws
              if(save.draws.est){
                  draws.est.overall[[index.draws.est.overall]] <- list(years = year_label[i], region = colnames(Amat)[j], draws = draws.sub.agg.sum[, k])
                  index.draws.est.overall <- index.draws.est.overall + 1
              }

              out2[index2[k], c("lower", "median", "upper")] <- stats::quantile(draws.sub.agg.sum[, k], c(lowerCI, 0.5, upperCI))
              out2[index2[k], "mean"] <- mean(draws.sub.agg.sum[, k])
              out2[index2[k], "variance"] <- var(draws.sub.agg.sum[, k])
            }

            # aggregate across frame
            if(!is.null(weight.frame)){
                index3 <- which(out3$area == j & out3$time == i)
                colnames(draws.sub.agg.sum) <- out2[index2, "frame"]
                draws.sub.agg.sum2 <- rep(0, nsim)
                this.weight <- weight.frame
                if("region" %in% colnames(weight.frame)) this.weight <- subset(this.weight, region == colnames(Amat)[j])
                if("years" %in% colnames(weight.frame)) this.weight <- subset(this.weight, years == year_label[i])

                for(k in 1:dim(draws.sub.agg.sum)[2]){
                  # aggregation on the probability scale, no longer used
                  # draws.sub.agg.sum2 <- draws.sub.agg.sum2 + draws.sub.agg.sum[, k] * as.numeric(this.weight[colnames(draws.sub.agg.sum)[k]])
                  draws.sub.agg.sum2 <- draws.sub.agg.sum2 + logit(draws.sub.agg.sum[, k]) * as.numeric(this.weight[colnames(draws.sub.agg.sum)[k]])
                }
                draws.sub.agg.sum2 <- expit(draws.sub.agg.sum2)
                out3[index3, c("lower", "median", "upper")] <- stats::quantile(draws.sub.agg.sum2, c(lowerCI, 0.5, upperCI))
                out3[index3, "mean"] <- mean(draws.sub.agg.sum2)
                out3[index3, "variance"] <- var(draws.sub.agg.sum2)
            }

            # Area-time specific done
          }
        }
      
      out1$is.yearly <- !(out1$years %in% year_label)
      out1$years.num <- suppressWarnings(as.numeric(as.character(out1$years)))
      out2$is.yearly <- !(out2$years %in% year_label)
      out2$years.num <- suppressWarnings(as.numeric(as.character(out2$years)))
      out3$is.yearly <- !(out3$years %in% year_label)
      out3$years.num <- suppressWarnings(as.numeric(as.character(out3$years)))
      out1$years <- factor(out1$years, year_label)
      out2$years <- factor(out2$years, year_label)
      out3$years <- factor(out3$years, year_label)
      class(out1) <- c("SUMMERproj", "data.frame")
      class(out2) <- c("SUMMERproj", "data.frame")
      class(out3) <- c("SUMMERproj", "data.frame")

      out <- list(overall = out2, stratified = out1)
      if(!is.null(weight.frame)){
        out$final = out3
      }
      if(save.draws){
        out$draws = sampAll
        msg <- paste0(msg, "\nPosterior draws are saved in the output. You can use 'getSmoothed(..., draws = ...$draws)' next time to speed up the call.")
      }
      if(save.draws.est){
        out$draws.est <- draws.est
        out$draws.est.overall <- draws.est.overall
      }
      out$msg <- msg
      out$nsim <- nsim
      class(out) <- "SUMMERprojlist"
      return(out) 


      ########################
      ## Mercer et al. methods
      ########################
    }else{

      if (!isTRUE(requireNamespace("INLA", quietly = TRUE))) {
        stop("You need to install the packages 'INLA'. Please run in your R terminal:\n  install.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/stable'), dep=TRUE)")
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
      if(!inla_mod$is.temporal) timelabel.yearly <-1
      results <- expand.grid(District = region_nums, Year = timelabel.yearly)
      results$median <- results$lower <- results$upper <- results$logit.median <- results$logit.lower <- results$logit.upper <- NA
      mod <- inla_mod$fit
      lincombs.info <- inla_mod$lincombs.info

      for(i in 1:length(timelabel.yearly)){
        for(j in 1:length(region_names)){
            index <- lincombs.info$Index[lincombs.info$District == region_nums[j] & lincombs.info$Year == i]
            tmp.logit <- INLA::inla.rmarginal(1e5, mod$marginals.lincomb.derived[[index]])
            # marg <- INLA::inla.tmarginal(expit, mod$marginals.lincomb.derived[[index]])
            # tmp <- INLA::inla.rmarginal(nsim, marg)
            tmp <- expit(tmp.logit)

            results$median[results$District == region_nums[j] & results$Year == timelabel.yearly[i]] <- stats::median(tmp)
            results$upper[results$District == region_nums[j] & results$Year == timelabel.yearly[i]] <- stats::quantile(tmp, upperCI)
            results$lower[results$District == region_nums[j] & results$Year == timelabel.yearly[i]] <- stats::quantile(tmp, lowerCI)
            results$logit.median[results$District == region_nums[j] & results$Year == timelabel.yearly[i]] <- stats::median(tmp.logit)
            results$logit.upper[results$District == region_nums[j] & results$Year == timelabel.yearly[i]] <- stats::quantile(tmp.logit, upperCI)
            results$logit.lower[results$District == region_nums[j] & results$Year == timelabel.yearly[i]] <- stats::quantile(tmp.logit, lowerCI)

        }
      }

      ## TODO: In the future, extract posterior draws with inla.posterior.sample and fit1$lincombs.info and fit1$fit$.args$lincomb 

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
      if(!inla_mod$is.temporal) results$years <- results$years.num <- results$is.yearly <- NA
      # Add S3 method
      class(results) <- c("SUMMERproj", "data.frame")
      return(results)
      }
  }


}


