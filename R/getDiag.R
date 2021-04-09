#' Extract posterior summaries of random effects
#' 
#' @param inla_mod output from \code{\link{smoothDirect}}
#' @param field which random effects to plot. It can be one of the following: space, time, and spacetime.
#' @param CI Desired level of credible intervals
#' @param draws Posterior samples drawn from the fitted model. This argument allows the previously sampled draws (by setting save.draws to be TRUE) be used in new aggregation tasks.  
#' @param ... Unused arguments, for users with fitted object from the package before v1.0.0, arguments including Amat, year_label, and year_range can still be specified manually.
#' 
#' @return List of diagnostic plots
#' @author Zehang Richard Li
#' @examples
#' \dontrun{
#'   data(DemoMap)
#'   years <- levels(DemoData[[1]]$time)
#'   
#'   # obtain direct estimates
#'   data <- getDirectList(births = DemoData, 
#'   years = years,  
#'   regionVar = "region", timeVar = "time", 
#'   clusterVar = "~clustid+id", 
#'   ageVar = "age", weightsVar = "weights", 
#'   geo.recode = NULL)
#'   # obtain direct estimates
#'   data_multi <- getDirectList(births = DemoData, years = years,
#'     regionVar = "region",  timeVar = "time", clusterVar = "~clustid+id",
#'     ageVar = "age", weightsVar = "weights", geo.recode = NULL)
#'   data <- aggregateSurvey(data_multi)
#'   
#'   #  national model
#'   years.all <- c(years, "15-19")
#'   fit1 <- smoothDirect(data = data, geo = DemoMap$geo, Amat = DemoMap$Amat, 
#'     year_label = years.all, year_range = c(1985, 2019), 
#'     rw = 2, is.yearly=FALSE, m = 5)
#' random.time <- getDiag(fit1, field = "time")
#'   random.space <- getDiag(fit1, field = "space")
#'   random.spacetime <- getDiag(fit1, field = "spacetime")
#' }
#' 
#' @export


getDiag <- function(inla_mod, field = c("space", "time", "spacetime")[1], CI = 0.95, draws = NULL, ...){
	lower <- (1 - CI) / 2
	upper <- 1 - lower
	if(!is.null(inla_mod$year_range)){
		year_range <- inla_mod$year_range
	}else{
		warning("The fitted object was from an old version of SUMMER, please specify 'year_range' argument when calling getDiag()")
	}
	if(!is.null(inla_mod$year_label)){
		year_label <- inla_mod$year_label
	}else{
		warning("The fitted object was from an old version of SUMMER, please specify 'year_label' argument when calling getDiag()")
	}
	if(!is.null(inla_mod$has.Amat)){
        Amat <- inla_mod$Amat
      }else{
        warning("The fitted object was from an old version of SUMMER, please specify 'Amat' argument when calling getDiag()")
      }
	getquants <- function(mlist, lower, upper){
		quants <- data.frame(matrix(NA, length(mlist), 3))
		for(i in 1:length(mlist)){
			quants[i, ] <- INLA::inla.qmarginal(c(lower, 0.5, upper), mlist[[i]])
		}
		colnames(quants) <- c("lower", "median", "upper")
		return(quants)
	}

	is.yearly = inla_mod$is.yearly
	has.slope <- sum(grepl("slope", rownames(inla_mod$fit$summary.fixed)))

	if(field == "time" && has.slope == 0){
		struct <- inla_mod$fit$marginals.random$time.struct
		unstruct <- inla_mod$fit$marginals.random$time.unstruct
		if(is.yearly){
		    label <- label.unstruct <- c(year_range[1] : year_range[2], year_label)
		  }else{
		    label <- label.unstruct <- year_label
		 }
		if(!is.null(inla_mod$age.rw.group)){
			group <- rep(inla_mod$age.groups, each = length(label))
			label <- rep(label, length(inla_mod$age.rw.group))
		} 
		expand <- 1
		if(!is.null(inla_mod$age.rw.group)) expand <-  max(inla_mod$age.rw.group) / length(inla_mod$age.rw.group) 
		if(length(struct) != length(label) * expand) stop("The input year_range or year_label does not match the fitted model. Please double check.")
		temp <- getquants(struct, lower = lower, upper = upper)	
		quants <- NULL
		if(!is.null(inla_mod$age.rw.group)){
			for(i in 1:length(inla_mod$age.groups)){
				where <- (inla_mod$age.rw.group[i] - 1) * length(year_label) + c(1:length(year_label))
				quants <- rbind(quants,  temp[where, ])
			}
		}else{
			quants <- getquants(struct, lower = lower, upper = upper)
		}
		n <- dim(quants)[1]
		m <- length(unstruct)
		quants <- rbind(quants, getquants(unstruct, lower = lower, upper = upper))
		quants$years <- c(label, label.unstruct)
		if(!is.null(inla_mod$age.rw.group)){
			quants$group <- c(group, rep(NA, m))
		}
	  	quants$years.num <- suppressWarnings(as.numeric(as.character(quants$years)))
		quants$label <- c(rep("RW", n), rep("IID", m))
		quants$is.yearly <- !(quants$years %in% year_label)

	}else if(field == "time" && has.slope > 0){

		fixed <- rownames(inla_mod$fit$summary.fixed)
		fixed <- fixed[grepl("slope", fixed)]
		fixed <- c(fixed, "time.struct")
        select <- list()
        for(i in 1:length(fixed)){
           select[[i]] <- 0
           names(select)[i] <- fixed[i]
        }
        if(is.null(draws)){
	        message("AR1 model diagnostics are still experimental, starting posterior sampling...")
	        sampAll <- INLA::inla.posterior.sample(n = 1e3, result = inla_mod$fit, intern = TRUE, selection = select, verbose = FALSE)        	
        }else{
        	sampAll <- draws
        }

	    #inla_mod$fit$marginals.random$time.struct 
	    re <- grep("time.struct", rownames(sampAll[[1]]$latent))
	    fe <- grep("time.slope.group", rownames(sampAll[[1]]$latent))
	    fe0 <- grep("time.slope:1", rownames(sampAll[[1]]$latent))
	    if(length(fe0) > 0){
	    	fe <- rep(fe0, length(inla_mod$age.rw.group))
	    } 

	    struct.all <- matrix(0, length(re), length(sampAll))
	    T <- length(re) / length(fe)
	    xx <- ((1:T) -  (T + 1)/2) / (T + 1)
	    for(j in 1:length(sampAll)){
	    	for(k in 1:length(inla_mod$age.rw.group)){
	    		group.index <- inla_mod$age.rw.group[k]
	    		where <- ((group.index-1) * T + 1) :(group.index * T )
		    	struct.all[where, j] <- sampAll[[j]]$latent[re[where], 1] + sampAll[[j]]$latent[fe[group.index], 1]	* xx    		
	    	}
	    }
	    temp <- data.frame(t(apply(struct.all, 1, stats::quantile, c(lower, 0.5, upper))))
	    colnames(temp) <- c("lower", "median", "upper")
	    rownames(temp) <- NULL
		unstruct <- inla_mod$fit$marginals.random$time.unstruct
		if(is.yearly){
		    label <- label.unstruct <- c(year_range[1] : year_range[2], year_label)
		  }else{
		    label <- label.unstruct <- year_label
		 }
		if(!is.null(inla_mod$age.rw.group)){
			group <- rep(inla_mod$age.groups, each = length(label))
			# Now in this version, do not repeat the same effects
			label <- rep(label, length(inla_mod$age.rw.group))
		} 
		expand <- 1

		if(!is.null(inla_mod$age.rw.group)) expand <-  max(inla_mod$age.rw.group) / length(inla_mod$age.rw.group) 
		if(nrow(temp) != length(label) * expand) stop("The input year_range or year_label does not match the fitted model. Please double check.")
		quants <- NULL
		if(!is.null(inla_mod$age.rw.group)){
			for(i in 1:length(inla_mod$age.groups)){
				where <- (inla_mod$age.rw.group[i] - 1) * length(year_label) + c(1:length(year_label))
				quants <- rbind(quants,  temp[where, ])
			}
		}else{
			quants <- temp
		}
		n <- dim(quants)[1]
		m <- length(unstruct)
		quants <- rbind(quants, getquants(unstruct, lower = lower, upper = upper))
		quants$years <- c(label, label.unstruct)
		if(!is.null(inla_mod$age.rw.group)){
			quants$group <- c(group, rep(NA, m))
		}
	  	quants$years.num <- suppressWarnings(as.numeric(as.character(quants$years)))
		quants$label <- c(rep("RW", n), rep("IID", m))
		quants$is.yearly <- !(quants$years %in% year_label)

	}else if(field == "space"){
		N <- dim(Amat)[1]
		if("region.unstruct" %in% names(inla_mod$fit$marginals.random)){
			struct <- inla_mod$fit$marginals.random$region.struct
			unstruct <- inla_mod$fit$marginals.random$region.unstruct
			group <- rep(c("Besag", "IID"), each = N)
		}else if("region.struct" %in% names(inla_mod$fit$marginals.random)){
			struct <- inla_mod$fit$marginals.random$region.struct[1:N]
			unstruct <- inla_mod$fit$marginals.random$region.struct[(N+1):(N*2)]	
			group <- rep(c("Total", "Besag"), each = N)
		}else{
			stop("No spatial term used in this model.")
		}
		quants <- rbind(getquants(struct, lower = lower, upper = upper), getquants(unstruct, lower = lower, upper = upper))
		quants$region <- rep(colnames(Amat), 2)
		quants$label <- group
	}else if(field == "spacetime"){
		if("region.struct" %in% names(inla_mod$fit$marginals.random)){
			N <- dim(Amat)[1]
		}else{
			stop("No spatial term used in this model.")
		}
		if(is.yearly){
		    label <- c(year_range[1] : year_range[2], year_label)
		  }else{
		    label <- year_label
		 }
		has.random.slope <- sum(grepl("slope", names(inla_mod$fit$summary.random)))
		if(has.random.slope){
			fixed <- c("st.slope.id", "region.int", "time.area")
	        select <- list()
	        for(i in 1:length(fixed)){
	           select[[i]] <- 0
	           names(select)[i] <- fixed[i]
	        }
        	if(is.null(draws)){
				message("AR1 model diagnostics are still experimental, starting posterior sampling...")
	        	sampAll <- INLA::inla.posterior.sample(n = 1000, result = inla_mod$fit, intern = TRUE, selection = select, verbose = FALSE)
	        }else{
	        	sampAll <- draws
	        }
	        fields <- rownames(sampAll[[1]]$latent)        
			A <- inla_mod$fit$.args$data
	        A <- A[, colnames(A) %in% c(fixed, "time.unstruct", "region.struct")]
	        A <- unique(A)
	        AA.loc <- A
	        AA.loc$time.area  <- match(paste0("time.area:", AA.loc$time.area), fields)
	        region.int.index <- (AA.loc$time.unstruct - 1) * max(AA.loc$region.struct) + AA.loc$region.struct
			AA.loc$region.int  <- match(paste0("region.int:", region.int.index), fields)
	        T <- max(A$time.unstruct)
	        AA.loc$ststar <- (AA.loc$time.unstruct - (T + 1)/2) / (T + 1)
	        AA.loc$st.slope  <- match(paste0("st.slope.id:", AA.loc$region.struct), fields)

	        sum <- matrix(NA, dim(AA.loc)[1], length(sampAll))
	        for(i in 1:length(sampAll)){
	    		s1 <- sampAll[[i]]$latent[AA.loc$time.area]  		
	    		s2 <- sampAll[[i]]$latent[AA.loc$region.int]  
	    		s3 <- sampAll[[i]]$latent[AA.loc$st.slope] * AA.loc$ststar
	    		if(sum(!is.na(s1)) == 0){
	    			sum[, i] <- s2 + s3
	    		}else{
	    			sum[, i] <- s1 + s3
	    		}
		    }
			group <- colnames(Amat)[AA.loc$region.struct]
		    label <- label[AA.loc$time.unstruct]
		    struct <- NULL
		    quants <- t(apply(sum, 1, function(x, lower, upper){as.numeric(stats::quantile(x, c(lower, 0.5, upper)))}, lower, upper))
		    quants <- data.frame(quants)
		    colnames(quants) <- c("lower", "median", "upper")

		}else if("region.int" %in% names(inla_mod$fit$marginals.random)){
			 group <- rep(colnames(Amat), length(label))
			 label <- rep(label, each = N)	
			 struct <- inla_mod$fit$marginals.random$region.int	 	
		 }else{
		 	group <- colnames(Amat)[inla_mod$time.area$region_number]	
		 	label  <- label[inla_mod$time.area$time.unstruct]	
			struct <- inla_mod$fit$marginals.random$time.area
		 }

		if(!is.null(struct) && length(struct) != length(label)) stop("The input year_range or year_label does not match the fitted model. Please double check.")
		if(!is.null(struct)) quants <- getquants(struct, lower = lower, upper = upper) 
		quants$years <- label
	  	quants$years.num <- suppressWarnings(as.numeric(as.character(quants$years)))
		quants$region <- group	
		quants$is.yearly <- !(quants$years %in% year_label)

	}else{
		stop("The field argument needs to be space, time, or spacetime.")
	}
  	class(quants) <- c("SUMMERproj", "data.frame")

	return(quants)
}
