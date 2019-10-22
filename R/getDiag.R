#' Make diagnostic plots
#' 
#' @param inla_mod output from \code{\link{fitINLA}}
#' @param field which random effects to plot. It can be one of the following: space, time, and spacetime.
#' @param year_range Entire range of the years (inclusive) defined in year_label. To be deprecated and imputed from the fitted object in the next version of SUMMER.
#' @param year_label vector of year string vector
#' @param Amat adjacency matrix
#' @param CI Desired level of credible intervals

#' 
#' @return List of diagnostic plots
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
#'   fit1 <- fitINLA(data = data, geo = DemoMap$geo, Amat = DemoMap$Amat, 
#'     year_label = years.all, year_range = c(1985, 2019), 
#'     rw = 2, is.yearly=FALSE, m = 5)
#' random.time <- getDiag(fit1, field = "time", year_label = years.all, #' year_range = c(1985, 2019))
#'   random.space <- getDiag(fit1, field = "space", Amat = DemoMap$Amat)
#'   random.spacetime <- getDiag(fit1, field = "spacetime",
#'    year_label = years, year_range = c(1985, 2019), 
#'    Amat = DemoMap$Amat)
#' }
#' 
#' @export


getDiag <- function(inla_mod, field = c("space", "time", "spacetime")[1], year_range = c(1985, 2019), year_label = c("85-89", "90-94", "95-99", "00-04", "05-09", "10-14", "15-19"), Amat = NULL, CI = 0.95){
	lower <- (1 - CI) / 2
	upper <- 1 - lower
	if(!is.null(inla_mod$year_range)){
		year_range <- inla_mod$year_range
	}
	getquants <- function(mlist, lower, upper){
		quants <- data.frame(matrix(NA, length(mlist), 3))
		for(i in 1:length(mlist)){
			quants[i, ] <- INLA::inla.qmarginal(c(0.025, 0.5, 0.975), mlist[[i]])
		}
		colnames(quants) <- c("lower", "median", "upper")
		return(quants)
	}

	is.yearly = inla_mod$is.yearly

	if(field == "time"){
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
		temp <- getquants(struct)	
		quants <- NULL
		if(!is.null(inla_mod$age.rw.group)){
			for(i in 1:length(inla_mod$age.groups)){
				where <- (inla_mod$age.rw.group[i] - 1) * length(year_label) + c(1:length(year_label))
				quants <- rbind(quants,  temp[where, ])
			}
		}else{
			quants <- getquants(struct)
		}
		n <- dim(quants)[1]
		m <- length(unstruct)
		quants <- rbind(quants, getquants(unstruct))
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
		quants <- rbind(getquants(struct), getquants(unstruct))
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

		 if("region.int" %in% names(inla_mod$fit$marginals.random)){
			 group <- rep(colnames(Amat), length(label))
			 label <- rep(label, each = N)	
			 struct <- inla_mod$fit$marginals.random$region.int	 	
		 }else{
		 	group <- colnames(Amat)[inla_mod$time.area$region_number]	
		 	label  <- label[inla_mod$time.area$time.unstruct]	
			struct <- inla_mod$fit$marginals.random$time.area
		 }
		if(length(struct) != length(label)) stop("The input year_range or year_label does not match the fitted model. Please double check.")
		quants <- getquants(struct) 
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