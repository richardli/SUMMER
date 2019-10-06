#' Fit space-time smoothing models for a generic outcome from complex surveys.
#'
#' This function calculates the direct estimates by region and fit a simple spatial smoothing model to the direct estimates adjusting for survey design.
#' 
#' Normal or binary variables are currently supported. For binary variables, the logit transformation is performed on the direct estimates of probabilities, and a Gaussian additive model is fitted on the logit scale using INLA.
#' 
#' @param data data frame with region and strata information.
#' @param geo Geo file
#' @param Amat Adjacency matrix for the regions
#' @param X Covariate matrix with the first column being the region names. Currently only supporting static region-level covariates.
#' @param responseType Type of the response variable, currently supports 'binary' (default with logit link function) or 'gaussian'. 
#' @param responseVar the response variable
#' @param strataVar the strata variable
#' @param weightVar the weights variable
#' @param regionVar Variable name for region, typically 'v024', for older surveys might be 'v101'
#' @param clusterVar Variable name for cluster, typically '~v001 + v002'
#' @param pc.u hyperparameter U for the PC prior on precisions.
#' @param pc.alpha hyperparameter alpha for the PC prior on precisions.
#' @param pc.u.phi hyperparameter U for the PC prior on the mixture probability phi in BYM2 model.
#' @param pc.alpha.phi hyperparameter alpha for the PC prior on the mixture probability phi in BYM2 model.
#' @param CI the desired posterior credible interval to calculate
#' @param FUN the function to transform the posterior draws. Default to be identify function for normal variable and inverse logit transformation for binomial variables
#' @param formula a string of user-specified random effects model to be used in the INLA call
#' @param  timeVar The variable indicating time period. If set to NULL then the temporal model and space-time interaction model are ignored.
#' @param time.model the model for temporal trends and interactions. It can be either "rw1" or "rw2".
#' @param type.st can take values 0 (no interaction), or 1 to 4, corresponding to the type I to IV space-time interaction.
#' 
#' 
#' @return \item{HT}{Direct estimates}
#' \item{smooth}{Spatially smoothed estimates}
#' \item{fit}{a fitted INLA object}
#' \item{geo}{input argument}
#' \item{Amat}{input argument}
#' \item{CI}{input argument}
#' \item{responseType}{input argument}
#' \item{FUN}{input argument}
#' @seealso \code{\link{getDirectList}}, \code{\link{fitINLA}}
#' @importFrom stats median quantile sd var aggregate as.formula
#' @examples
#' \dontrun{
#' data(DemoData2)
#' data(DemoMap2)
#' fit <- fitGeneric(data=DemoData2, geo=DemoMap2$geo, 
#' Amat=DemoMap2$Amat, responseType="binary", 
#' responseVar="tobacco.use", strataVar="strata", 
#' weightVar="weights", regionVar="region", 
#' clusterVar = "~clustid+id", CI = 0.95)
#' 
#' # Example with region-level covariates
#'  Xmat <- aggregate(age~region, data = DemoData2, FUN = mean)
#'  fit <- fitGeneric(data=DemoData2, geo=DemoMap2$geo, 
#'   Amat=DemoMap2$Amat, responseType="binary", 
#'   X = Xmat,
#'   responseVar="tobacco.use", strataVar="strata", 
#'   weightVar="weights", regionVar="region", 
#'   clusterVar = "~clustid+id", CI = 0.95)
#' }
#' @export


fitGeneric <- function(data, geo, Amat, X = NULL, responseType = c("binary", "gaussian")[1], responseVar, strataVar="strata", weightVar="weights", regionVar="region", clusterVar = "~v001+v002", pc.u = 1, pc.alpha = 0.01, pc.u.phi = 0.5, pc.alpha.phi = 2/3, CI = 0.95, FUN=NULL, formula = NULL, timeVar = NULL, time.model = c("rw1", "rw2")[1], type.st = 1){

    svy <- TRUE
	if(!is.data.frame(data)){
		stop("Input data needs to be a data frame")
	}
	if (is.null(strataVar)) {
        message("Strata not defined. Ignoring sample design")
        svy <- FALSE
    }
    if(is.null(responseVar)){
    	stop("Response variable not specified")
    }
    if(is.null(responseType)){
    	stop("responseType not specified")
    }
    if(is.null(rownames(Amat))){
        stop("Row names of Amat needs to be specified to region names.")
    }
    if(is.null(colnames(Amat))){
        stop("Column names of Amat needs to be specified to region names.")
    }
    if(sum(rownames(Amat) != colnames(Amat)) > 0){
        stop("Row and column names of Amat needs to be the same.")
    }
    if(!is.null(X)){
        if(sum(X[,1] %in% colnames(Amat) == FALSE) > 0 ||
           sum(colnames(Amat) %in% X[,1] == FALSE) > 0) stop("Regions in the X matrix does not match the region names in Amat.")
    }
    
    if (is.null(clusterVar)){
        message("cluster not specified. Ignoring sample design")
        svy <- FALSE
    }
    if(is.null(FUN)){
        if(responseType == "binary"){
            message("FUN is not specified, default to be expit()")
            FUN <- expit
        }else if(responseType == "gaussian"){
            message("FUN is not specified, default to be no transformation")
            FUN <- function(x){x}
        }
    }
    if (!isTRUE(requireNamespace("INLA", quietly = TRUE))) {
    stop("You need to install the packages 'INLA'. Please run in your R terminal:\n install.packages('INLA', repos='https://www.math.ntnu.no/inla/R/stable')")
  }
  # If INLA is installed, then attach the Namespace (so that all the relevant functions are available)
  if (isTRUE(requireNamespace("INLA", quietly = TRUE))) {
    if (!is.element("INLA", (.packages()))) {
      attachNamespace("INLA")
    }
  }


    data$response0 <- data[, responseVar]
    data$region0 <- data[, regionVar]        
    if(svy){
        data$weights0 <- data[, weightVar]
        data$strata0 <- data[, strataVar]
    }
    if(!is.null(timeVar)) data$time0 <- data[, timeVar]

    hyperpc1 <- list(prec = list(prior = "pc.prec", param = c(pc.u , pc.alpha)))
    hyperpc2 <- list(prec = list(prior = "pc.prec", param = c(pc.u , pc.alpha)),  phi = list(prior = 'pc', param = c(pc.u.phi , pc.alpha.phi)))


    if(is.null(colnames(Amat)) || is.null(rownames(Amat))){
        stop("column and row names for Amat needs to be specified to region names.")
    }
    if(sum(colnames(Amat) != rownames(Amat)) > 0){
        stop("column names and row names do not agree in Amat.")
    }
    if(sum(!data$region0 %in% colnames(Amat)) > 0){
       stop("Exist regions in data but not in the Amat.")
    }

    if(svy){
        design <- survey::svydesign(
        				ids = stats::formula(clusterVar),
                        weights = ~weights0,
                        strata = ~strata0,
                        data = data)
        if(!is.null(timeVar)){
            mean <- survey::svyby(formula=~response0, by=~region0+time0, design=design, survey::svymean)
            time.i <- mean$time0
        }else{
            mean <- survey::svyby(formula=~response0, by=~region0, design=design, survey::svymean)
        }
        name.i <- mean$region0            
        p.i <- mean$response0
        var.i <- mean$se^2
        if(responseType == "binary"){
            ht <- log(p.i/(1-p.i))
            ht.v <- var.i/(p.i^2*(1-p.i)^2)
            ht.prec <- 1/ht.v
        }else if(responseType == "gaussian"){
            ht <- p.i
            ht.v <- var.i
            ht.prec <- 1/ht.v
        }else{
            stop("responseType argument only supports binary or gaussian at the time.")
        }
        n <- y <- NA
    }else{
        if(!is.null(timeVar)){
            mean <- aggregate(response0 ~ region0+time0, data = data, FUN = function(x){c(mean(x), length(x), sum(x))})    
              name.i <- mean$region0
              time.i <- mean$time0
              mean <- data.frame(mean[, -c(1:2)])
        }else{
            mean <- aggregate(response0 ~ region0, data = data, FUN = function(x){c(mean(x), length(x), sum(x))})
            name.i <- mean$region0
            mean <- data.frame(mean[, -1])
        }

        colnames(mean) <- c("mean", "n", "y") 
        p.i <- mean$mean
        var.i <- p.i * (1-p.i)/mean$n
        if(responseType == "binary"){
            ht <- log(p.i/(1-p.i))
            ht.v <- var.i/(p.i^2*(1-p.i)^2)
            ht.prec <- 1/ht.v
        }else if(responseType == "gaussian"){
            ht <- p.i
            ht.v <- var.i
            ht.prec <- 1/ht.v
        }else{
            stop("responseType argument only supports binary or gaussian at the time.")
        }
        n <- mean$n
        y <- mean$y
    }

    dat <- data.frame(HT.est = ht, 
                      HT.sd = ht.v ^ 0.5,
                      HT.variance = ht.v,
                      HT.prec = ht.prec,
                      HT.est.original = p.i,
                      HT.variance.original = var.i, 
                      n = n, 
                      y = y)
    if(!is.null(timeVar)){
        dat$time <- time.i
        dat$time.struct <- dat$time.unstruct <- time.int <- dat$time - min(dat$time) + 1
    }
    # make it consistent with map
    regnames <- as.character(name.i)
    # dat <- dat[match(rownames(Amat), regnames), ]
    dat$region <- as.character(name.i)
    dat$region.unstruct <- dat$region.struct <- dat$region.int <-  match(dat$region, rownames(Amat))
    dat$HT.est[is.infinite(abs(dat$HT.est))] <- NA
    if(!is.null(timeVar)){
        dat <- dat[order(dat[, "time.struct"], dat[, "region.struct"]), ]
    }else{
        dat <- dat[order(dat[, "region.struct"]), ]
    }
    

    if(!svy && responseType == "binary"){
        formulatext <- "y ~ 1"
    }else{
        formulatext <- "HT.est ~ 1"
    }
    if(!is.null(X)){
        X <- data.frame(X)        
        fixed <- colnames(X)[-1]
        colnames(X)[1] <- "region"
        dat <- merge(dat, X, by = "region", all = TRUE)
        formulatext <- paste(formulatext, " + ", paste(fixed, collapse = " + "))
    }

    if(is.null(formula)){
        formula <- paste(formulatext, "f(region.struct, graph=Amat, model='bym2', hyper = hyperpc2, scale.model = TRUE)", sep = "+")
       if(!is.null(timeVar)){
            formula <- paste(formula, "f(time.unstruct, model = 'iid', hyper = hyperpc1) + f(time.struct, model=tolower(time.model), hyper = hyperpc1, scale.model = TRUE)", sep = "+")
            if(type.st == 1){
                formula <- paste(formula, "f(region.int, model = 'iid', hyper = hyperpc1, group = time.int, control.group = list(model ='iid', hyper = hyperpc1))", sep = "+")
            }else if(type.st == 2){
                formula <- paste(formula, "f(region.int, model = 'besag', graph = Amat, scale.model = TRUE, param=hyperpc1, group = time.int, control.group = list(model ='iid'))", sep = "+")
            }else if(type.st == 3){
                formula <- paste(formula, "f(region.int, model = 'iid', hyper = hyperpc1, group = time.int, control.group = list(model =tolower(time.model), scale.model = TRUE))", sep = "+")
            }else{
                formula <- paste(formula, "f(region.int, model = 'besag', graph = Amat, scale.model = TRUE, hyper = hyperpc1, group = time.int, control.group = list(model =tolower(time.model), hyper = hyperpc1, scale.model = TRUE))", sep = "+")
            }
        }
        formula <- as.formula(formula)
    }else{
        formula <- as.formula(paste(formulatext, formula, sep = "+"))
    }   

    if(!svy && responseType == "binary"){
         fit <- INLA::inla(formula, family="binomial", Ntrials=n, control.compute = list(dic = T, mlik = T, cpo = T), data = dat, control.predictor = list(compute = TRUE),  lincomb = NULL, quantiles = c((1-CI)/2, 0.5, 1-(1-CI)/2)) 
    
    # Weighted estimates
    }else if(!svy && responseType == "gaussian"){
        fit <- INLA::inla(formula, family="gaussian", control.compute = list(dic = T, mlik = T, cpo = T), data = dat, control.predictor = list(compute = TRUE), control.family = list(hyper= list(prec = list(initial= log(1), fixed= TRUE))), lincomb = NULL, quantiles = c((1-CI)/2, 0.5, 1-(1-CI)/2)) 
    }else{
         fit <- INLA::inla(formula, family = "gaussian", control.compute = list(dic = T, mlik = T, cpo = T), data = dat, control.predictor = list(compute = TRUE), control.family = list(hyper= list(prec = list(initial= log(1), fixed= TRUE))), scale = dat$HT.prec,  lincomb = NULL, quantiles = c((1-CI)/2, 0.5, 1-(1-CI)/2)) 
    }
    

    
    n <- max(dat$region.struct) 
    temp <- NA
    if(!is.null(timeVar)) {
        n <- n * max(dat$time.struct)
        temp <- dat$time[1:n]
    }
    proj <- data.frame(region=dat$region[1:n], time = temp, mean=NA, variance=NA, median=NA, lower=NA, upper=NA, mean.original=NA, variance.original=NA, median.original=NA, lower.original=NA, upper.original=NA)
    for(i in 1:n){
        tmp <- matrix(INLA::inla.rmarginal(1e5, fit$marginals.linear.predictor[[i]]))
        if(!svy && responseType == "binary"){
           tmp2 <- tmp
            proj[i, "mean"] <- mean(tmp2)
            proj[i, "variance"] <- var(tmp2)
            proj[i, "lower"] <- quantile(tmp2, (1-CI)/2)
            proj[i, "upper"] <- quantile(tmp2, 1-(1-CI)/2)
            proj[i, "median"] <- median(tmp2)

            proj[i, "mean.original"] <- fit$summary.fitted.values[i, "mean"]
            proj[i, "variance.original"] <- fit$summary.fitted.values[i, "sd"]^2
            proj[i, "lower.original"] <- fit$summary.fitted.values[i, 3]
            proj[i, "upper.original"] <- fit$summary.fitted.values[i, 5]
            proj[i, "median.original"] <- fit$summary.fitted.values[i, "0.5quant"]  
        }else{
           tmp2 <- apply(tmp, 2, FUN)
            proj[i, "mean.original"] <- mean(tmp2)
            proj[i, "variance.original"] <- var(tmp2)
            proj[i, "lower.original"] <- quantile(tmp2, (1-CI)/2)
            proj[i, "upper.original"] <- quantile(tmp2, 1-(1-CI)/2)
            proj[i, "median.original"] <- median(tmp2)

            proj[i, "mean"] <- fit$summary.fitted.values[i, "mean"]
            proj[i, "variance"] <- fit$summary.fitted.values[i, "sd"]^2
            proj[i, "lower"] <- fit$summary.fitted.values[i, 3]
            proj[i, "upper"] <- fit$summary.fitted.values[i, 5]
            proj[i, "median"] <- fit$summary.fitted.values[i, "0.5quant"]  
        }
      
    }
   return(list(HT = dat[, !colnames(dat) %in% c("region.struct", "region.unstruct", "region.int", "time.struct", "time.unstruct", "time.int")],
               smooth = proj, 
               fit = fit, 
               CI = CI,
               Amat = Amat,
               geo = geo,
               responseType = responseType,
               FUN = FUN, 
               formula = formula))
}