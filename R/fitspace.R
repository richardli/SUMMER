#' Fit space-time smoothing models for a generic outcome from complex surveys.
#'
#' This function calculates the direct estimates by region and fit a simple spatial smoothing model to the direct estimates adjusting for survey design.
#' Normal or binary variables are currently supported. For binary variables, the logit transformation is performed on the direct estimates of probabilities, and a Gaussian additive model is fitted on the logit scale using INLA.
#' 
#' The function \code{smoothSurvey} replaces the previous function name \code{fitGeneric} (before version 1.0.0).
#' 
#' @param data data frame with region and strata information.
#' @param geo Deprecated argument from early versions.
#' @param Amat Adjacency matrix for the regions.
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
#' @param formula a string of user-specified random effects model to be used in the INLA call
#' @param timeVar The variable indicating time period. If set to NULL then the temporal model and space-time interaction model are ignored.
#' @param time.model the model for temporal trends and interactions. It can be either "rw1" or "rw2".
#' @param type.st can take values 0 (no interaction), or 1 to 4, corresponding to the type I to IV space-time interaction.
#' @param direct.est data frame of direct estimates, with column names of response and region specified by \code{responseVar}, \code{regionVar}, and \code{timeVar}.  When \code{direct.est} is specified, it overwrites the \code{data} input. 
#' @param direct.est.var the column name corresponding to the variance of direct estimates
#' @param counts.bb data frame of counts of binary responses by cluster, with column names of response, region, stratification within region, and cluster ID specified by \code{responseVar}, \code{strataVar.bb}, and \code{clusterVar.bb}.  When \code{counts} is specified, it overwrites the \code{data} input. 
#' @param strataVar.bb the variable specifying within region stratification variable. This is only used for the BetaBinomial model. 
#' @param clusterVar.bb the variable specifying within cluster ID variable. This is only used for the BetaBinomial model. 
#' @param totalVar.bb the variable specifying total observations in \code{counts}. This is only used for the BetaBinomial model when \code{counts} is specified. 
#' @param weight.strata a data frame with one column corresponding to \code{regionVar}, and columns specifying proportion of each strata for each region. This argument specifies the weights for strata-specific estimates on the probability scale. This is only used for the BetaBinomial model. 
#' @param nsim number of posterior draws to take. This is only used for the BetaBinomial model when \code{weight.strata} is provided. 
#' @param ... additional arguments passed to \code{svydesign} function.
#' 
#' 
#' @return \item{HT}{Direct estimates}
#' \item{smooth}{Smoothed direct estimates}
#' \item{fit}{a fitted INLA object}
#' \item{CI}{input argument}
#' \item{Amat}{input argument}
#' \item{responseType}{input argument}
#' \item{formula}{INLA formula}
#' @seealso \code{\link{getDirectList}}, \code{\link{smoothDirect}}
#' @importFrom stats median quantile sd var aggregate as.formula
#' @author Zehang Richard Li 
#' @examples
#' \dontrun{
#' data(DemoData2)
#' data(DemoMap2)
#' fit0 <- smoothSurvey(data=DemoData2,  
#' Amat=DemoMap2$Amat, responseType="binary", 
#' responseVar="tobacco.use", strataVar="strata", 
#' weightVar="weights", regionVar="region", 
#' clusterVar = "~clustid+id", CI = 0.95)
#' 
#' # Example with region-level covariates
#'  Xmat <- aggregate(age~region, data = DemoData2, FUN = mean)
#'  fit1 <- smoothSurvey(data=DemoData2,  
#'   Amat=DemoMap2$Amat, responseType="binary", 
#'   X = Xmat,
#'   responseVar="tobacco.use", strataVar="strata", 
#'   weightVar="weights", regionVar="region", 
#'   clusterVar = "~clustid+id", CI = 0.95)
#' 
#' # Example with using only direct estimates as input instead of the full data
#' direct <- fit$HT[, c("region", "HT.est", "HT.var")]
#' fit2 <- smoothSurvey(data=NULL, direct.est = direct, 
#'                     Amat=DemoMap2$Amat, regionVar="region",
#'                     responseVar="HT.est", direct.est.var = "HT.var", 
#'                     responseType = "binary")
#' # Check it is the same as fit0
#' plot(fit2$smooth$mean, fit0$smooth$mean)
#' 
#' # Example with using only direct estimates as input, 
#' #   and after transformation into a Gaussian smoothing model
#' # Notice: the output are on the same scale as the input 
#' #   and in this case, the logit estimates.    
#' direct.logit <- fit$HT[, c("region", "HT.logit.est", "HT.logit.var")]
#' fit3 <- smoothSurvey(data=NULL, direct.est = direct.logit, 
#'                Amat=DemoMap2$Amat, regionVar="region",
#'                responseVar="HT.logit.est", direct.est.var = "HT.logit.var",
#'                responseType = "gaussian")
#' # Check it is the same as fit0
#' plot(fit3$smooth$mean, fit0$smooth$logit.mean)
#' 
#' }
#' @export


smoothSurvey <- function(data, geo = NULL, Amat, X = NULL, responseType = c("binary", "gaussian")[1], responseVar, strataVar="strata", weightVar="weights", regionVar="region", clusterVar = "~v001+v002", pc.u = 1, pc.alpha = 0.01, pc.u.phi = 0.5, pc.alpha.phi = 2/3, CI = 0.95, formula = NULL, timeVar = NULL, time.model = c("rw1", "rw2")[1], type.st = 1, direct.est = NULL, direct.est.var = NULL, counts.bb = NULL, strataVar.bb = NULL, clusterVar.bb = NULL, totalVar.bb = NULL, weight.strata = NULL, nsim = 1000, ...){

    svy <- TRUE
    if(is.null(responseVar)){
        stop("Response variable not specified")
    }
    if(is.null(responseType)){
            stop("responseType not specified")
    }  
    is.bb <- FALSE

    if(!is.null(strataVar.bb) || !is.null(clusterVar.bb)){
        message("Fitting beta-binomial model.")
        if(responseType != "binary") stop("beta-binomial model only implemented for binary response.")
        if(!is.null(counts.bb)) data <- counts.bb

        if(is.null(strataVar.bb) || is.na(strataVar.bb)){
            message("Within region stratification variable (strataVar.bb) not defined. Unstratified model is fitted.")
            strataVar.bb <- "strata0"
            data$strata0 <- 1
        }
       
        if(is.null(clusterVar.bb)) stop("Cluster variable (clusterVar.bb) not defined.")

        if(is.null(data)) stop("Survey dataset not defined.")
        data$response0 <- data[, responseVar]
        data$region0 <- as.character(data[, regionVar])
        if(sum(!data$region0 %in% colnames(Amat)) > 0) stop("Exist regions in the data frame but not in Amat.")
        data$region0 <- factor(data$region0, levels = colnames(Amat))
        data$strata0 <- data[, strataVar.bb]
        data$cluster0 <- data[, clusterVar.bb]

        if(is.null(counts.bb)){
            vars <- c("cluster0", "region0", "strata0")

            if(!is.null(timeVar)){
                data$time0 <- data[, timeVar]
                vars <- c(vars, "time0")
            }
            counts <- getCounts(data[, c(vars, "response0")], variables = 'response0', by = vars, drop=TRUE)   
        }else{
            if(!is.null(timeVar)){
                data$time0 <- data[, timeVar]
            }
            if(is.null(totalVar.bb)){
                stop("Which column correspond to cluster total is not specified")
            }
            data$total <- data[, totalVar.bb]
            counts <- data
        } 
        
        is.bb <- TRUE   
        
        if(!is.null(weight.strata)){
            if(sum(!data$strata0 %in% colnames(weight.strata)) > 0) stop("Exist within-area strata (strataVar.bb) not in the weight.strata data frame.")
            stratalist <- unique(data$strata0)
        }else{
            stratalist <- NULL
        }

    }else if(!is.null(direct.est)){
        # TODO: this new input type later can be extended to Fay-Herriot
        message("Using direct estimates as input instead of survey data.")
        data <- direct.est
        data$region0 <- as.character(direct.est[, regionVar])
        if(sum(!data$region0 %in% colnames(Amat)) > 0) stop("Exist regions in the data frame but not in Amat.")
        data$region0 <- factor(data$region0, levels = colnames(Amat))        
        data$response0 <- direct.est[, responseVar]
        if(is.null(direct.est.var)){
            stop("Need to specify column for the variance of direct estimates")
        }
        data$var0 <- direct.est[, direct.est.var]
    }else{
        if(!is.data.frame(data)){
            stop("Input data needs to be a data frame")
        }
        if (is.null(strataVar)) {
            message("Strata not defined. Ignoring sample design")
            svy <- FALSE
        }
       
        if (is.null(clusterVar)){
            message("cluster not specified. Ignoring sample design")
            svy <- FALSE
        } 
        data$response0 <- data[, responseVar]
        data$region0 <- as.character(data[, regionVar])
        if(sum(!data$region0 %in% colnames(Amat)) > 0) stop("Exist regions in the data frame but not in Amat.")
        data$region0 <- factor(data$region0, levels = colnames(Amat))   
        if(svy){
            data$weights0 <- data[, weightVar]
            data$strata0 <- data[, strataVar]
        }
    }
    if(!is.null(timeVar)) data$time0 <- data[, timeVar]     


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
    
 
    # if(is.null(FUN)){
    if(responseType == "binary"){
        # message("FUN is not specified, default to be expit()")
        FUN <- expit
    }else if(responseType == "gaussian"){
        # message("FUN is not specified, default to be no transformation")
        FUN <- function(x){x}
    }
    # }
    if (!isTRUE(requireNamespace("INLA", quietly = TRUE))) {
    stop("You need to install the packages 'INLA'. Please run in your R terminal:\n  install.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/stable'), dep=TRUE)")
  }
  # If INLA is installed, then attach the Namespace (so that all the relevant functions are available)
  if (isTRUE(requireNamespace("INLA", quietly = TRUE))) {
    if (!is.element("INLA", (.packages()))) {
      attachNamespace("INLA")
    }
  }
  

    hyperpc1 <- list(prec = list(prior = "pc.prec", param = c(pc.u , pc.alpha)))
    hyperpc2 <- list(prec = list(prior = "pc.prec", param = c(pc.u , pc.alpha)),  phi = list(prior = 'pc', param = c(pc.u.phi , pc.alpha.phi)))

    if(sum(!data$region0 %in% colnames(Amat)) > 0){
       stop("Exist regions in data but not in the Amat.")
    }

    if(!is.bb){
        if(svy && is.null(direct.est)){
            design <- survey::svydesign(
                            ids = stats::formula(clusterVar),
                            weights = ~weights0,
                            strata = ~strata0,
                            data = data, 
                            ...)
            if(!is.null(timeVar)){
                mean <- survey::svyby(formula=~response0, by=~region0+time0, design=design, survey::svymean, drop.empty.groups=FALSE)
                time.i <- as.numeric(as.character(mean$time0))
            }else{
                mean <- survey::svyby(formula=~response0, by=~region0, design=design, survey::svymean, drop.empty.groups=FALSE)
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
        }else if(is.null(direct.est)){
            if(!is.null(timeVar)){
                mean <- aggregate(response0 ~ region0+time0, data = data, FUN = function(x){c(mean(x), length(x), sum(x))})    
                  name.i <- mean$region0
                  time.i <- as.numeric(as.character(mean$time0))
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
        }else{
            p.i <- data$response0
            var.i <- data$var0
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
            n <- NA
            y <- NA
            time.i <- as.numeric(as.character(data$time0))
            name.i <- as.character(data$region0)
        }

        dat <- data.frame(
                          HT.est = p.i,
                          HT.var = var.i, 
                          HT.logit.est = ht,
                          HT.logit.var = ht.v, 
                          HT.logit.prec = ht.prec,  
                          n = n, 
                          y = y)   
    }else{
        dat <- counts
        dat$y <- dat$response0
        name.i <- dat$region0
        if(!is.null(timeVar)) time.i <- as.numeric(as.character(dat$time0))
    }

    if(!is.null(timeVar)){
        dat$time <- time.i
        dat$time.struct <- dat$time.unstruct <- time.int <- dat$time - min(dat$time) + 1
    }   
    # make it consistent with map
    regnames <- as.character(name.i)
    # dat <- dat[match(rownames(Amat), regnames), ]
    dat$region <- as.character(name.i)
    dat$region.unstruct <- dat$region.struct <- dat$region.int <-  match(dat$region, rownames(Amat))
    if(!is.bb) dat$HT.logit.est[is.infinite(abs(dat$HT.logit.est))] <- NA  


    if(!is.null(timeVar)){
        dat <- dat[order(dat[, "time.struct"], dat[, "region.struct"]), ]
    }else{
        dat <- dat[order(dat[, "region.struct"]), ]
    }
    

    if(!svy && responseType == "binary" && !is.bb){
        formulatext <- "y ~ 1"
    }else if(!is.bb){
        formulatext <- "HT.logit.est ~ 1"
    }else if(length(unique(data$strata0)) == 1){
        formulatext <- "y ~ 1"
    }else{
        formulatext <- "y ~ strata0 - 1"        
    }
    if(!is.null(X)){
        X <- data.frame(X)        
        fixed <- colnames(X)[-1]
        colnames(X)[1] <- "region"
        if(fixed %in% colnames(dat)){
            message("The following covariates exist in the input data frame. They are replaced with region-level covariates provided in X: ", fixed[fixed %in% colnames(dat)])
            dat <- dat[, !colnames(dat) %in% fixed]
        }
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

    if(is.bb){
        fit <- INLA::inla(formula, family="betabinomial", Ntrials=dat$total, control.compute = list(dic = T, mlik = T, cpo = T, config = TRUE), data = dat, control.predictor = list(compute = TRUE),  lincomb = NULL, quantiles = c((1-CI)/2, 0.5, 1-(1-CI)/2)) 

    }else if(!svy && responseType == "binary"){
         fit <- INLA::inla(formula, family="binomial", Ntrials=n, control.compute = list(dic = T, mlik = T, cpo = T), data = dat, control.predictor = list(compute = TRUE),  lincomb = NULL, quantiles = c((1-CI)/2, 0.5, 1-(1-CI)/2)) 
    
    # Weighted estimates
    }else if(!svy && responseType == "gaussian"){
        fit <- INLA::inla(formula, family="gaussian", control.compute = list(dic = T, mlik = T, cpo = T), data = dat, control.predictor = list(compute = TRUE), control.family = list(hyper= list(prec = list(initial= log(1), fixed= TRUE))), lincomb = NULL, quantiles = c((1-CI)/2, 0.5, 1-(1-CI)/2)) 
    }else{
         fit <- INLA::inla(formula, family = "gaussian", control.compute = list(dic = T, mlik = T, cpo = T), data = dat, control.predictor = list(compute = TRUE), control.family = list(hyper= list(prec = list(initial= log(1), fixed= TRUE))), scale = dat$HT.logit.prec,  lincomb = NULL, quantiles = c((1-CI)/2, 0.5, 1-(1-CI)/2)) 
    }
    

    
    n <- max(dat$region.struct) 
    temp <- NA
    if(!is.bb){
        if(!is.null(timeVar)) {
            n <- n * max(dat$time.struct)
            temp <- dat$time[1:n]
        }
        proj <- data.frame(region=dat$region[1:n], time = temp, mean=NA, var=NA, median=NA, lower=NA, upper=NA, logit.mean=NA, logit.var=NA, logit.median=NA, logit.lower=NA, logit.upper=NA)        
    }else{
        if(!is.null(timeVar)) {
            temp <- 1:max(dat$time.struct)
        }
        proj <- expand.grid(region = colnames(Amat), time = temp, strata = unique(dat$strata0))
        proj <- cbind(proj, data.frame(mean=NA, var=NA, median=NA, lower=NA, upper=NA, logit.mean=NA, logit.var=NA, logit.median=NA, logit.lower=NA, logit.upper=NA)) 
    }
  

  
    for(i in 1:dim(proj)[1]){
        if(!is.bb){
            tmp <- matrix(INLA::inla.rmarginal(1e5, fit$marginals.linear.predictor[[i]]))
        }else{
            if(is.null(timeVar)){
                which <- which(dat$region == proj$region[i] & dat$strata0 == proj$strata[i])[1]
            }else{
                which <- which(dat$region == proj$region[i] & dat$strata0 == proj$strata[i] &
                               dat$time == proj$time[i])[1] 
            }
            if(is.na(which)) next
            tmp <- matrix(INLA::inla.rmarginal(1e5, fit$marginals.linear.predictor[[which]]))
        }
        if(!svy && responseType == "binary"){
           tmp2 <- tmp
            proj[i, "logit.mean"] <- mean(tmp2)
            proj[i, "logit.var"] <- var(tmp2)
            proj[i, "logit.lower"] <- quantile(tmp2, (1-CI)/2)
            proj[i, "logit.upper"] <- quantile(tmp2, 1-(1-CI)/2)
            proj[i, "logit.median"] <- median(tmp2)

            proj[i, "mean"] <- fit$summary.fitted.values[i, "mean"]
            proj[i, "var"] <- fit$summary.fitted.values[i, "sd"]^2
            proj[i, "lower"] <- fit$summary.fitted.values[i, 3]
            proj[i, "upper"] <- fit$summary.fitted.values[i, 5]
            proj[i, "median"] <- fit$summary.fitted.values[i, "0.5quant"]  
        }else{
            tmp2 <- apply(tmp, 2, FUN)
            proj[i, "mean"] <- mean(tmp2)
            proj[i, "var"] <- var(tmp2)
            proj[i, "lower"] <- quantile(tmp2, (1-CI)/2)
            proj[i, "upper"] <- quantile(tmp2, 1-(1-CI)/2)
            proj[i, "median"] <- median(tmp2)

            if(responseType == "binary"){
                proj[i, "logit.mean"] <- fit$summary.fitted.values[i, "mean"]
                proj[i, "logit.var"] <- fit$summary.fitted.values[i, "sd"]^2
                proj[i, "logit.lower"] <- fit$summary.fitted.values[i, 3]
                proj[i, "logit.upper"] <- fit$summary.fitted.values[i, 5]
                proj[i, "logit.median"] <- fit$summary.fitted.values[i, "0.5quant"]  
            }
        }
      
    }

   if(is.bb && !is.null(weight.strata) && length(unique(data$strata0)) > 1){
     proj.agg <- expand.grid(region = colnames(Amat), time = temp)
     proj.agg <- cbind(proj.agg, data.frame(mean=NA, var=NA, median=NA, lower=NA, upper=NA, logit.mean=NA, logit.var=NA, logit.median=NA, logit.lower=NA, logit.upper=NA)) 
     sampAll <- INLA::inla.posterior.sample(n = nsim, result = fit, intern = TRUE)

     for(i in 1:dim(proj.agg)[1]){
        if(is.null(timeVar)){
            which <- which(weight.strata[, regionVar] == proj.agg$region[i])
        }else{
            which <-  which(weight.strata[, regionVar] == proj.agg$region[i] &
                            weight.strata[, timeVar] == proj.agg$time[i]) 
        }
        draws <- rep(0, nsim)
        for(j in 1:length(sampAll)){
            r <- sampAll[[j]]$latent[paste0("region.struct:", match(proj.agg$region[i], colnames(Amat))), ]
            if(!is.null(timeVar)) r <- r + sampAll[[j]]$latent[paste0("time.struct:", proj.agg$time[i]), ]

            # only handling region-level covariates  
            if(!is.null(X)){
                for(xx in fixed){
                    slope = sampAll[[j]]$latent[paste0(xx, ":1"), ]
                    r <- r + slope * dat[which(dat$region == proj.agg$region[i])[1], xx]
                }
            } 
            for(s in stratalist){
                intercept = sampAll[[j]]$latent[paste0("strata0", s, ":1"), ]
                draws[j] <- draws[j] + expit(r + intercept) * weight.strata[which, s]
            } 
        }
        proj.agg[i, "mean"] <- mean(draws)
        proj.agg[i, "var"] <- var(draws)
        proj.agg[i, "lower"] <- quantile(draws, (1 - CI)/2)
        proj.agg[i, "upper"] <- quantile(draws, 1 - (1 - CI)/2)
        proj.agg[i, "median"] <- median(draws)
        proj.agg[i, "logit.mean"] <- mean(logit(draws))
        proj.agg[i, "logit.var"] <- var(logit(draws))
        proj.agg[i, "logit.lower"] <- quantile(logit(draws), (1-CI)/2)
        proj.agg[i, "logit.upper"] <- quantile(logit(draws), 1-(1-CI)/2)
        proj.agg[i, "logit.median"] <- median(logit(draws))

     }
    }else{
        proj.agg <- NULL
    }

   if(!is.bb){
       # organize output nicer 
       HT <- dat[, !colnames(dat) %in% c("region.struct", "region.unstruct", "region.int", "time.struct", "time.unstruct", "time.int")]
       HT <- cbind(region=HT[, "region"], HT[, colnames(HT) != "region"])
       if("time" %in% names(HT)){
           HT <- cbind(region=HT[, "region"], 
                       time=HT[, "time"], 
                       HT[, !(colnames(HT) %in% c("region", "time"))])
       }
       for(i in 1:dim(HT)[2]) HT[is.nan(HT[,i]), i] <- NA
       if(sum(!is.na(HT$n)) == 0) HT <- HT[, colnames(HT) != "n"]
       if(sum(!is.na(HT$y)) == 0) HT <- HT[, colnames(HT) != "y"]
       if(sum(!is.na(proj$time)) == 0) proj <- proj[, colnames(proj) != "time"]
       if(responseType == "gaussian"){
            HT <- HT[, !colnames(HT) %in% c("HT.logit.est", "HT.logit.var", "HT.logit.prec")]
            proj <- proj[, !colnames(HT) %in% c("logit.mean", "logit.var", "logit.median", "logit.lower", "logit.upper")]
       }    
   }else{
        HT <- NULL
   }


   return(list(HT = HT,
               smooth = proj, 
               smooth.overall = proj.agg, 
               fit = fit, 
               CI = CI,
               Amat = Amat,
               responseType = responseType,
               formula = formula))
}


#' @export
#' @rdname smoothSurvey
fitGeneric <-  smoothSurvey