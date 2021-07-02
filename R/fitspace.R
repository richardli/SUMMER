#' Fit space-time smoothing models for a generic outcome from complex surveys.
#'
#' This function calculates the direct estimates by region and fit a simple spatial smoothing model to the direct estimates adjusting for survey design.
#' Normal or binary variables are currently supported. For binary variables, the logit transformation is performed on the direct estimates of probabilities, and a Gaussian additive model is fitted on the logit scale using INLA.
#' 
#' The function \code{smoothSurvey} replaces the previous function name \code{fitGeneric} (before version 1.0.0).
#' 
#' @param data The input data frame. The input data  with column of the response variable (\code{responseVar}), region ID (\code{regionVar}), stratification within region (\code{strataVar}), and cluster ID (\code{clusterVar}).
#' \itemize{ 
#' \item For area-level model, the data frame consist of survey observations and corresponding survey weights (\code{weightVar}). 
#' \item For unit-level model and \code{is.agg = FALSE}, the data frame should consist of aggregated counts by clusters (for binary responses), or any cluster-level response (for continuous response). For binary response (\code{responseType = 'binary'}), the beta-binomial model will be fitted for cluster-level counts. For continuous response (\code{responseType = 'gaussian'}), a Gaussian smoothing model will be fitted on the cluster-level response. 
#' \item For unit-level model and \code{is.agg = TRUE}, the data frame should be the same as in the area-level model. For binary response (\code{responseType = 'binary'}), the beta-binomial model will be fitted for cluster-level counts aggregated internally. For continuous response (\code{responseType = 'gaussian'}), the nested error model will be fitted on unit-level response.
#' }
#' @param geo Deprecated argument from early versions.
#' @param Amat Adjacency matrix for the regions. If set to NULL, the IID spatial effect will be used.
#' @param X Areal covariates data frame. One of the column name needs to match the \code{regionVar} specified in the function call, in order to be linked to the data input. Currently only supporting time-invariant region-level covariates.
#' @param X.unit Column names of unit-level covariates. When \code{X.unit} is specified, a nested error model will be fitted with unit-level IID noise, and area-level predictions are produced by plugging in the covariate specified in the \code{X} argument. When \code{X} is not specified, the empirical mean of each covariate will be used. This is only implemented for continuous response with the Gaussian likelihood model and unit-level model. 
#' @param responseType Type of the response variable, currently supports 'binary' (default with logit link function) or 'gaussian'. 
#' @param responseVar the response variable
#' @param strataVar the strata variable used in the area-level model. 
#' @param weightVar the weights variable
#' @param regionVar Variable name for region.
#' @param clusterVar Variable name for cluster. For area-level model, this should be a formula for cluster in survey design object, e.g., '~clusterID + householdID'. For unit-level model, this should be the variable name for cluster unit.
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
#' @param direct.est.var the column name corresponding to the variance of direct estimates.
#' @param is.unit.level logical indicator of whether unit-level model is fitted instead of area-level model. 
#' @param is.agg logical indicator of whether the input is at the aggregated counts by cluster. Only used for unit-level model and binary response variable.
#' @param strataVar.within the variable specifying within region stratification variable. This is only used for the unit-level model. 
#' @param totalVar the variable specifying total observations in \code{counts}. This is only used for the unit-level model when \code{counts} is specified. 
#' @param weight.strata a data frame with one column corresponding to \code{regionVar}, and columns specifying proportion of each strata for each region. This argument specifies the weights for strata-specific estimates. This is only used for the unit-level model. 
#' @param nsim number of posterior draws to take. This is only used for the unit-level model when \code{weight.strata} is provided. 
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
#' ##
#' ## 1. Area-level model with binary response
#' ##
#' 
#' data(DemoData2)
#' data(DemoMap2)
#' fit0 <- smoothSurvey(data=DemoData2,  
#' Amat=DemoMap2$Amat, responseType="binary", 
#' responseVar="tobacco.use", strataVar="strata", 
#' weightVar="weights", regionVar="region", 
#' clusterVar = "~clustid+id", CI = 0.95)
#' summary(fit0)
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
#' direct <- fit0$HT[, c("region", "HT.est", "HT.var")]
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
#' direct.logit <- fit0$HT[, c("region", "HT.logit.est", "HT.logit.var")]
#' fit3 <- smoothSurvey(data=NULL, direct.est = direct.logit, 
#'                Amat=DemoMap2$Amat, regionVar="region",
#'                responseVar="HT.logit.est", direct.est.var = "HT.logit.var",
#'                responseType = "gaussian")
#' # Check it is the same as fit0
#' plot(fit3$smooth$mean, fit0$smooth$logit.mean)
#' 
#' # Example with non-spatial smoothing using IID random effects
#' fit4 <- smoothSurvey(data=DemoData2, responseType="binary", 
#'        responseVar="tobacco.use", strataVar="strata", 
#'        weightVar="weights", regionVar="region", 
#'        clusterVar = "~clustid+id", CI = 0.95)
#' 
#' # Using the formula argument, further customizations can be added to the 
#' #  model fitted. For example, we can fit the Fay-Harriot model with 
#' #  IID effect instead of the BYM2 random effect as follows.
#' #  The "region.struct" and "hyperpc1" are picked to match internal object 
#' #  names. Other object names can be inspected from the source of smoothSurvey.
#' fit5 <- smoothSurvey(data=DemoData2,  
#'        Amat=DemoMap2$Amat, responseType="binary", 
#'        formula = "f(region.struct, model = 'iid', hyper = hyperpc1)",
#'        pc.u = 1, pc.alpha = 0.01,
#'        responseVar="tobacco.use", strataVar="strata", 
#'        weightVar="weights", regionVar="region", 
#'        clusterVar = "~clustid+id", CI = 0.95)
#' # Check it is the same as fit4, notice the region order may be different
#' regions <- fit5$smooth$region
#' plot(fit4$smooth[match(regions, fit4$smooth$region),]$logit.mean, fit5$smooth$logit.mean)
#' 
#' ##
#' ## 2. Unit-level model with binary response  
#' ##
#' 
#' # For unit-level models, we need to create stratification variable within regions
#' data <- DemoData2
#' data$urbanicity <- "rural"
#' data$urbanicity[grep("urban", data$strata)] <- "urban"
#' 
#' # Beta-binomial likelihood is used in this model
#' fit6 <- smoothSurvey(data=data, 
#'   Amat=DemoMap2$Amat, responseType="binary", 
#'   X = Xmat, is.unit.level = TRUE,
#'   responseVar="tobacco.use", strataVar.within = "urbanicity", 
#'   regionVar="region", clusterVar = "clustid", CI = 0.95)
#' 
#' # We may use aggregated PSU-level counts as input as well
#' #    in the case of modeling a binary outcome 
#' data.agg <- aggregate(tobacco.use~region + urbanicity + clustid, 
#'                       data = data, FUN = sum)
#' data.agg.total <- aggregate(tobacco.use~region + urbanicity + clustid, 
#'                       data = data, FUN = length)
#' colnames(data.agg.total)[4] <- "total"
#' data.agg <- merge(data.agg, data.agg.total)
#' head(data.agg)
#' 
#' fit7 <- smoothSurvey(data=data.agg, 
#'   Amat=DemoMap2$Amat, responseType="binary", 
#'   X = Xmat, is.unit.level = TRUE, is.agg = TRUE,
#'   responseVar = "tobacco.use", strataVar.within = "urbanicity", 
#'   totalVar = "total", regionVar="region", clusterVar = "clustid", CI = 0.95)
#' 
#' # Check it is the same as fit6
#' plot(fit6$smooth$mean, fit7$smooth$mean)  
#' 
#' ##
#' ## 3. Area-level model with continuous response
#' ##
#' 
#' # The smoothing model is the same as area-level model with binary response
#' #  the continuous direct estimates are smoothed instead of 
#' #  their logit-transformed versions for binary response.
#' fit8 <- smoothSurvey(data=DemoData2, Amat=DemoMap2$Amat, 
#'        responseType="gaussian", responseVar="age", strataVar="strata", 
#'        weightVar="weights", regionVar="region", 
#'        pc.u.phi = 0.5, pc.alpha.phi = 0.5,
#'        clusterVar = "~clustid+id", CI = 0.95)
#' 
#' ##
#' ## 4. Unit-level model with continuous response  
#' ##    (or nested error models)
#' 
#' # The unit-level model assumes for each of the i-th unit,
#' #    Y_{i} ~ intercept + region_effect + IID_i
#' #    where IID_i is the error term specific to i-th unit
#' 
#' # When more than one level of cluster sampling is carried out, 
#' #   they are ignored here. Only the input unit is considered.
#' #   So here we do not need to specify clusterVar any more. 
#' fit9 <- smoothSurvey(data= data, 
#'   Amat=DemoMap2$Amat, responseType="gaussian", 
#'   is.unit.level = TRUE, responseVar="age", strataVar.within = NULL,
#'   regionVar="region", clusterVar = NULL, CI = 0.95)
#' 
#' # To compare, we may also model PSU-level responses. As an illustration, 
#' data.median <- aggregate(age~region + urbanicity + clustid, 
#'                       data = data, FUN = median)
#' 
#' fit10 <- smoothSurvey(data= data.median, 
#'   Amat=DemoMap2$Amat, responseType="gaussian", 
#'   is.unit.level = TRUE, responseVar="age", strataVar.within = NULL,
#'   regionVar="region", clusterVar = "clustid", CI = 0.95)
#' 
#' 
#' # To further incorporate within-area stratification
#' 
#' fit11 <- smoothSurvey(data = data, 
#'   Amat = DemoMap2$Amat, responseType = "gaussian", 
#'   is.unit.level = TRUE, responseVar="age", strataVar.within = "urbanicity",
#'   regionVar = "region", clusterVar = NULL, CI = 0.95)  
#' 
#' # Notice the usual output is now stratified within each region
#' # The aggregated estimates require strata proportions for each region
#' # For illustration, we set strata population proportions below
#' prop <- data.frame(region = unique(data$region), 
#'                             urban = 0.3, 
#'                             rural = 0.7)
#' fit12 <- smoothSurvey(data=data, 
#'   Amat=DemoMap2$Amat, responseType="gaussian", 
#'   is.unit.level = TRUE, responseVar="age", strataVar.within = "urbanicity",
#'   regionVar="region", clusterVar = NULL, CI = 0.95,
#'   weight.strata = prop)  
#' 
#' # aggregated outcome
#' head(fit12$smooth.overall)
#' 
#' # Compare aggregated outcome with direct aggregating posterior means. 
#' # There could be small differences if only 1000 posterior draws are taken.
#' est.urb <- subset(fit11$smooth, strata == "urban")
#' est.rural <- subset(fit11$smooth, strata == "rural")
#' est.mean.post <- est.urb$mean * 0.3 + est.rural$mean * 0.7
#' plot(fit12$smooth.overall$mean, est.mean.post)
#' 
#' 
#' ##
#' ## 6. Unit-level model with continuous response and unit-level covariate 
#' ## 
#' 
#' # For area-level prediction, area-level covariate mean needs to be  
#' #   specified in X argument. And unit-level covariate names are specified
#' #   in X.unit argument.
#' 
#' set.seed(1)
#' sim <- data.frame(region = rep(c(1, 2, 3, 4), 1000),
#'                    X1 = rnorm(4000), X2 = rnorm(4000))
#' Xmean <- aggregate(.~region, data = sim, FUN = sum)
#' sim$Y <- rnorm(4000, mean = sim$X1 + 0.3 * sim$X2 + sim$region)
#' samp <- sim[sample(1:4000, 20), ]
#' fit.sim <- smoothSurvey(data=samp , 
#'                   X.unit = c("X1", "X2"),
#'                   X = Xmean, Amat=NULL, responseType="gaussian", 
#'                   is.unit.level = TRUE, responseVar="Y", regionVar = "region",  
#'                   pc.u = 1, pc.alpha = 0.01, CI = 0.95) 
#' 
#' }
#' @export


smoothSurvey <- function(data, geo = NULL, Amat = NULL, X = NULL, X.unit = NULL, responseType = c("binary", "gaussian")[1], responseVar, strataVar="strata", weightVar="weights", regionVar="region", clusterVar = "~v001+v002", pc.u = 1, pc.alpha = 0.01, pc.u.phi = 0.5, pc.alpha.phi = 2/3, CI = 0.95, formula = NULL, timeVar = NULL, time.model = c("rw1", "rw2")[1], type.st = 1, direct.est = NULL, direct.est.var = NULL, is.unit.level = FALSE, is.agg = FALSE, strataVar.within = NULL,  totalVar = NULL, weight.strata = NULL, nsim = 1000, ...){

    svy <- TRUE
    if(is.null(responseVar)){
        stop("Response variable not specified")
    }
    if(is.null(responseType)){
        stop("responseType not specified")
    }  
    if(is.null(Amat)){
        message("No spatial adjacency matrix is specified. IID area random effect is used.")
    }
    is.iid.space <- FALSE
    responseType <- tolower(responseType)
    
    # there is no aggregated input for Gaussian models
    if(responseType == "gaussian") is.agg <- FALSE

    if(responseType == "binary" && !is.null(X.unit)){
        X.unit <- NULL
        warning("Unit-level covariates not implemented for binary response variable. Set X.unit = NULL.")
    }
    if(is.unit.level == FALSE && !is.null(X.unit)){
        X.unit <- NULL
        warning("Area-level model is fitted. Set X.unit = NULL.")
    }
    if(is.agg && !is.null(X.unit)){
        X.unit <- NULL
        warning("Unit-level covariates cannot be used when input data is aggregated to cluster level (is.agg = TRUE). Set X.unit = NULL.")
    }

    # whether we are fitting nested error model
    # For future update, we can also include two-fold nested error model by adding PSU-level effect
    if(responseType == "gaussian" && is.unit.level){
        is.nested <- TRUE
        clusterVar <- NULL
    }else{
        is.nested <- FALSE
    }

    if(is.nested && !is.null(X.unit) && !is.null(timeVar)){
        stop("Unit-level nested error model with covariates is not implemented with temporal components yet.")
    }
    if(is.nested && !is.null(X.unit) && !is.null(strataVar.within)){
        stop("Unit-level nested error model with covariates is not implemented with stratification components yet.")
    } 

    if(is.unit.level){
        
        if(responseType == "binary"){
            message("Fitting unit-level model.")
            msg <- "Unit-level model"
        }else{
            message("Fitting unit-level nested error model.")
            msg <- "Unit-level nested error model"
        }
        if(is.null(strataVar.within)){
            message("Within region stratification variable (strataVar.within) not defined. Unstratified model is fitted.")
            strataVar.within <- "strata0"
            data$strata0 <- 1
        }
        
        if(is.null(clusterVar) && !is.nested) stop("Cluster variable (clusterVar) not defined.")

        if(is.null(data)) stop("Survey dataset not defined.")
        data$response0 <- data[, responseVar]
        data$region0 <- as.character(data[, regionVar])
        if(is.null(Amat)){
            regions <- unique(data$region0)
            Amat <- matrix(0, length(regions), length(regions))
            colnames(Amat) <- rownames(Amat) <- regions
            is.iid.space <- TRUE
        }
        if(sum(!data$region0 %in% colnames(Amat)) > 0){
            stop("Exist regions in the data frame but not in Amat.")
        }
        data$region0 <- factor(data$region0, levels = colnames(Amat))
        data$strata0 <- data[, strataVar.within]
        if(!is.null(clusterVar)) data$cluster0 <- data[, clusterVar]

        # check column names of covariates
        if(!is.null(X.unit)){
            if(sum(X.unit %in% colnames(data)) != length(X.unit)){
                stop("Exist columns specified X.unit but not in the data.")
            }
            if(!is.null(X) && sum(X.unit %in% colnames(X)) != length(X.unit)){
                stop("Exist columns specified X.unit but not in X.")
            }
            # if X not specified, using empirical mean
            if(is.null(X)){
                data.sub <- data[, c(X.unit, "region0")]
                X <- aggregate(.~region0, data = data.sub, FUN = function(x){mean(x, na.rm = TRUE)})
            }
        }

        # If input is not aggregated, and not fitting a nested model, then aggregate
        if(!is.agg && !is.nested){
            vars <- c("cluster0", "region0", "strata0")

            if(!is.null(timeVar)){
                data$time0 <- data[, timeVar]
                vars <- c(vars, "time0")
            }
            counts <- getCounts(data[, c(vars, "response0")], variables = 'response0', by = vars, drop=TRUE) 
            # if(responseType == "gaussian"){
            #     stop("Response variable is at SSU level. The function currently only supports PSU-level response.")
            #     counts[, "response0"] <- counts[, "response0"] / counts[, "total"]
            # }  
        # If input is aggregated then define total
        }else if(is.agg){
            if(!is.null(timeVar)){
                data$time0 <- data[, timeVar]
            }
            if(is.null(totalVar) && responseType == "binary"){
                stop("Which column correspond to cluster total is not specified")
            }
            data$total <- data[, totalVar]
            counts <- data
        # if fitting nested model
        }else{
            counts <- data
        }
                
        if(!is.null(weight.strata)){
            if(sum(!data$strata0 %in% colnames(weight.strata)) > 0) stop("Exist within-area strata (strataVar.within) not in the weight.strata data frame.")
            stratalist <- unique(data$strata0)
        }else{
            stratalist <- NULL
        }

    }else if(!is.null(direct.est)){
        msg <- "Area-level model using direct estimates as input"
        message("Using direct estimates as input instead of survey data.")
        data <- direct.est
        data$region0 <- as.character(direct.est[, regionVar])
        if(is.null(Amat)){
            regions <- unique(data$region0)
            Amat <- matrix(0, length(regions), length(regions))
            colnames(Amat) <- rownames(Amat) <- regions
            is.iid.space <- TRUE
        }
        if(sum(!data$region0 %in% colnames(Amat)) > 0){
            stop("Exist regions in the data frame but not in Amat.")
        }
        data$region0 <- factor(data$region0, levels = colnames(Amat))        
        data$response0 <- direct.est[, responseVar]
        if(is.null(direct.est.var)){
            stop("Need to specify column for the variance of direct estimates")
        }
        data$var0 <- direct.est[, direct.est.var]
    }else{
        msg <- "Area-level model using survey data as input"
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
        if(is.null(Amat)){
            regions <- unique(data$region0)
            Amat <- matrix(0, length(regions), length(regions))
            colnames(Amat) <- rownames(Amat) <- regions
            is.iid.space <- TRUE
        }
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
        if(regionVar %in% colnames(X)){
            region.col <- which(colnames(X) == regionVar)

        # for backward compatibility, default to first column being region ID
        }else{
            region.col <- 1
        }
        if(sum(X[, region.col] %in% colnames(Amat) == FALSE) > 0 ||
           sum(colnames(Amat) %in% X[,region.col] == FALSE) > 0){
            stop(paste0(colnames(X)[region.col], "is used as region identifier in the input covariate X. It does not match the region names in the data."))
        }
    }
    
 
    # if(is.null(FUN)){
    if(responseType == "binary"){
        FUN <- expit
    }else if(responseType == "gaussian"){
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

    if(!is.unit.level){
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
                mean <- stats::aggregate(response0 ~ region0+time0, data = data, FUN = function(x){c(mean(x), length(x), sum(x))})    
                  name.i <- mean$region0
                  time.i <- as.numeric(as.character(mean$time0))
                  mean <- data.frame(mean[, -c(1:2)])
            }else{
                mean <- stats::aggregate(response0 ~ region0, data = data, FUN = function(x){c(mean(x), length(x), sum(x))})
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
    if(!is.unit.level) dat$HT.logit.est[is.infinite(abs(dat$HT.logit.est))] <- NA  


    if(!is.null(timeVar)){
        dat <- dat[order(dat[, "time.struct"], dat[, "region.struct"]), ]
    }else{
        dat <- dat[order(dat[, "region.struct"]), ]
    }
    
    # binary non-survey area-level model
    if(!svy && responseType == "binary" && !is.unit.level){
        formulatext <- "y ~ 1"

    # binary and continuous survey area-level model  
    # and continuous non-survey area-level model  
    }else if(!is.unit.level){
        formulatext <- "HT.logit.est ~ 1"

    # unit-level model    
    }else if(length(unique(data$strata0)) == 1){
        formulatext <- "y ~ 1"
    }else{
        formulatext <- "y ~ strata0 - 1"        
    }

    # area-level covariates only
    fixed <- NULL
    if(!is.null(X) && is.null(X.unit)){
        X <- data.frame(X)        
        fixed <- colnames(X)[colnames(X) != regionVar]
        colnames(X)[colnames(X) == regionVar] <- "region"
        if(fixed %in% colnames(dat)){
            message("The following covariates exist in the input data frame. They are replaced with region-level covariates provided in X: ", fixed[fixed %in% colnames(dat)])
            dat <- dat[, !colnames(dat) %in% fixed]
        }
        dat <- merge(dat, X, by = "region", all = TRUE)
        formulatext <- paste(formulatext, " + ", paste(fixed, collapse = " + "))

    # unit-level covariates  
    }else if(is.nested && !is.null(X.unit)){
        formulatext <- paste(formulatext, " + ", paste(X.unit, collapse = " + "))

        X <- data.frame(X)        
        fixed <- colnames(X)[colnames(X) != regionVar]
        colnames(X)[colnames(X) == regionVar] <- "region"
        X <- X[, colnames(X) %in% c("region", X.unit)]
        for(j in colnames(dat)){
            if(!j %in% c("region", X.unit)){
                X[, j] <- dat[match(X$region, dat$region), j]
            }
        }

        X$y <- NA
        dat <- rbind(X, dat)
        out.index <- 1:dim(X)[1]

    # unit-level model without covariates   
    }else if(is.nested){
        X <- data.frame(region = colnames(Amat))        
        for(j in colnames(dat)){
            if(!j %in% c("region")){
                X[, j] <- dat[match(X$region, dat$region), j]
            }
        }

        X$y <- NA
        dat <- rbind(X, dat)
        out.index <- 1:dim(X)[1]

    }

    if(is.null(formula)){
       if(is.iid.space){
            formula <- paste(formulatext, "f(region.struct, model = 'iid', hyper = hyperpc1)", sep = "+")    
       }else{
            formula <- paste(formulatext, "f(region.struct, graph=Amat, model='bym2', hyper = hyperpc2, scale.model = TRUE)", sep = "+")    
       }
       ##
       ##  Constraints are not specified for the interactions.
       ##
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

    if(is.unit.level && responseType == "binary"){
        fit <- INLA::inla(formula, family="betabinomial", Ntrials=dat$total, control.compute = list(dic = T, mlik = T, cpo = T, config = TRUE), data = dat, control.predictor = list(compute = TRUE),  lincomb = NULL, quantiles = c((1-CI)/2, 0.5, 1-(1-CI)/2)) 
    
    }else if(is.unit.level && responseType == "gaussian"){
        fit <- INLA::inla(formula, family="gaussian", control.compute = list(dic = T, mlik = T, cpo = T, config = TRUE), data = dat, control.predictor = list(compute = TRUE),  lincomb = NULL, quantiles = c((1-CI)/2, 0.5, 1-(1-CI)/2))  

    }else if(!svy && responseType == "binary"){
         fit <- INLA::inla(formula, family="binomial", Ntrials=n, control.compute = list(dic = T, mlik = T, cpo = T), data = dat, control.predictor = list(compute = TRUE),  lincomb = NULL, quantiles = c((1-CI)/2, 0.5, 1-(1-CI)/2)) 
    
    }else if(!svy && responseType == "gaussian"){
        fit <- INLA::inla(formula, family="gaussian", control.compute = list(dic = T, mlik = T, cpo = T), data = dat, control.predictor = list(compute = TRUE), control.family = list(hyper= list(prec = list(initial= log(1), fixed= TRUE))), lincomb = NULL, quantiles = c((1-CI)/2, 0.5, 1-(1-CI)/2)) 
    }else{
         fit <- INLA::inla(formula, family = "gaussian", control.compute = list(dic = T, mlik = T, cpo = T), data = dat, control.predictor = list(compute = TRUE), control.family = list(hyper= list(prec = list(initial= log(1), fixed= TRUE))), scale = dat$HT.logit.prec,  lincomb = NULL, quantiles = c((1-CI)/2, 0.5, 1-(1-CI)/2)) 
    }
    

    
    n <- max(dat$region.struct) 
    temp <- NA
    if(!is.unit.level){
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
        if(!is.unit.level){
            tmp <- matrix(INLA::inla.rmarginal(1e5, fit$marginals.linear.predictor[[i]]))
        }else{
            if(is.nested && !is.null(X.unit)){
                which <- which(dat[out.index, ]$region == proj$region[i] & dat[out.index, ]$strata0 == proj$strata[i])[1]
                which <- out.index[which]
            }else if(is.null(timeVar)){
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
                proj[i, "logit.mean"] <-mean(tmp)
                proj[i, "logit.var"] <- var(tmp)
                proj[i, "logit.lower"] <- quantile(tmp, (1-CI)/2)
                proj[i, "logit.upper"] <- quantile(tmp, (1-CI)/2) 
                proj[i, "logit.median"] <- median(tmp)
            }
        }
      
    }

   # Aggregation with posterior samples
   if(is.unit.level && !is.null(weight.strata) && length(unique(data$strata0)) > 1){
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
                if(responseType == "binary"){
                    draws[j] <- draws[j] + expit(r + intercept) * weight.strata[which, s]
                }else if(responseType == "gaussian"){
                    draws[j] <- draws[j] + (r + intercept) * weight.strata[which, s]
                }
            } 
        }
        proj.agg[i, "mean"] <- mean(draws)
        proj.agg[i, "var"] <- var(draws)
        proj.agg[i, "lower"] <- quantile(draws, (1 - CI)/2)
        proj.agg[i, "upper"] <- quantile(draws, 1 - (1 - CI)/2)
        proj.agg[i, "median"] <- median(draws)
        if(responseType == "binary"){
            proj.agg[i, "logit.mean"] <- mean(logit(draws))
            proj.agg[i, "logit.var"] <- var(logit(draws))
            proj.agg[i, "logit.lower"] <- quantile(logit(draws), (1-CI)/2)
            proj.agg[i, "logit.upper"] <- quantile(logit(draws), 1-(1-CI)/2)
            proj.agg[i, "logit.median"] <- median(logit(draws))
        }
     }
    }else{
        proj.agg <- NULL
    }

   if(!is.unit.level){
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

   out <- list(HT = HT,
               smooth = proj, 
               smooth.overall = proj.agg, 
               fit = fit, 
               CI = CI,
               Amat = Amat,
               responseType = responseType,
               formula = formula, 
               msg = msg)
   class(out) <-  "SUMMERmodel.svy"

   return(out)
}


#' @export
#' @rdname smoothSurvey
fitGeneric <-  smoothSurvey