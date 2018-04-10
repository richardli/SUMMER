#' Fit INLA models to perform simple space smoothing.
#'
#' This function calculates the direct estimates by region and fit a simple spatial smoothing model to the direct estimates adjusting for survey design.
#' 
#' Normal or binary variables are currently supported. For binary variables, the logit transformation is performed on the direct estimates of probabilities, and a Gaussian additive model is fitted on the logit scale using INLA.
#' 
#' @param data data frame with region and strata information.
#' @param geo Geo file
#' @param Amat Adjacency matrix for the regions
#' @param family ‘family’ object currently supports 'binomial' (default with logit link function) or 'gaussian'. 
#' @param responseVar the response variable
#' @param strataVar the strata variable
#' @param weightVar the weights variable
#' @param regionVar Variable name for region, typically 'v024', for older surveys might be 'v101'
#' @param clusterVar Variable name for cluster, typically '~v001 + v002'
#' @param hyper the vector of two hyper parameters if specified by user
#' @param CI the desired posterior credible interval to calculate
#' @param FUN the function to transform the posterior draws. Default to be identify function for normal variable and inverse logit transformation for binomial variables
#' 
#' 
#' @return \item{HT}{Direct estimates}
#' \item{smooth}{Spatially smoothed estimates}
#' \item{fit}{a fitted inla object}
#' \item{geo}{input argument}
#' \item{Amat}{input argument}
#' \item{CI}{input argument}
#' \item{family}{input argument}
#' \item{FUN}{input argument}
#' @seealso \code{\link{countrySummary_mult}}, \code{\link{fitINLA}}
#' @importFrom stats median quantile sd
#' @examples
#' \dontrun{
#' data(DemoData2)
#' data(DemoMap2)
#' fit <- fitSpace(data=DemoData2, geo=DemoMap2$geo, 
#' Amat=DemoMap2$Amat, family="binomial", 
#' responseVar="tobacco.use", strataVar="strata", 
#' weightVar="weights", regionVar="region", 
#' clusterVar = "~clustid+id", 
#' hyper=NULL, CI = 0.95)
#' }
#' @export


fitSpace <- function(data, geo, Amat, family, responseVar, strataVar="strata", weightVar="weights", regionVar="region", clusterVar = "~v001+v002", hyper=NULL, CI = 0.95, FUN=NULL){

	if(!is.data.frame(data)){
		stop("Input data needs to be a data frame")
	}
	if (is.null(strataVar)) {
        stop("Strata not defined.")
    }
    if(is.null(responseVar)){
    	stop("Response variable not specified")
    }
    if(is.null(family)){
    	stop("family not specified")
    }

    if (is.null(clusterVar)){
        stop("cluster not specified")
    }
    if(is.null(FUN)){
        if(family == "binomial"){
            FUN <- expit
        }else if(family == "gaussian"){
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
    data$weights0 <- data[, weightVar]
    data$strata0 <- data[, strataVar]
    data$region0 <- data[, regionVar]
    if(is.null(hyper)){
        hyper <- simhyper(R = 2, nsamp = 1e+05, nsamp.check = 5000, Amat = Amat, only.iid = TRUE)
        a.iid <- hyper$a.iid
        b.iid <- hyper$b.iid
    }else{
        a.iid <- hyper[1]
        b.iid <- hyper[2]
    }

    if(is.null(colnames(Amat)) || is.null(rownames(Amat))){
        stop("column and row names for Amat needs to be specified to region names.")
    }
    if(sum(colnames(Amat) != rownames(Amat)) > 0){
        stop("column names and row names do not agree in Amat.")
    }
    if(sum(!data$region0 %in% colnames(Amat)) > 0){
       stop("Exist regions in data but not in the Amat.")
    }


    design <- survey::svydesign(
    				ids = stats::formula(clusterVar),
                    weights = ~weights0,
                    strata = ~strata0,
                    data = data)
  
    # svyglm for simple model gives the same results
    # for future extension with more families 
    # tmp <- svyglm(response0~1, design=subset(design, region0 == "northeastern"), family=stats::quasibinomial)
    # summary(tmp)$coef

    mean <- survey::svyby(formula=~response0, by=~region0, design=design, survey::svymean)
    name.i <- mean$region0
    p.i <- mean$response0
    var.i <- mean$se^2
    if(family == "binomial"){
        ht <- log(p.i/(1-p.i))
        ht.v <- var.i/(p.i^2*(1-p.i)^2)
        ht.prec <- 1/ht.v
    }else if(family == "gaussian"){
        ht <- p.i
        ht.v <- var.i
        ht.prec <- 1/ht.v
    }else{
        stop("family argument only supports binomial or gaussian at the time.")
    }

    regnames <- as.character(name.i)
    reg.unstruct <- reg.struct <- 1:length(regnames)
    dat <- data.frame(HT.est = ht, 
                      HT.sd = ht.v ^ 0.5,
                      HT.variance = ht.v,
                      HT.prec = ht.prec,
                      HT.est.original = p.i,
                      HT.variance.original = var.i,
                      reg.unstruct = reg.unstruct, 
                      reg.struct = reg.struct, 
                      region=regnames)
    formula = HT.est ~ 1 + f(reg.unstruct, model = 'iid', param=c(a.iid,b.iid)) + f(reg.struct, graph=Amat, model="besag", param=c(a.iid,b.iid), scale.model = TRUE)

    fit <- INLA::inla(formula, family = "gaussian", control.compute = list(dic = T, mlik = T, cpo = T), data = dat, control.predictor = list(compute = TRUE), control.family = list(hyper= list(prec = list(initial= log(1), fixed= TRUE))), scale = dat$HT.prec,  lincomb = NULL)
    proj <- data.frame(region=regnames, mean.trans=NA, sd.trans=NA, median.trans=NA, lower.trans=NA, upper.trans=NA, mean=NA, sd=NA, median=NA, lower=NA, upper=NA)
    for(i in 1:length(regnames)){
        tmp <- matrix(INLA::inla.rmarginal(1e5, fit$marginals.linear.predictor[[i]]))
        tmp2 <- apply(tmp, 2, FUN)
        proj[i, "mean.trans"] <- mean(tmp2)
        proj[i, "sd.trans"] <- sd(tmp2)
        proj[i, "lower.trans"] <- quantile(tmp2, (1-CI)/2)
        proj[i, "upper.trans"] <- quantile(tmp2, 1-(1-CI)/2)
        proj[i, "median.trans"] <- median(tmp2)

        proj[i, "mean"] <- mean(tmp)
        proj[i, "sd"] <- sd(tmp)
        proj[i, "lower"] <- quantile(tmp, (1-CI)/2)
        proj[i, "upper"] <- quantile(tmp, 1-(1-CI)/2)
        proj[i, "median"] <- median(tmp)
    }

   return(list(HT = dat,
               smooth = proj, 
               fit = fit, 
               CI = CI,
               Amat = Amat,
               geo = geo,
               family = family,
               FUN = FUN))
}