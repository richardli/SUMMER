#' Fit INLA models to direct estimators.
#' 
#' 
#' 
#'
#' @param data Combined dataset
#' @param geo Geo file
#' @param Amat Adjacency matrix for the regions
#' @param formula INLA formula. Defaults to RW2, ICAR, IID time, IID, region, IID survey effect, IID time-region interaction, IID survey-region interaction, and IID survey-time-region interaction. 
#' @param year_names string vector of year names
#' @param na.rm Logical indicator of whether to remove rows with NA values in the data. Default set to TRUE.
#' @param redo.prior Logical indicator of whether to re-estimate hyper parameters
#' @param priors priors from simhyper
#' @param useHyper option to manually set all hyper priors
#' @param a.iid hyper parameter for i.i.d random effects, only need if \code{useHyper = TRUE}
#' @param b.iid hyper parameter for i.i.d random effects, only need if \code{useHyper = TRUE}
#' @param a.rw1 hyper parameter for RW1 random effects, only need if \code{useHyper = TRUE}
#' @param b.rw1 hyper parameter for RW1 random effects, only need if \code{useHyper = TRUE}
#' @param a.rw2 hyper parameter for RW2 random effects, only need if \code{useHyper = TRUE}
#' @param b.rw2 hyper parameter for RW2 random effects, only need if \code{useHyper = TRUE}
#' @param a.icar hyper parameter for ICAR random effects, only need if \code{useHyper = TRUE}
#' @param b.icar hyper parameter for ICAR random effects, only need if \code{useHyper = TRUE}
#' @seealso \code{\link{countrySummary}}
#' @return INLA model fit using the provided formula, country summary data, and geographic data
#' @examples
#' \dontrun{
#' data(Uganda)
#' data(UgandaMap)
#' geo <- UgandaMap$geo
#' mat <- UgandaMap$Amat
#' years <- c("85-89", "90-94", "95-99", "00-04", "05-09", "10-14")
#' 
#' # Get direct estimates
#' u5m <- countrySummary_mult(births = Uganda, years = years, idVar = "id", 
#' regionVar = "region", timeVar = "time", clusterVar = "~clustid+id", 
#' ageVar = "age", weightsVar = "weights", geo.recode = NULL)
#' 
#' # Get hyper priors
#' priors <- simhyper(R = 2, nsamp = 1e+05, nsamp.check = 5000, Amat = mat)
#' 
#' # Fit INLA models
#' data <- data[data$region %in% c("central","eastern","northern","western"),]
#' inla_model <- fitINLA(data = data, geo = geo, Amat = mat, year_names = years, priors = priors)
#' }
#' 
#' @export
fitINLA <- function(data, Amat, geo, formula = NULL, year_names, na.rm = TRUE, redo.prior = FALSE, priors = NULL, useHyper = FALSE, a.iid = NULL, b.iid = NULL, a.rw1 = NULL, b.rw1 = NULL, a.rw2 = NULL, b.rw2 = NULL, a.icar = NULL, b.icar = NULL) {
    
    
    data$region_num <- match(data$region, c("All", geo$NAME_final)) - 1
    
    
    data <- data[which(data$region != "All"), ]
    #################################################################### Re-calculate hyper-priors
    if (redo.prior) {
        priors <- simhyper(R = 2, nsamp = 1e+05, nsamp.check = 5000, Amat = Amat, nperiod = length(year_names))
    }
    if (!useHyper) {
      a.iid <- priors$a.iid
      b.iid <- priors$b.iid
      a.rw1 <- priors$a.rw1
      b.rw1 <- priors$b.rw1
      a.rw2 <- priors$a.rw2
      b.rw2 <- priors$b.rw2
      a.icar <- priors$a.icar
      b.icar <- priors$b.icar
    }

    
    #################################################################### # remove NA rows? e.g. if no 10-14 available
    if (na.rm) {
        na.count <- apply(data, 1, function(x) {
            length(which(is.na(x)))
        })
        to_remove <- which(na.count == 6)
        if (length(to_remove) > 0) 
            data <- data[-to_remove, ]
    }
    #################################################################### get the list of region and numeric index in one data frame
    region_names <- colnames(Amat)  #names(table(data$region))
    # reorder as in Amat region_names <- region_names[match(colnames(Amat), region_names)]
    region_count <- length(region_names)
    regions <- data.frame(region = region_names, region_number = seq(1, region_count))
    
    # -- merging in the alphabetical region number 1:21 -- #
    dat <- merge(data, regions, by = "region")
    
    # -- creating IDs for the spatial REs -- #
    dat$region.struct <- dat$region.unstruct <- dat$region_number
    
    ################################################################### get the lsit of region and numeric index in one data frame
    year_count <- length(year_names)
    years <- data.frame(year = year_names, year_number = seq(1, year_count))
    
    # -- creating IDs for the temporal REs -- #
    dat$time.unstruct <- dat$time.struct <- years[match(dat$years, years[, 1]), 2]
    
    ################################################################## get the number of surveys
    survey_count <- length(table(data$survey))
    ################################################################## -- these are the time X survey options -- #
    x <- expand.grid(1:year_count, 1:survey_count)
    survey.time <- data.frame(time.unstruct = x[, 1], survey = x[, 2], survey.time = c(1:nrow(x)))
    
    # -- these are the area X survey options -- #
    x <- expand.grid(1:region_count, 1:survey_count)
    survey.area <- data.frame(region_number = x[, 1], survey = x[, 2], survey.area = c(1:nrow(x)))
    
    # -- these are the area X time options -- #
    x <- expand.grid(1:region_count, 1:year_count)
    time.area <- data.frame(region_number = x[, 1], time.unstruct = x[, 2], time.area = c(1:nrow(x)))
    
    # -- these are the area X time X survey options -- #
    x <- expand.grid(1:region_count, 1:year_count, 1:survey_count)
    survey.time.area <- data.frame(region_number = x[, 1], time.unstruct = x[, 2], survey = x[, 3], survey.time.area = c(1:nrow(x)))
    
    # -- merge these all into the data sets -- #
    newdata <- dat
    if (sum(!is.na(dat$survey)) > 0) {
        newdata <- merge(newdata, survey.time, by = c("time.unstruct", "survey"))
        newdata <- merge(newdata, survey.area, by = c("region_number", "survey"))
        newdata <- merge(newdata, survey.time.area, by = c("region_number", "time.unstruct", "survey"))
    }
    newdata <- merge(newdata, time.area, by = c("region_number", "time.unstruct"))
    
    # -- interactions for structured space-time effects -- #
    newdata$idII <- newdata$time.unstruct
    newdata$groupII <- newdata$region.unstruct
    
    newdata$idIII <- newdata$region.unstruct
    newdata$groupIII <- newdata$time.unstruct
    
    newdata$idIV <- newdata$region.unstruct
    newdata$groupIV <- newdata$time.unstruct
    
    #-------------------------------------------------------#
    # ---------------- Beginning of Step 6 -----------------#
    #-------------------------------------------------------#
    # Fit space-time models and compare using LCPO.  this section includes one example all of the other models you might want to
    # consider are included at the end of the doc
    
    ########################## Model Selection ######
    
    # -- subset of not missing and not direct estimate of 0 -- #
    exdat <- newdata
    exdat <- exdat[!is.na(exdat$logit.est) & exdat$logit.est > (-20), ]
    
    
    # -- prior distributions -- #
    
    # -- fitting the model in INLA -- #
    if (is.null(formula)) {
        formula <- logit.est ~ f(survey, model = "iid", param = c(a.iid, b.iid)) + f(survey.time.area, model = "iid", param = c(a.iid, 
            b.iid)) + f(survey.area, model = "iid", param = c(a.iid, b.iid)) + f(survey.time, model = "iid", param = c(a.iid, 
            b.iid)) + f(region.unstruct, model = "iid", param = c(a.iid, b.iid)) + f(region.struct, graph = Amat, model = "besag", 
            param = c(a.icar, b.icar)) + f(time.struct, model = "rw2", param = c(a.rw2, b.rw2)) + f(time.unstruct, model = "iid", 
            param = c(a.iid, b.iid)) + f(time.area, model = "iid", param = c(a.iid, b.iid))
    }
    mod <- formula
    

    if (!isTRUE(requireNamespace("INLA", quietly = TRUE))) {
      stop("You need to install the packages 'INLA'. Please run in your R terminal:\n install.packages('INLA', repos='https://www.math.ntnu.no/inla/R/stable')")
    }
    # If INLA is installed, then attach the Namespace (so that all the relevant functions are available)
    if (isTRUE(requireNamespace("INLA", quietly = TRUE))) {
      if (!is.element("INLA", (.packages()))) {
        attachNamespace("INLA")
      }
      inla11 <- INLA::inla(mod, family = "gaussian", control.compute = list(dic = T, mlik = T, cpo = T), data = exdat, control.predictor = list(compute = TRUE), 
                           control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))), scale = exdat$logit.prec)
    }
    
    return(list(model = mod, fit = inla11, Amat = Amat, newdata = exdat, time = seq(0, year_count - 1), area = seq(0, region_count - 
        1), survey.time = survey.time, survey.area = survey.area, time.area = time.area, survey.time.area = survey.time.area, 
        a.iid = a.iid, b.iid = b.iid, a.rw1 = a.rw1, b.rw1 = b.rw1, a.rw2 = a.rw2, b.rw2 = b.rw2, a.icar = a.icar, b.icar = b.icar))
    
}
