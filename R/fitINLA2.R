#' Cluster-level space-time smoothing models for mortality rates 
#' 
#' The function \code{smoothCluster} replace the previous function name \code{fitINLA2} (before version 1.0.0).
#' 
#' @param data count data of person-months with the following columns
#' \itemize{
#'  \item cluster: cluster ID
#'  \item years: time period
#'  \item region: region of the cluster
#' \item strata: stratum of the cluster
#' \item age: age group corresponding to the row
#' \item total: total number of person-month in this age group, stratum, cluster, and period
#' \item Y: total number of deaths in this age group, stratum, cluster, and period
#' }
#' @param X Covariate matrix. It must contain either a column with name "region", or a column with name "years", or both. The covariates must not have missing values for all regions (if varying in space) and all time periods (if varying in time). The rest of the columns are treated as covariates in the mean model.
#' @param age.groups a character vector of age groups in increasing order.
#' @param age.n number of months in each age groups in the same order.
#' @param age.rw.group vector indicating grouping of the ages groups. For example, if each age group is assigned a different random walk component, then set age.rw.group to c(1:length(age.groups)); if all age groups share the same random walk component, then set age.rw.group to a rep(1, length(age.groups)). The default for 6 age groups is c(1,2,3,3,3,3), which assigns a separate random walk to the first two groups and a common random walk for the rest of the age groups. The vector should contain values starting from 1.
#' @param age.strata.fixed.group vector indicating grouping of the ages groups for different strata. The default is c(1:length(age.groups)), which correspond to each age group within each stratum receives a separate intercept. If several age groups are specified to be the same value in this vector, the stratum specific deviation from the baseline is assumed to be the same for these age groups. For example, if \code{age.strata.fixed.group = c(1, 2, 3, 3, 3, 3)}, then the fixed effect part of the linear predictor consists of 6 overall age-specific intercepts and 3 set of strata effects (where a base stratum is chosen internally), for age groups 1, 2, and the rest respectively.
#' 
#' For example, if each age group is assigned a different intercept, then set age.strata.fixed.group to c(1:length(age.groups)); if all age groups share the same intercept, then set age.strata.fixed.group to a rep(1, length(age.groups)). The default for 6 age groups is the former. It can also be set to be the same as \code{age.rw.group}. The vector should contain values starting from 1.
#' @param family family of the model. This can be either binomial (with logistic normal prior), betabiniomial.
#' @param time.model Model for the main temporal trend, can be rw1, rw2, ar1, or NULL (for spatial-only smoothing). Default to be rw2. For ar1 main effect, a linear slope is also added with time scaled to be between -0.5 to 0.5, i.e., the slope coefficient represents the total change between the first year and the last year in the projection period on the logit scale. 
#' @param st.time.model Temporal component model for the interaction term, can be rw1, rw2, or ar1. Default to be the same as time.model unless specified otherwise. The default does not include region-specific random slopes. They can be added to the interaction term by specifying \code{pc.st.slope.u} and \code{pc.st.slope.alpha}.  
#' @param Amat Adjacency matrix for the regions
#' @param bias.adj the ratio of unadjusted mortality rates or age-group-specific hazards to the true rates or hazards. It needs to be a data frame that can be merged to thee outcome, i.e., with the same column names for time periods (for national adjustment), or time periods and region (for subnational adjustment). The column specifying the adjustment ratio should be named "ratio".
#' @param bias.adj.by vector of the column names specifying how to merge the bias adjustment to the count data. For example, if bias adjustment factor is provided in bias.adj for each region and time, then bias.adj.by should be `c("region", "time")`.
#' @param formula INLA formula.  See vignette for example of using customized formula.
#' @param year_label string vector of year names
#' @param type.st type for space-time interaction
#' @param survey.effect logical indicator whether to include a survey fixed effect. If this is set to TRUE, there needs to be a column named 'survey' in the input data frame. In prediction, this effect term will be set to 0.
#' @param strata.time.effect logical indicator whether to include strata specific temporal trends.  
#' @param linear.trend logical indicator whether a linear trend is added to the temporal main effect. If the temporal main effect is RW2, then it will be forced to FALSE. Default is TRUE.
#' @param common.trend logical indicator whether all age groups and/or strata share the same linear trend in the temporal main effect.  
#' @param hyper Deprecated. which hyperpriors to use. Only supports PC prior ("pc"). 
#' @param pc.u hyperparameter U for the PC prior on precisions.
#' @param pc.alpha hyperparameter alpha for the PC prior on precisions.
#' @param pc.u.phi hyperparameter U for the PC prior on the mixture probability phi in BYM2 model.
#' @param pc.alpha.phi hyperparameter alpha for the PC prior on the mixture probability phi in BYM2 model.
#' @param pc.u.cor hyperparameter U for the PC prior on the autocorrelation parameter in the AR prior, i.e. Prob(cor > pc.u.cor) = pc.alpha.cor.
#' @param pc.alpha.cor hyperparameter alpha for the PC prior on the autocorrelation parameter in the AR prior.
#' @param pc.st.u hyperparameter U for the PC prior on precisions for the interaction term.
#' @param pc.st.alpha hyperparameter alpha for the PC prior on precisions for the interaction term.
#' @param pc.st.slope.u hyperparameter U for the PC prior on precisions for the area-level random slope. If both pc.st.slope.u and pc.st.slope.alpha are not NA, an area-level random slope with iid prior will be added to the model. The parameterization of the random slope is so that Prob(|beta| > pc.st.slope.u) = pc.st.slope.alpha, where time covariate is rescaled to be -0.5 to 0.5, so that the random slope can be interpreted as the total deviation from the main trend from the first year to the last year to be projected, on the logit scale. 
#' @param pc.st.slope.alpha hyperparameter alpha for the PC prior on precisions for the area-level random slope. See above for the parameterization.
#' @param overdisp.mean hyperparameter for the betabinomial likelihood. Mean of the over-dispersion parameter on the logit scale. 
#' @param overdisp.prec hyperparameter for the betabinomial likelihood. Precision of the over-dispersion parameter on the logit scale. 
#' @param options list of options to be passed to control.compute() in the inla() function.
#' @param control.inla list of options to be passed to control.inla() in the inla() function. Default to the "adaptive" integration strategy.
#' @param verbose logical indicator to print out detailed inla() intermediate steps.
#' @param rw Deprecated. Take values 0, 1 or 2, indicating the order of random walk. If rw = 0, the autoregressive process is used instead of the random walk in the main trend. See the description of the argument ar for details.
#' @param ar Deprecated. Order of the autoregressive component. If ar is specified to be positive integer, the random walk components will be replaced by AR(p) terms in the interaction part. The main temporal trend remains to be random walk of order rw unless rw = 0.
#' @param st.rw Deprecated. Take values 1 or 2, indicating the order of random walk for the interaction term. If not specified, it will take the same order as the argument rw in the main effect. Notice that this argument is only used if ar is set to 0.
#' @param geo Deprecated. Spatial polygon file, legacy parameter from previous versions of the package.
#' @param ... arguments to be passed to the inla() function call.
#' @seealso \code{\link{getDirect}}
#' @import Matrix
#' @importFrom stats dgamma na.pass
#' @importFrom Matrix Diagonal 
#' @return INLA model fit using the provided formula, country summary data, and geographic data
#' @author Zehang Richard Li 
#' @examples
#' \dontrun{
#' library(dplyr)
#' data(DemoData)
#' # Create dataset of counts
#' counts.all <- NULL
#' for(i in 1:length(DemoData)){
#'   counts <- getCounts(DemoData[[i]][, c("clustid", "time", "age", "died",
#'                                         "region", "strata")],
#'             variables = 'died', by = c("age", "clustid", "region", 
#'                                          "time", "strata"))
#'   counts <- counts %>% mutate(cluster = clustid, years = time, Y=died)
#'   counts$strata <- gsub(".*\\.","",counts$strata)
#'   counts$survey <- names(DemoData)[i] 
#'   counts.all <- rbind(counts.all, counts)
#' }
#' 
#' # fit cluster-level model on the periods
#' periods <- levels(DemoData[[1]]$time)
#' fit <- smoothCluster(data = counts.all, 
#'      Amat = DemoMap$Amat, 
#'      time.model = "rw2", 
#'      st.time.model = "rw1",
#'      strata.time.effect =  TRUE, 
#'      survey.effect = TRUE,
#'      family = "betabinomial",
#'      year_label = c(periods, "15-19"))
#' summary(fit)
#' est <- getSmoothed(fit, nsim = 1000)
#' plot(est$stratified, plot.CI=TRUE) + ggplot2::facet_wrap(~strata) 
#' 
#' # fit cluster-level space-time model with covariate
#' # notice without projected covariates, we use periods up to 10-14 only
#' # construct a random covariate matrix for illustration
#' periods <- levels(DemoData[[1]]$time)
#' X <- expand.grid(years = periods, 
#'        region = unique(counts.all$region))
#' X$X1 <- rnorm(dim(X)[1])
#' X$X2 <- rnorm(dim(X)[1])
#' fit.covariate <- smoothCluster(data = counts.all, 
#'    X = X,
#'      Amat = DemoMap$Amat, 
#'      time.model = "rw2", 
#'      st.time.model = "rw1",
#'      strata.time.effect =  TRUE, 
#'      survey.effect = TRUE,
#'      family = "betabinomial",
#'      year_label = c(periods))
#' est <- getSmoothed(fit.covariate, nsim = 1000)
#' 
#' # fit cluster-level model for one time point only
#' # i.e., space-only model
#' fit.sp <- smoothCluster(data = subset(counts.all, time == "10-14"), 
#'      Amat = DemoMap$Amat, 
#'      time.model = NULL, 
#'      survey.effect = TRUE,
#'      family = "betabinomial")
#' summary(fit.sp)
#' est <- getSmoothed(fit.sp, nsim = 1000)
#' plot(est$stratified, plot.CI = TRUE) + ggplot2::facet_wrap(~strata) 
#' 
#' # fit cluster-level model for one time point and covariate
#' # construct a random covariate matrix for illustration
#' X <- data.frame(region = unique(counts.all$region),
#'       X1 = c(1, 2, 2, 1), 
#'       X2 = c(1, 1, 1, 2))
#' fit.sp.covariate <- smoothCluster(data = subset(counts.all, time == "10-14"), 
#'      X = X, 
#'      Amat = DemoMap$Amat, 
#'      time.model = NULL, 
#'      survey.effect = TRUE,
#'      family = "betabinomial")
#' summary(fit.sp.covariate)
#' est <- getSmoothed(fit.sp.covariate, nsim = 1000)
#' }
#' 
#' @export
#' 
#' 

smoothCluster <- function(data, X = NULL, family = c("betabinomial", "binomial")[1], age.groups = c("0", "1-11", "12-23", "24-35", "36-47", "48-59"), age.n = c(1,11,12,12,12,12), age.rw.group = c(1,2,3,3,3,3), age.strata.fixed.group = c(1,2,3,4,5,6), time.model = c("rw1", "rw2", "ar1")[2], st.time.model = NULL, Amat, bias.adj = NULL, bias.adj.by = NULL, formula = NULL, year_label, type.st = 4, survey.effect = FALSE, linear.trend = TRUE, common.trend = FALSE, strata.time.effect = FALSE, hyper = "pc", pc.u = 1, pc.alpha = 0.01, pc.u.phi = 0.5, pc.alpha.phi = 2/3, pc.u.cor = 0.7, pc.alpha.cor = 0.9,  pc.st.u = NA, pc.st.alpha = NA, pc.st.slope.u = NA, pc.st.slope.alpha = NA, overdisp.mean = 0, overdisp.prec = 0.4, options = list(config = TRUE), control.inla = list(strategy = "adaptive", int.strategy = "auto"), verbose = FALSE, geo = NULL, rw = NULL, ar = NULL, st.rw = NULL, ...){

  # if(family == "betabinomialna") stop("family = betabinomialna is still experimental.")
  # check region names in Amat is consistent

  # add check: if user input options and forgot config, set it to TRUE
  #            if user turned if explicitly, then leave it along

  msg <- NULL

  if(!"config" %in% names(options)){
    message("config = TRUE is added to options so that posterior draws can be taken.")
    options$config <- TRUE
  }
  if(!is.null(geo)){
    message("Argument geo is deprecated in the smoothCluster function. Only Amat is needed.")
  }

  # get around CRAN check of using un-exported INLA functions
  rate0 <- shape0 <- my.cache <- inla.as.sparse <- type <- strata <- rescale_U <- sim_alpha <- pc.st.slope.prec.u <- NULL

  if (!isTRUE(requireNamespace("INLA", quietly = TRUE))) {
    stop("You need to install the packages 'INLA'. Please run in your R terminal:\n  install.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/stable'), dep=TRUE)")
  }
  if (!is.element("Matrix", (.packages()))) {
    attachNamespace("Matrix")
  }
  # If INLA is installed, then attach the Namespace (so that all the relevant functions are available)
  if (isTRUE(requireNamespace("INLA", quietly = TRUE))) {
    if (!is.element("INLA", (.packages()))) {
      attachNamespace("INLA")
    }

  # Backward compatibility
  if(hyper != "pc"){
    hyper <- "pc"
    message("Argument 'hyper' have been deprecated in version 1.1.0. The hyper parameters will be specified with PC priors.")
  }
  # define whether there is the temporal component
  is.temporal <- TRUE

  if(!is.null(rw) || !is.null(ar) || !is.null(st.rw)){

    if(is.null(ar)) ar <- 0
    if(is.null(st.rw)) st.rw <- rw
    st.ar <- ar
    is.ar <- ar > 0 
    is.main.ar <- rw == 0

    # if(is.ar) message("ar > 0: using the AR(p) process for the space-time interaction component.")
    if(rw %in% c(0, 1, 2) == FALSE) stop("Random walk only support rw = 1 or 2.")
    # if(rw == 0) message("rw = 0: using the AR(p) process for the main temporal trend component.")
    # if(!is.ar  && st.rw != rw) message("The main effect is random walk of order ", rw, ", and the interaction effects are random walks of order ", st.rw)

    time.model <- paste0("rw", rw)
    if(is.main.ar) time.model <- paste0("ar", ar)
    st.time.model <- paste0("rw", st.rw)
    if(is.ar) st.time.model <- paste0("ar", ar)  
    message(paste0("Argument 'rw' and 'ar' have been deprecated in version 1.0.0. They still work as intended for now. In the future, use time.model and st.time.model to specify temporal components\n e.g., time.model = '", time.model, "', st.time.model = '", st.time.model, "'"))

  # New time mode definition...    
  }else{
    if(any(is.null(st.time.model)) || any(is.na(st.time.model))){
      st.time.model <- time.model
    }
    if(any(is.na(time.model)) || any(is.null(time.model))){
      time.model <- NULL
    }
    if(!is.null(time.model)){
      time.model <- tolower(time.model)
    }else{
      is.temporal <- FALSE
      strata.time.effect <- FALSE
      year_label <- NULL
    }

    st.time.model <- tolower(st.time.model)

    is.main.ar = FALSE
    is.ar  = FALSE
    ar = 0
    st.ar = 0

    # main effect
    if(is.temporal){
      if(time.model == "rw1"){
        rw = 1
      }else if(time.model == "rw2"){
        rw = 2
      }else if(time.model == "ar1"){ # for future: allow AR(p) here
        rw = 1 # replaced later
        ar = 1
        is.main.ar = TRUE
      }else{
        stop("Temporal main effect only support rw1, rw2, ar1, or NULL.")
      }      
    }else{
      rw = 0
      ar = 0
    }

    # interaction effect
    if(!is.temporal){
      st.rw = 0
    }else{
      if(st.time.model == "rw1"){
        st.rw = 1
      }else if(st.time.model == "rw2"){
        st.rw = 2
      }else if(st.time.model == "ar1"){ # for future: allow AR(p) here
        st.rw = NULL
        st.ar = 1
        is.ar = TRUE
      }else{
        stop("Temporal interaction effect only support rw1, rw2, and ar1.")
      }
    }
  }
  if(time.model == "rw2"){
    linear.trend = FALSE
    common.trend = FALSE
  }
  if(!is.temporal){
    message("----------------------------------",
            "\nCluster-level model",
            "\n  No temporal components", appendLF = FALSE)
    msg <- paste0(msg, "\nCluster-level model",
                  "\n  No temporal components")

  }else{
    message("----------------------------------",
            "\nCluster-level model",
            "\n  Main temporal model:        ", time.model, appendLF = FALSE)
    msg <- paste0(msg, "\nCluster-level model",
                       "\n  Main temporal model:        ", time.model)
    if(linear.trend && !common.trend){
        message("\n  Additional linear trends:   stratum-age-specific", appendLF = FALSE)
        msg <- paste0(msg, "\n  Additional linear trends:   stratum-age-specific")
    }
    if(linear.trend && common.trend){
        message("\n  Additional linear trends:   shared", appendLF = FALSE)
        msg <- paste0(msg, "\n  Additional linear trends:   shared")
    }
  }

  
  if(!is.null(Amat)){
    if(is.null(rownames(Amat))){
        stop("Row names of Amat needs to be specified to region names.")
    }
    if(is.null(colnames(Amat))){
        stop("Column names of Amat needs to be specified to region names.")
    }
    if(sum(rownames(Amat) != colnames(Amat)) > 0){
        stop("Row and column names of Amat needs to be the same.")
    }
    is.spatial <- TRUE
    if(sum(!data$region %in% colnames(Amat)) > 0){
        stop("Exist regions in the data frame but not in Amat.")
    }
  }else{
    Amat <- matrix(1,1,1)
    colnames(Amat) <- rownames(Amat) <- "All"
    data$region <- "All"
    is.spatial <- FALSE
    if(!is.temporal){
      stop("Neither spatial or temporal component is specified.")
    }
  }

  if(is.spatial){ 
    message("\n  Spatial effect:             bym2", appendLF=FALSE) 
    msg <- paste0(msg, "\n  Spatial effect:             bym2")
  }
  if(is.spatial && is.temporal){
    message("\n  Interaction temporal model: ", st.time.model, 
            "\n  Interaction type:           ", type.st, appendLF=FALSE)
    msg <- paste0(msg, "\n  Interaction temporal model: ", st.time.model, 
            "\n  Interaction type:           ", type.st)

  }
  if(!is.na(pc.st.slope.u) && !is.na(pc.st.slope.alpha)){
    message("\n  Interaction random slopes:  yes", appendLF=FALSE)
    msg <- paste0(msg, "\n  Interaction random slopes:  yes")
  }else{
    message("\n  Interaction random slopes:  no", appendLF=FALSE)
    msg <- paste0(msg, "\n  Interaction random slopes:  no")
  }   


  if(is.null(age.groups)){
    age.n <- 1
    age.rw.group <- 1
  }
  age.strata.fixed.group <- age.strata.fixed.group[1:length(age.n)]

  message("\n  Number of age groups: ", length(age.n), appendLF=FALSE)
  msg <- paste0(msg, "\n  Number of age groups: ", length(age.n))


  if(is.ar && hyper=="gamma"){
    stop("AR1 model only implemented with PC priors for now.")
  }
  


  multi.frame <- FALSE
  if("frame" %in% colnames(data)){
    message("\n  Column specifying sampling frame(s): yes", appendLF=FALSE)
    msg <- paste0(msg, "\n  Column specifying sampling frame(s): yes")

    if(length(unique(data$frame)) > 1){
      multi.frame <- TRUE
      if("strata" %in% colnames(data) == FALSE || sum(!is.na(data$strata)) == 0){
        data$strata <- data$frame 
        message("\n  Strata renaming: yes, renamed to the same as frame", appendLF=FALSE)
        msg <- paste0(msg, "\n  Strata renaming: yes, renamed to the same as frame")

      }else{
        data$strata <- paste(data$frame, data$strata, sep = "-") 
        message("\n  Strata renaming: yes, renamed frame-strata", appendLF=FALSE)
        msg <- paste0(msg, "\n  Strata renaming: yes, renamed frame-strata")
     }
      
    }else{
      message("\n  Strata renaming: no, only one frame", appendLF=FALSE)
      msg <- paste0(msg, "\n  Strata renaming: no, only one frame")
    }
  }


  has.strata <- TRUE
  if("strata" %in% colnames(data) == FALSE || sum(!is.na(data$strata)) == 0){
    has.strata <- FALSE
    data$strata <- ""
    strata.time.effect <- FALSE
    message("\n  Stratification: no, strata variable not in the input", appendLF=FALSE)
    msg <- paste0(msg, "\n  Stratification: no, strata variable not in the input")

  }else if(length(unique(data$strata[data$total > 0])) == 1){
    # if only one strata has observations
    has.strata <- FALSE
    message("\n  Stratification: no, all data in the same stratum", appendLF=FALSE)
    msg <- paste0(msg, "\n  Stratification: no, all data in the same stratum")

  }else if(length(unique(data$strata[data$total > 0])) < length(unique(data$strata))){
    # If there are more than 2 strata and we need to drop levels
    data <- subset(data, strata %in% unique(data$strata[data$total > 0])) 
    message("\n  Stratification: yes, but some frame-strata combination does not exist", appendLF=FALSE)
    msg <- paste0(msg, "\n  Stratification: yes, but some frame-strata combination does not exist")

  }else{
    message("\n  Stratification: yes", appendLF=FALSE)
    msg <- paste0(msg, "\n  Stratification: yes")
  }

  if(has.strata){
    message("\n  Number of age-specific fixed effect per stratum: ", length(unique(age.strata.fixed.group)), appendLF=FALSE)
    msg <- paste0(msg, "\n  Number of age group fixed effect per stratum: ", length(unique(age.strata.fixed.group)))
    message("\n  Number of age-specific trends per stratum: ", length(unique(age.rw.group)), appendLF=FALSE)
    msg <- paste0(msg, "\n  Number of age-specific trends per stratum: ", length(unique(age.rw.group)))    
  }


  ## Survey fixed effect should only exist if they are nested within the same strata.
  ## This is checked here
  survey.message <- NULL
  if(survey.effect){
      if("survey" %in% colnames(data) == FALSE){
          survey.message <- "\n  Survey effect: no, survey column does not exist"
          survey.effect <- FALSE
          survey.table <- NULL
      }else if(!multi.frame){
        # If there is no frame variable in the input data
          data$survey2 <- data$survey
          survey_count <- length(unique(data$survey2))
          survey.table <- data.frame(survey = unique(data$survey2), 
                                     survey.id = 1:survey_count)
          data$survey.id <- match(data$survey2, survey.table$survey)
          survey.A <- matrix(1, 1, survey_count)
          survey.e <- 0
          if(survey_count == 1){
            survey.message <- "\n  Survey effect: no, only one survey"
            survey.effect <- FALSE
          }

      }else{
          # If there are frames, may need to remove survey effects
          # Then impose sum-to-zero constraints for each frame
          tab <- table(data$survey, data$frame)
          if(max(apply(tab, 1, function(x){sum(x>0)})) > 1){
            stop("There exist surveys associated with more than one frames. Please check the 'survey' column and 'frame' column of the input data.")
          }
          # for each frame, check if multiple survey exist
          nested <- apply(tab, 2, function(x){sum(x!=0) > 1})
          if(sum(nested) == 0){
            survey.message <- "\n  Survey effect: no, no multiple surveys nested within frames"
            survey.effect <- FALSE
            survey.table <- NULL
          }else{
            tt <- NULL
            for(j in which(nested)) tt <- c(tt, which(tab[, j] > 0)) 
            tt <- rownames(tab)[sort(unique(tt))]
            data$survey2 <- as.character(data$survey)
            data$survey2[data$survey2 %in% tt == FALSE] <- NA
            ss <- unique(data$survey2)
            ss <- ss[!is.na(ss)]
            survey_count <- length(ss)
            # Obtain the correct sum to zero constraints within Frame
            survey.table <- data.frame(survey = ss, 
                                       survey.id = 1:survey_count)
            survey.table$frame <- data[match(survey.table$survey, data$survey), "frame"]
            nested.frame <- unique(survey.table$frame)
            survey.A <- matrix(0, length(nested.frame), survey_count)
            for(j in 1:dim(survey.A)[1]){
                same_frame <- which(survey.table$frame == nested.frame[j])
                survey.A[j, same_frame] <- 1
            }
            survey.e <- rep(0, length(nested.frame))
            data$survey.id <-match(data$survey2, survey.table$survey)
          }
      }
  }else{
    survey.table <- NULL
  }
  if(survey.effect) survey.message <- "\n  Survey effect: yes"
  stratalevels <- unique(data$strata)


age.groups.vec <- age.groups
data$age.intercept <- data$age
data$age.orig <- data$age

# If only one strata, this becomes the overall trend
if(strata.time.effect){
    if(is.null(data$strata)) stop("No strata column in the data.")
    ## Update here age -> age x strata to allow stratum-specific trends
    ## Then redefine age.n and age.rw.group.
    strata.n <- length(stratalevels) 

    age.n <- rep(age.n, strata.n)

    tmp <- NULL
    ngroup <- max(age.rw.group)
    for(i in 1:strata.n){
        tmp <- c(tmp, age.rw.group + (i-1) * ngroup)
    }
    age.rw.group <- tmp

    if(is.null(age.groups)){
      age.groups <- stratalevels
      data$age <- data$strata
      data$age.intercept <- data$strata
    }else{
      age.groups <- expand.grid(age.groups.vec, stratalevels)
      age.groups <- paste(age.groups[,1], age.groups[,2], sep = ":")
      data$age <- paste(data$age, data$strata, sep = ":")
      # deal with age.intercept in the next section
    }
  }

  # make new age.intercept variable
  is.age.diff <- FALSE
  age.diff.levels <- NULL
  if(!is.null(age.groups)){
      # reformat age variable to correspond to intercepts
      age.groups.new <- NULL
      for(k in 1:max(age.strata.fixed.group)){
          age.groups.new <- c(age.groups.new, paste0(age.groups.vec[which(age.strata.fixed.group == k)], collapse = ", "))
      }
      # make the list the same length as the original age.groups
      age.groups.new <- age.groups.new[age.strata.fixed.group]

      # age-stratum effect modeled as difference from baseline?
      is.age.diff <- ifelse(length(unique(age.strata.fixed.group)) == length(age.strata.fixed.group), 0, 1)
      if(is.age.diff){
        # new age.intercept variable as just original age
        data$age.intercept <- data$age.orig
        # add new age.diff variable as the difference 
        data$age.diff <- age.groups.new[match(data$age.orig, age.groups.vec)]
        data$age.diff <- paste(data$age.diff, data$strata, sep = ":") 
        ## 
        ##  This artificial first level (instead of NA) is needed to 
        ##    make inla specify the model.matrix correctly.
        ##
        data$age.diff[data$strata == stratalevels[1]] <- "base"
        age.diff.levels <- sort(unique(data$age.diff))
        age.diff.levels <- c("base", age.diff.levels[age.diff.levels != "base"])
      }else{          
        # new age.intercept variable as age * strata
        data$age.intercept <- age.groups.new[match(data$age.orig, age.groups.vec)]
        data$age.intercept <- paste(data$age.intercept, data$strata, sep = ":") 
        data$age.diff <- NA

        age.groups.new <- expand.grid(age.groups.new, stratalevels)
        age.groups.new <- paste(age.groups.new[,1], age.groups.new[,2], sep = ":")
      }
  }
  if(strata.time.effect){
    message("\n  Strata-specific temporal trends: yes", appendLF=FALSE)
    msg <- paste0(msg, "\n  Strata-specific temporal trends: yes")
  }
  if(!is.null(survey.message)){
    message(survey.message, appendLF=FALSE)
    msg <- paste0(msg, survey.message)
  }

    
    ## ---------------------------------------------------------
    ## Common Setup
    ## --------------------------------------------------------- 
    if(dim(Amat)[1] != 1){
      data <- data[which(data$region != "All"), ]
    }  
    ####  Re-calculate hyper-priors (for gamma only)     
    # priors <- simhyper(R = 2, nsamp = 1e+05, nsamp.check = 5000, only.iid = TRUE)
    # a.iid <- priors$a.iid
    # b.iid <- priors$b.iid
    # a.rw <- priors$a.iid
    # b.rw <- priors$b.iid
    # a.icar <- priors$a.iid
    # b.icar <- priors$b.iid
    

    #### get the list of region and numeric index in one data frame
    if(!is.spatial){
      region_names <- regions <- "All"
      region_count <- S <- 1
      dat <- cbind(data, region_number = 0)
    }else{
      region_names <- colnames(Amat) 
      region_count <- S <- length(region_names)
      regions <- data.frame(region = region_names, region_number = seq(1, region_count))      
      # -- merging in the alphabetical region number -- #
      dat <- merge(data, regions, by = "region")
    }
    # -- creating IDs for the spatial REs -- #
    dat$region.struct <- dat$region.unstruct <- dat$region.int <- dat$region_number
    
    #### get the list of region and numeric index in one data frame    
    n <- 0
    N <- length(year_label)
    if(N > 0){
      years <- data.frame(year = year_label, year_number = seq(1, N))   
      # -- creating IDs for the temporal REs -- #
      dat$time.unstruct <- dat$time.struct <- dat$time.int <- years[match(dat$years, years[, 1]), 2]
      if(sum(is.na(dat$time.unstruct)) > 0){
        dat <- dat[!is.na(dat$time.unstruct), ]
        message("Data contains time periods not in the specified range. These observations are removed.")
      }
      x <- expand.grid(1:N, 1:region_count)
      time.area <- data.frame(region_number = x[, 2], time.unstruct = x[, 1], time.area = c(1:nrow(x)))
      # fix for 0 instead of 1 when no Amat file provided
      if(!is.spatial){
        time.area$region_number <- 0
      }
      # -- merge these all into the data sets -- #
      newdata <- dat
      if(dim(Amat)[1] != 1){
        newdata <- merge(newdata, time.area, 
          by = c("region_number", "time.unstruct"))
      }else{
        newdata$time.area <- NA
      }
    }else{
      years <- NULL
      dat$time.unstruct <- dat$time.struct <- dat$time.int <- dat$time.area <- NA
      newdata <- dat
      time.area <- NULL
    }
   
    if(!is.null(X)){
      by <- NULL
      if("region" %in% colnames(X)){
        by <- c(by, "region")
      }
      if("years" %in% colnames(X)){
        by <- c(by, "years")
      }
      covariate.names <- colnames(X)[colnames(X) %in% by == FALSE]
      if(length(covariate.names) == 0){
        warning("The X argument is specified but no covariate identified.")
        X <- NULL
      }
      if("region" %in% by && (!"years" %in% by)){
        for(ii in unique(newdata$region)){
          which <- which(X$region == ii)
          if(length(which) == 0){
            stop(paste("Missing region in the covariate matrix:", ii))
          }
          if(length(which) > 1) stop(paste("Duplicated covariates for region", ii))
          if(sum(is.na(X[which, ])) > 0) stop("NA in covariate matrix.")
        }
      }else if((!"region" %in% by) && "years" %in% by){
        for(tt in unique(newdata$years)){
          which <- which(X$years == tt)
          if(length(which) == 0){
            stop(paste("Missing years in the covariate matrix:", tt))
          }
          if(length(which) > 1) stop(paste("Duplicated covariates for period", tt))
          if(sum(is.na(X[which, ])) > 0) stop("NA in covariate matrix.")
        }
      }else if("region" %in% by && "years" %in% by){
        for(tt in unique(newdata$years)){
          for(ii in unique(newdata$region)){
            which <- intersect(which(X$years == tt), which(X$region == ii))
            if(length(which) == 0){
              stop(paste("Missing region-years in the covariate matrix:", ii, tt))
            }
           if(length(which) > 1) stop(paste("Duplicated covariates for region-years", ii, tt))
           if(sum(is.na(X[which, ])) > 0) stop("NA in covariate matrix.")
          }
        }
      }else{
        stop("Covariate need to contain column 'region' or 'years'.")
      }
      # check if prediction year exist in covariates
      if("years" %in% by){
        for(tt in unique(year_label)){
          which <- which(X$years == tt)
          if(length(which) == 0){
            stop(paste("Missing years in the covariate matrix:", tt))
          }
        }
      }
    }else{
      covariate.names <- NULL
    }

    if(!is.null(X)){
      message(paste0("\n  Covariates: ", paste(covariate.names, collapse = ", ")), appendLF=FALSE)
      message(paste0("\n  Covariates by: ", paste(by, collapse = ", ")), appendLF=FALSE)
      msg <- paste0(msg, paste0("\n  Covariates: ", paste(covariate.names, collapse = ", ")), paste0("\n  Covariates by: ", paste(by, collapse = ", ")))
    }

    message("\n----------------------------------")
    msg <- paste0(msg, "\n")
    
    #### format additional columns in the input
    exdat <- newdata
    if("survey" %in% colnames(exdat)){
      exdat$cluster <- paste(exdat$survey, exdat$cluster)
    }
    clusters <- unique(exdat$cluster)
    exdat$cluster.id <- match(exdat$cluster, clusters)
    exdat$nugget.id <- exdat$cluster.id
    # cluster.time <- expand.grid(cluster = clusters, time = 1:N)
    # cluster.time$nugget.id <- 1:dim(cluster.time)[1]
    # exdat <- merge(exdat, cluster.time, by.x = c("cluster", "time.struct"), by.y = c("cluster", "time"))
    # # exdat$nugget.id <- 1:dim(exdat)[1]

    if(is.ar || is.main.ar || linear.trend){
      exdat$time.slope <- exdat$time.struct 
      center <- N/2 + 1e-5 # avoid exact zero in the lincomb creation
      exdat$time.slope <- (exdat$time.slope  - center) / ((N - 1))
      slope.fixed.names <- NULL

      for(aa in unique(age.rw.group)){
          tmp <- paste0("time.slope.group", aa)
          exdat[, tmp] <- exdat$time.slope
          which.age <- age.groups[which(age.rw.group == aa)]
          this.group <- exdat$age %in% which.age
          exdat[!this.group, tmp] <- NA
          slope.fixed.names <- c(slope.fixed.names, tmp)
      }
    }
    replicate.rw <- length(unique(age.rw.group)) > 1

  if(is.null(formula)){
   
    ## ---------------------------------------------------------
    ## Setup PC prior model
    ## ---------------------------------------------------------
    # if(tolower(hyper) == "pc"){
        hyperpc1 <- list(prec = list(prior = "pc.prec", param = c(pc.u , pc.alpha)))
        hyperpc2 <- list(prec = list(prior = "pc.prec", param = c(pc.u , pc.alpha)), 
                         phi = list(prior = 'pc', param = c(pc.u.phi , pc.alpha.phi)))
        hyperar1 = list(prec = list(prior = "pc.prec", param = c(pc.u , pc.alpha)), 
                        theta2 = list(prior = "pc.cor1", param = c(pc.u.cor, pc.alpha.cor)))
        hyperar2 = list(theta2 = list(prior = "pc.cor1", param = c(pc.u.cor, pc.alpha.cor)))
        pc.st.u <- ifelse(is.na(pc.st.u), pc.u, pc.st.u)
        pc.st.alpha <- ifelse(is.na(pc.st.alpha), pc.alpha, pc.st.alpha)
        hyperpc1.interact <- list(prec = list(prior = "pc.prec", param = c(pc.st.u , pc.st.alpha)))
        ## -----------------------
        ##  National + PC
        ## ----------------------- 
        if(!is.spatial && is.temporal){

          if(replicate.rw && (!is.main.ar)){
              # Replicated RW
              formula <- Y ~ f(time.struct,model=paste0("rw", rw), constr = TRUE,  extraconstr = NULL, hyper = hyperpc1, replicate =  age.rep.idx)  
          }else if(!is.main.ar){
              # Single RW
              formula <- Y ~ f(time.struct,model=paste0("rw", rw), constr = TRUE,  extraconstr = NULL, hyper = hyperpc1)          
          }else if(replicate.rw && is.main.ar){
            # Replicated AR1
            formula <- Y ~  f(time.struct, model="ar", order=ar, hyper = hyperar1, constr = TRUE, extraconstr = NULL, values = 1:N, replicate =  age.rep.idx) 
          }else{
            # Single AR1
           formula <- Y ~  f(time.struct, model="ar", order=ar, hyper = hyperar1, constr = TRUE, extraconstr = NULL, values = 1:N) 

          }
         
          # AR1 add slope
           if(linear.trend && !common.trend){
              tmp <- paste(slope.fixed.names, collapse = " + ")
              formula <- as.formula(paste(c(formula, tmp), collapse = "+"))
           }
           if(linear.trend && common.trend){
              formula <- as.formula(paste(c(formula, "time.slope"), collapse = "+"))
           }

          # Add IID
          formula <- update(formula, ~. + 
                  f(time.unstruct,model="iid", hyper = hyperpc1, values = 1:N))

        ## -------------------------
        ## Subnational + PC
        ## ------------------------- 
        }else if(dim(Amat)[1] != 1 && is.temporal){
             if(replicate.rw && (!is.main.ar)){
              # Replicated RW
              formula <- Y ~ 
                f(time.struct, model=paste0("rw", rw), hyper = hyperpc1, scale.model = TRUE, extraconstr = NULL, values = 1:N, replicate =  age.rep.idx) 

             }else if(!is.main.ar){
              # Single RW
               formula <- Y ~ 
                f(time.struct, model=paste0("rw", rw), hyper = hyperpc1, scale.model = TRUE, extraconstr = NULL, values = 1:N) 
             }else if(replicate.rw && is.main.ar){
              # Replicated AR1
               formula <- Y ~ 
                f(time.struct, model="ar", order=ar, hyper = hyperar1, constr = TRUE, extraconstr = NULL, values = 1:N, replicate =  age.rep.idx) 
             }else{
              # Single AR1
              formula <- Y ~ 
                f(time.struct, model="ar", order=ar, hyper = hyperar1, constr = TRUE, extraconstr = NULL, values = 1:N) 
             }

            # AR1 add slope
             if(linear.trend && !common.trend){
                tmp <- paste(slope.fixed.names, collapse = " + ")
                formula <- as.formula(paste(c(formula, tmp), collapse = "+"))
             }
             if(linear.trend && common.trend){
                formula <- as.formula(paste(c(formula, "time.slope"), collapse = "+"))
             }


             formula <- update(formula, ~. + 
                  f(time.unstruct,model="iid", hyper = hyperpc1, values = 1:N) + 
                  f(region.struct, graph=Amat,model="bym2", hyper = hyperpc2, scale.model = TRUE, adjust.for.con.comp = TRUE)) 

            # Interaction models 
            if(type.st == 1){
                formula <- update(formula, ~. + 
                    f(time.area,model="iid", hyper = hyperpc1.interact))
            }else if(type.st == 2){
                if(!is.ar){
                  formula <- update(formula, ~. + 
                    f(region.int,model="iid", group=time.int,control.group=list(model=paste0("rw", st.rw), scale.model = TRUE), hyper = hyperpc1.interact,  adjust.for.con.comp = TRUE))
                }else{
                    formula <- update(formula, ~. + 
                    f(region.int,model="iid", hyper = hyperpc1.interact, group=time.int,control.group=list(model="ar", order = st.ar, hyper = hyperar2), scale.model = TRUE, adjust.for.con.comp = TRUE))
                }
            }else if(type.st == 3){
                formula <- update(formula, ~. + 
                    f(region.int, model="besag", graph = Amat, group=time.int,control.group=list(model="iid"), hyper = hyperpc1.interact, scale.model = TRUE, adjust.for.con.comp = TRUE))
            }else{
                
             # Interaction with AR1  
             if(is.ar){
              formula <- update(formula, ~. + 
                   f(region.int, model="besag", hyper = hyperpc1.interact, graph = Amat, group=time.int,control.group=list(model="ar", order = st.ar, hyper = hyperar2), scale.model = TRUE, adjust.for.con.comp = TRUE)) 
             }else{
                # defines type IV explicitly with constraints
                # Use time.area as index
                # S blocks each with time 1:T (in this code, 1:N)
                 inla.rw = utils::getFromNamespace("inla.rw", "INLA")
                 inla.scale.model.bym = utils::getFromNamespace("inla.scale.model.bym", "INLA")
                 inla.bym.constr.internal = utils::getFromNamespace("inla.bym.constr.internal", "INLA")
                 R2 <- inla.rw(N, order = st.rw, scale.model=TRUE, sparse=TRUE)
                 R4 = Amat
                if(sum(R4 > 0 & R4 < 1) != 0){
                  for(row in 1:nrow(R4)){
                    idx <- which(R4[row,] > 0 & R4[row,] < 1)
                    R4[row,idx] <- 1
                  }
                }
                 diag(R4) <- 0
                 diag <- apply(R4, 1, sum)
                 R4[R4 != 0] <- -1
                 diag(R4) <- diag
                 Q1 <- inla.bym.constr.internal(R4, adjust.for.con.comp = TRUE)
                 R4 <- inla.scale.model.bym(R4, adjust.for.con.comp = TRUE)
                 R <- R4 %x% R2
                 tmp <- matrix(0, S, N * S)
                 for(i in 1:S){
                   tmp[i, ((i-1)*N + 1) : (i*N)] <- 1
                 }
                 # tmp2 <- matrix(0, N, N * S)
                 # for(i in 1:N){
                 #    tmp2[i , which((1:(N*S)) %% N == i-1)] <- 1
                 # }
                 ## Sum-to-zero over connected components
                 tmp2 <- matrix(0, N * Q1$rankdef, N * S)
                 for(j in  1:Q1$rankdef){
                   for(i in 1:N){
                      this_t_i <- which((1:(N*S)) %% N == i-1)
                      this_t_i <- this_t_i[Q1$constr$A[j, ] == 1]
                      tmp2[i + (j-1)*N , this_t_i] <- 1
                   }
                 }
                 tmp <- rbind(tmp[-1, ], tmp2)
                 constr.st <- list(A = tmp, e = rep(0, dim(tmp)[1]))
               
                  formula <- update(formula, ~. + 
                    f(time.area,model="generic0", Cmatrix = R, extraconstr = constr.st, rankdef = N*S -(N - st.rw)*(S - Q1$rankdef), hyper = hyperpc1.interact))
             }
             # END of type IV specification
            }
          
            if(!is.na(pc.st.slope.u) && !is.na(pc.st.slope.alpha)){
          
              # Find u* so that P(|beta| * t.range > u0) = alpha
              #                 where beta ~ Normal(0, tau)
              #                        tau ~ PC(u*, alpha)
              rescale_U <- function(u0, alpha, t.range, tol = 1e-4){

                  sim_alpha <- function(ustar, u0, alpha, t.range){
                        precVal = INLA::inla.pc.rprec(1e6, u = ustar, alpha = alpha)
                        stdDev = 1/sqrt(precVal)
                        betaSim = stats::rnorm(1e6, mean = 0, sd = stdDev)
                        betaSim = abs(betaSim) * t.range 
                        alpha1 <- sum(betaSim > u0) / length(betaSim)
                        return(alpha1)
                  }

                  ustar <- u0 / t.range 
                  ustar.min <- 0
                  ustar.max <- u0 * 2
                  itr = 1
                  while(itr < 10){
                    alpha.max <-  sim_alpha(ustar.max, u0, alpha, t.range)
                    if(alpha.max > alpha){
                      break
                    }else{
                      ustar.max <- ustar.max * 2
                      itr <- itr + 1
                    }
                  }

                  # ss <- al <- seq(ustar.min, ustar.max, length = 100)
                  # for(i in 1:length(ss)) al[i] = sim_alpha(ss[i], u0, alpha, t.range)
                  # plot(ss, al, type = 'l')
                  # abline(v=ustar)

                  # A stupid search...
                  alpha1 = NA
                  while(abs(alpha - alpha1) > tol || is.na(alpha1)){
                        
                        alpha1 <- sim_alpha(ustar, u0, alpha, t.range)

                        if(is.na(alpha1)) stop("Cannot identify the PC prior for the random slope.")

                        if(alpha1 < alpha && abs(alpha - alpha1) > tol){
                           ustar.min <- ustar
                           ustar <- (ustar + ustar.max) / 2
                        }else if(alpha1 > alpha && abs(alpha - alpha1) > tol){
                           ustar.max <- ustar
                           ustar <- (ustar + ustar.min) / 2
                        }
                        itr <- itr + 1
                        if(itr > 1000) stop("Failed to find the PC prior for the random slope.")
                  }

                  return(ustar)
              }

              center <- (N+1)/2 + 1e-5 
              exdat$st.slope <- exdat$time.struct 
              exdat$st.slope <- (exdat$st.slope  - center) / ((N - 1))
              exdat$st.slope.id <- exdat$region.struct

              t.range <- 1
              pc.st.slope.prec.u <- rescale_U(pc.st.slope.u, pc.st.slope.alpha, t.range, tol = 1e-4)
              hyperpc.slope <-  list(prec = list(prior = "pc.prec", param = c(pc.st.slope.prec.u, pc.st.slope.alpha)))
              formula <- update(formula, ~. + f(st.slope.id, st.slope, model = "iid", hyper = hyperpc.slope))
            }
          
        # END of subnational model specification
        # spatial-only model specification
        }else{
          formula <- Y ~ f(region.struct, graph=Amat,model="bym2", hyper = hyperpc2, scale.model = TRUE, adjust.for.con.comp = TRUE) 
        }

        if(survey.effect){
          formula <- update(formula, ~. + f(survey.id, model = "iid", extraconstr = list(A = survey.A, e = survey.e), hyper = list(theta = list(initial=log(0.001), fixed=TRUE))))
        }



    ## ---------------------------------------------------------
    ## Setup Gamma prior model
    ## ---------------------------------------------------------
    # }else if(tolower(hyper) == "gamma"){
    #     ## ------------------- 
    #     ## Period + National
    #     ## ------------------- 
    #     if(!is.spatial){
    #         if(replicate.rw){
    #           formula <- Y ~
    #             f(time.struct,model=paste0("rw", rw),param=c(a.rw,b.rw), constr = TRUE, extraconstr = NULL, hyper = hyperpc1, replicate =  age.rep.idx)  + 
    #             f(time.unstruct,model="iid",param=c(a.iid,b.iid)) 

    #         }else{
    #           formula <- Y ~
    #             f(time.struct,model=paste0("rw", rw),param=c(a.rw,b.rw), constr = TRUE)  + 
    #             f(time.unstruct,model="iid",param=c(a.iid,b.iid))               
    #         }
            
    #     }else{
    #       if(replicate.rw){
    #         formula <- Y ~
    #               f(time.struct,model=paste0("rw", rw), param=c(a.rw,b.rw), scale.model = TRUE, extraconstr = NULL)  
    #       }else{
    #         formula <- Y ~
    #               f(time.struct,model=paste0("rw", rw), param=c(a.rw,b.rw), scale.model = TRUE, extraconstr = NULL)  
    #       }
    #       formula <- update(formula,  ~. +
    #               f(time.unstruct,model="iid",param=c(a.iid,b.iid)) + 
    #               f(region.struct, graph=Amat,model="besag",param=c(a.icar,b.icar), scale.model = TRUE, adjust.for.con.comp = TRUE) + 
    #               f(region.unstruct,model="iid",param=c(a.iid,b.iid))) 
                  
    #         if(type.st == 1){
    #             formula <- update(formula, ~. + f(time.area,model="iid", param=c(a.iid,b.iid)))
    #         }else if(type.st == 2){
    #             formula <- update(formula, ~. + f(region.int,model="iid", group=time.int,control.group=list(model=paste("rw", st.rw), scale.model = TRUE), param=c(a.iid,b.iid)))
    #         }else if(type.st == 3){
    #             formula <- update(formula, ~. + f(region.int,model="besag", graph = Amat, group=time.int,control.group=list(model="iid"),param=c(a.iid,b.iid), scale.model = TRUE, adjust.for.con.comp = TRUE))
    #         }else{
                

    #              # defines type IV explicitly with constraints
    #              # Use time.area as index
    #              # S blocks each with time 1:T (in this code, 1:N)
    #              # UPDATE for connected components:
    #              # nc2 sum-to-zero constraints for each of the connected components of size >= 2. Scaled so that the geometric mean of the marginal variances in each connected component of size >= 2 is 1, and modified so that singletons have a standard Normal distribution.
    #              inla.rw = utils::getFromNamespace("inla.rw", "INLA")
    #              inla.scale.model.bym = utils::getFromNamespace("inla.scale.model.bym", "INLA")
    #              inla.bym.constr.internal = utils::getFromNamespace("inla.bym.constr.internal", "INLA")
    #              R2 <- inla.rw(N, order = st.rw, scale.model=TRUE, sparse=TRUE)
    #              R4 = Amat
    #             if(sum(R4 > 0 & R4 < 1) != 0){
    #               for(row in 1:nrow(R4)){
    #                 idx <- which(R4[row,] > 0 & R4[row,] < 1)
    #                 R4[row,idx] <- 1
    #               }
    #             }
    #              diag(R4) <- 0
    #              diag <- apply(R4, 1, sum)
    #              R4[R4 != 0] <- -1
    #              diag(R4) <- diag
    #              Q1 <- inla.bym.constr.internal(R4, adjust.for.con.comp = TRUE)
    #              R4 <- inla.scale.model.bym(R4, adjust.for.con.comp = TRUE)
    #             R <- R4 %x% R2
    #             tmp <- matrix(0, S, N * S)
    #             for(i in 1:S){
    #               tmp[i, ((i-1)*N + 1) : (i*N)] <- 1
    #             }
    #             # tmp2 <- matrix(0, N, N * S)
    #              # for(i in 1:N){
    #              #    tmp2[i , which((1:(N*S)) %% N == i-1)] <- 1
    #              # }
    #              ## Sum-to-zero over connected components
    #              tmp2 <- matrix(0, N * Q1$rankdef, N * S)
    #              for(j in  1:Q1$rankdef){
    #                for(i in 1:N){
    #                   this_t_i <- which((1:(N*S)) %% N == i-1)
    #                   this_t_i <- this_t_i[Q1$constr$A[j, ] == 1]
    #                   tmp2[i + (j-1)*N , this_t_i] <- 1
    #                }
    #              }
    #             tmp <- rbind(tmp[-1, ], tmp2)
    #             constr.st <- list(A = tmp, e = rep(0, dim(tmp)[1]))
                
    #             if(family == "betabinomialna"){
    #              formula <- update(formula, ~. + 
    #                     f(time.area,model="generic0", Cmatrix = R, extraconstr = constr.st, rankdef = N*S -(N - st.rw)*(S - Q1$rankdef), param=c(a.iid,b.iid), initial=10))
    #              }else{
    #              formula <- update(formula, ~. + 
    #                     f(time.area,model="generic0", Cmatrix = R, extraconstr = constr.st, rankdef = N*S -(N - st.rw)*(S - Q1$rankdef), param=c(a.iid,b.iid)))            
    #              }
    #         }
         
    #     if(survey.effect){
    #       formula <- update(formula, ~. + f(survey.id, model = "iid", extraconstr = list(A = survey.A, e = survey.e), hyper = list(theta = list(initial=log(0.001), fixed=TRUE))))
    #     }
    #     ## ------------------- 
    #     ## Yearly + Subnational
    #     ## ------------------- 
    #     }
    # }else{
    #   stop("hyper needs to be either pc or gamma.")
    # }

    if(family == "binomial"){
      # if(tolower(hyper) == "gamma"){
      #     formula <- update(formula, ~.+ f(nugget.id,model="iid",model="iid", param=c(a.iid,b.iid)))
      # }else if(tolower(hyper) == "pc"){
          formula <- update(formula, ~.+ f(nugget.id,model="iid", hyper = hyperpc1))
      # }else{
      #     stop("hyper needs to be either pc or gamma.")
      # }
    }
    if(strata.time.effect){
        # In this case, age is age x strata already
        if(is.age.diff){
           formula <- update(formula, ~. -1 + age.intercept + age.diff) 
         }else{
            # In this case, age is age x strata already
            formula <- update(formula, ~. -1 + age.intercept)           
         } 
    }else{
        if(has.strata){
          # In this case, age is age x strata already
          if(is.age.diff){
             formula <- update(formula, ~. -1 + age.intercept + strata + age.diff) 
           }else{
              # In this case, age is age x strata already
              formula <- update(formula, ~. -1 + age.intercept + strata)
           } 
        }else{
          formula <- update(formula, ~. -1 + age.intercept) 
        }
    }
    if(!is.null(bias.adj)){
      if(is.null(bias.adj.by)) stop("bias.adj.by argument is not specified. Please specify which bias adjustment factors are specified by which columns. ")
      bias.adj <- bias.adj[, c(bias.adj.by, "ratio")]
      exdat$years <- as.character(exdat$years)
      exdat <- merge(exdat, bias.adj, all.x = TRUE, by = bias.adj.by)
      if("ratio" %in% colnames(exdat) == FALSE){
        stop("bias.adj argument is misspecified. It require the following column: ratio.")
      }
    }else{
      exdat$ratio <- 1
    }
    exdat$logoffset <- log(exdat$ratio)

    formula <- update(formula, ~. + offset(logoffset))
}




  total <- NA
  exdat <- subset(exdat, total != 0)
  
  # Create filler data frame for spatial-only models
  if(N == 0){
    tmp <- exdat[rep(1, dim(Amat)[1]), ]
    tmp$region.struct <- 1:dim(Amat)[1]
    tmp$region_number <- tmp$region.int <- tmp$region.unstruct <- tmp$region.struct
    tmp$region <- colnames(Amat)[tmp$region.struct]
    created = c("region.struct", "region_number", "region.unstruct", "region.int", "region", "years", "age", "age.intercept", "age.diff", "strata")
    tmp[, colnames(tmp) %in% created == FALSE] <- NA
    exdat <- rbind(exdat, tmp)
  }
  # Create filler data frame for all space-time pairs
  if(N >= 1){
    tmp <- exdat[rep(1, dim(Amat)[1]*N), ]
    tmp[, c("region.struct", "time.struct")] <- expand.grid(region.struct = 1:dim(Amat)[1], time.struct = 1:N)
    tmp$time.unstruct <- tmp$time.int <- tmp$time.struct
    created <- NULL
    if("time.slope" %in% colnames(exdat)){
      tmp$time.slope <- (tmp$time.unstruct  - center) / ((N - 1))
      created <- c(created, "time.slope")
    }
    if("st.slope" %in% colnames(exdat)){
      tmp$st.slope <- (tmp$time.unstruct  - center) / ((N - 1))   
      tmp$st.slope.id <- tmp$region.struct
      created <- c(created, "st.slope", "st.slope.id") 
    } 
    tmp$years <- years$year[match(tmp$time.struct, years$year_number)]
    tmp$region_number <- tmp$region.int <- tmp$region.unstruct <- tmp$region.struct
    tmp$region <- colnames(Amat)[tmp$region.struct]
    if(dim(Amat)[1] != 1 && is.temporal){
      tmp <- tmp[, colnames(tmp) != "time.area"]
      tmp <- merge(tmp, time.area, by = c("region_number", "time.unstruct"))
      tmp <- subset(tmp, tmp$time.area %in% exdat$time.area == FALSE)
    }
    # If need filler data
    if(dim(tmp)[1] > 0){
        # remove contents in other columns
        if(dim(Amat)[1] != 1){
          created <- c(created, "region.struct", "region_number", "region.unstruct", "region.int", "region", "time.struct", "time.unstruct", "time.int", "time.area", "years", "age", "age.intercept", "age.diff", "strata")
        }else{
            created <- c(created, "time.struct", "time.unstruct", "time.int",  "years", "age", "age.intercept", "age.diff", "strata")
        }
        tmp[, colnames(tmp) %in% created == FALSE] <- NA
        exdat <- rbind(exdat, tmp)
    }
  }

  if(!is.null(X)){
    exdat <- exdat[, colnames(exdat) %in% covariate.names == FALSE]
    if(!is.temporal){
          Xnew <- data.frame(region.struct = 1:S)
    }else{
      Xnew <- expand.grid(time.unstruct = 1:N, region.struct = 1:S)
      if(!is.spatial) Xnew$region.struct <- 0
    } 
    if("region" %in% by) X$region.struct <- match(X$region, colnames(Amat))
    if("years" %in% by) X$time.unstruct <- match(X$years, year_label)
    Xnew <- merge(Xnew, X, all.x = TRUE)
    if("region" %in% by && "years" %in% by){
      exdat <- merge(exdat, Xnew[, c(covariate.names, "time.unstruct", "region.struct")], by = c("time.unstruct", "region.struct"), all.x = TRUE)
    }else if("region" %in% by){
      exdat <- merge(exdat, Xnew[, c(covariate.names,"region.struct")], by = c("region.struct"), all.x = TRUE)
    }else if("years" %in% by){
      exdat <- merge(exdat, Xnew[, c(covariate.names,"region.struct")], by = c("time.unstruct"), all.x = TRUE)
    }
   }

  if(has.strata) exdat$strata <- factor(exdat$strata, levels = stratalevels)
  if(!is.null(age.groups)){
      exdat$age <- factor(exdat$age, levels = age.groups)
      if(is.age.diff){
        exdat$age.intercept <- factor(exdat$age.intercept, levels = unique(age.groups.vec))
        exdat$age.diff <- factor(exdat$age.diff, levels = age.diff.levels)
      }else{
          exdat$age.intercept <- factor(exdat$age.intercept, levels = unique(age.groups.new))
      }   
      
  }else{
     formula <- update(formula, ~. - age.intercept)  
     if(is.age.diff)  formula <- update(formula, ~. - age.diff)  
  }
  # if only one level, use the default intercept instead of age
  if(length(age.groups) == 1){
    formula <- update(formula, ~. - age.intercept + 1)   
    if(is.age.diff)  formula <- update(formula, ~. - age.diff)  
  }
  exdat$age.idx <- match(exdat$age, age.groups)
  exdat$age.rep.idx <- age.rw.group[exdat$age.idx]


if(!is.null(X)){
  formula <- as.formula(paste("Y~", as.character(formula)[3], "+", paste(covariate.names, collapse = " + ")))
}
if(family == "betabinomialna"){
    # ## Prepare compact version of the data frame
    # if(sum(!is.na(exdat$age.idx)) == 0) exdat$age.idx <- 1
    # if(sum(!is.na(exdat$age.rep.idx)) == 0) exdat$age.rep.idx <- 1
    # exdat$logoffset[is.na(exdat$logoffset)] <- 0

    # ## get list of covariates
    # if(survey.effect){
    #   rhs <- "years + region_number + time.unstruct + region + age + strata + region.int + region.unstruct + region.struct + time.int + time.struct + time.area + logoffset + age.idx + age.rep.idx + survey2" 
    # }else{
    #   rhs <- "years + region_number + time.unstruct + region + age + strata + region.int + region.unstruct + region.struct + time.int + time.struct + time.area + logoffset + age.idx + age.rep.idx"
    # }


    # # Total Y
    # tmp1 <- aggregate(as.formula(paste0("Y~", rhs)), data = exdat, FUN = function(x){sum(x, na.rm=TRUE)}, na.action = na.pass)
    # # Total N
    # tmp2 <- aggregate(as.formula(paste0("total~", rhs)), data = exdat, FUN = function(x){sum(x, na.rm=TRUE)}, na.action = na.pass)
    # # Scaling factor (\sum n(n-1))/N/N-1 
    # tmp3 <- aggregate(as.formula(paste0("total~", rhs)), data = exdat, FUN = function(x){sum(x*(x-1), na.rm=TRUE)/(sum(x, na.rm=TRUE)*(sum(x, na.rm=TRUE)-1)) }, na.action = na.pass)
    # colnames(tmp3)[which(colnames(tmp3) == "total")] <- "s"
    # exdat <- merge(merge(tmp1, tmp2), tmp3)
    # # Handle division by 0
    # if(sum(is.na(exdat$s)) > 0) exdat$s[is.na(exdat$s)] <- 1
    # exdat$Y[exdat$total == 0] <- NA
    # exdat$total[exdat$total == 0] <- NA

    # fit <- INLA::inla(formula, family = family, control.compute = options, data = exdat, control.predictor = list(compute = FALSE), Ntrials = exdat$total, scale = exdat$s, lincomb = NULL, control.inla =control.inla, verbose = verbose)

  }else{
    control.family <- NULL
    if(family == "betabinomial"){
      control.family <- list(hyper = list(rho = list(param = c(overdisp.mean, overdisp.prec))))
      fit <- INLA::inla(formula, family = family, control.compute = options, control.family = control.family, data = exdat, control.predictor = list(compute = FALSE, link=1), Ntrials = exdat$total, lincomb = NULL, control.inla = control.inla, verbose = verbose, ...)
    }else{
      fit <- INLA::inla(formula, family = family, control.compute = options, data = exdat, control.predictor = list(compute = FALSE, link=1), Ntrials = exdat$total, lincomb = NULL, control.inla = control.inla, verbose = verbose, ...)
    } 
  }  

 # find the name for baseline strata
 levels <- grep("strata", rownames(fit$summary.fixed))   
 levels <- gsub("strata", "", rownames(fit$summary.fixed)[levels])
 strata.all <- ""
 if(has.strata) strata.all <- as.character(unique(exdat$strata))
 strata.base <- strata.all[strata.all%in%levels == FALSE]

 priors <- list(pc.u = pc.u, pc.alpha = pc.alpha, pc.u.phi = pc.u.phi, pc.alpha.phi = pc.alpha.phi, pc.u.cor = pc.u.cor, pc.alpha.cor = pc.alpha.cor,  pc.st.u = pc.st.u, pc.st.alpha = pc.st.alpha, pc.st.slope.u = pc.st.slope.u, pc.st.slope.prec.u = pc.st.slope.prec.u, pc.st.slope.alpha = pc.st.slope.alpha, overdisp.mean = overdisp.mean, overdisp.prec = overdisp.prec)

 out <- list(model = formula, fit = fit, family= family, Amat = Amat, newdata = exdat, time = seq(0, N - 1), area = seq(0, region_count - 1), time.area = time.area, survey.table = survey.table, is.yearly = FALSE, type.st = type.st, year_label = year_label, age.groups = age.groups, age.groups.new = age.groups.new, age.n = age.n, age.rw.group = age.rw.group, age.strata.fixed.group = age.strata.fixed.group, strata.base = strata.base, rw = rw, ar = ar, strata.time.effect = strata.time.effect,  priors = priors, year_range = NA, Amat = Amat, has.Amat = TRUE, is.temporal = is.temporal, covariate.names = covariate.names, msg = msg)
 class(out) <- "SUMMERmodel"
 return(out)
    
  }
}

#' @export
#' @rdname smoothCluster
fitINLA2 <- smoothCluster