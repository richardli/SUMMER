#' Fit INLA models to direct estimators with a yearly model.
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
#' @param redo.prior Logical indicator of whether to re-estimate hyperparameters
#' @param priors priors from \code{\link{simhyper}}
#' @param useHyper option to manually set all hyperpriors
#' @param rw Take values 1 or 2, indicating the order of random walk.
#' @param is.yearly Logical indicator for fitting yearly or period model.
#' @param year_range Entire range of the years (inclusive) defined in year_names.
#' @param m Number of years in each period.
#' @param type.st type for space-time interaction
#' @param a.iid hyperparameter for i.i.d random effects, only need if \code{useHyper = TRUE}
#' @param b.iid hyperparameter for i.i.d random effects, only need if \code{useHyper = TRUE}
#' @param a.rw1 hyperparameter for RW1 random effects, only need if \code{useHyper = TRUE}
#' @param b.rw1 hyperparameter for RW1 random effects, only need if \code{useHyper = TRUE}
#' @param a.rw2 hyperparameter for RW2 random effects, only need if \code{useHyper = TRUE}
#' @param b.rw2 hyperparameter for RW2 random effects, only need if \code{useHyper = TRUE}
#' @param a.icar hyperparameter for ICAR random effects, only need if \code{useHyper = TRUE}
#' @param b.icar hyperparameter for ICAR random effects, only need if \code{useHyper = TRUE}
#' @seealso \code{\link{countrySummary}}
#' @importFrom stats dgamma
#' @importFrom Matrix Diagonal 
#' @return INLA model fit using the provided formula, country summary data, and geographic data
#' @examples
#' \dontrun{
#' 
#' data(DemoData)
#' data(DemoMap)
#' years <- levels(DemoData[[1]]$time)
#' 
#' # obtain direct estimates
#' data <- countrySummary_mult(births = DemoData, 
#' years = years, idVar = "id", 
#' regionVar = "region", timeVar = "time", 
#' clusterVar = "~clustid+id", 
#' ageVar = "age", weightsVar = "weights", 
#' geo.recode = NULL)
#' 
#' # obtain maps
#' geo <- DemoMap$geo
#' mat <- DemoMap$Amat
#' 
#' # Simulate hyperpriors
#' priors <- simhyper(R = 2, nsamp = 1e+05, nsamp.check = 5000, Amat = mat, only.iid = TRUE)
#' 
#' # combine data from multiple surveys
#' data <- aggregateSurvey(data)
#' 
#' # Model fitting with INLA
#' years.all <- c(years, "15-19")
#' fit <- fitINLA(data = data, geo = geo, Amat = mat, 
#' year_names = years.all, year_range = c(1985, 2019),
#'  priors = priors, rw = 2,
#'  is.yearly=TRUE, m = 5, type.st = 4)
#' # Projection
#' out <- projINLA(fit, Amat = mat, is.yearly = TRUE)
#' plot(out, is.yearly=TRUE, is.subnational=TRUE) + ggplot2::ggtitle("Subnational yearly model")
#' 
#' }
#' 
#' @export
fitINLA <- function(data, Amat, geo, formula = NULL, rw = 2, is.yearly = TRUE, year_names, year_range = c(1980, 2014), m = 5, na.rm = TRUE, redo.prior = FALSE, priors = NULL, type.st = 1, useHyper = FALSE, a.iid = NULL, b.iid = NULL, a.rw1 = NULL, b.rw1 = NULL, a.rw2 = NULL, b.rw2 = NULL, a.icar = NULL, b.icar = NULL){

  # check region names in Amat is consistent
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
  }

  # get around CRAN check of using un-exported INLA functions
  rate0 <- shape0 <- my.cache <- inla.as.sparse <- type <- NULL

  if (!isTRUE(requireNamespace("INLA", quietly = TRUE))) {
    stop("You need to install the packages 'INLA'. Please run in your R terminal:\n install.packages('INLA', repos='https://www.math.ntnu.no/inla/R/stable')")
  }
  # If INLA is installed, then attach the Namespace (so that all the relevant functions are available)
  if (isTRUE(requireNamespace("INLA", quietly = TRUE))) {
    if (!is.element("INLA", (.packages()))) {
      attachNamespace("INLA")
    }


    tau = exp(10)

    ## ---------------------------------------------------------
    ## New definition of the yearly + multi-year Q structure
    ## ---------------------------------------------------------
     rw.new = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
      ## assume 'tau', 'order', 'n' and 'm' 'n' is the dim of RW and 'm' is the aggregated length,
      ## averaging over n/m variables, non-overlapping
      
      ## the environment of this function which holds the variables and we can store 'my.cache'
      ## there.
      envir = environment(sys.call()[[1]]) 
      if(is.null(envir)) envir <- environment()
      inla.rw = utils::getFromNamespace("inla.rw", "INLA")
      inla.ginv = utils::getFromNamespace("inla.ginv", "INLA")
      
      if (!exists("my.cache", envir = envir, mode = "list")) {
        nn = n %/% m
        stopifnot (nn == as.integer(n/m))
        R = inla.rw(n, order = order,  scale.model=TRUE, sparse=TRUE)
        A = matrix(0, nn, n)
        j = 1
        for(i in 1:nn) {
          A[i, j:(j+m-1)] = 1/m
          j = j + m
        }
        A = inla.as.sparse(A)
        D = Matrix::Diagonal(nn, x=1)
        assign("my.cache", list(R=R, A=A, D=D, nn=nn), envir = envir)
      } 
      
      interpret.theta = function() {
        return(list(kappa = exp(theta[1L])))
      }
      
      graph = function() {
        return (Q())
      }
      
      Q = function() {
        QQ = rbind(cbind(p$kappa * my.cache$R + tau * t(my.cache$A) %*% my.cache$A,
                         -tau * t(my.cache$A)),
                   cbind(-tau * my.cache$A, tau * my.cache$D))
        return(QQ)
      }
      
      mu = function() {
        return(numeric(0))
      }
      
      log.norm.const = function() {
        val = (n-order) * (-0.5 * log(2 * pi) + 0.5 * log(p$kappa)) +
          (my.cache$nn * (-0.5 * log(2 * pi) + 0.5 * log(tau)))
        return(val)
      }
      
      log.prior = function() {
        val = dgamma(p$kappa, shape = shape0, rate = rate0, log = TRUE) + theta[1]
        return(val)
      }
      
      initial = function() {
        return(4)
      }
      
      quit = function() {
        return(invisible())
      }
      
      ## as some calls to this function does not define 'theta',  its convenient to have to
      ## defined still (like in the graph-function)
      if (is.null(theta))
        theta = initial()
      
      p = interpret.theta()
      val = do.call(match.arg(cmd), args = list())
      return(val)
     }   
    
    ## ---------------------------------------------------------
    ## New definition of the yearly + multi-year Q structure
    ## ---------------------------------------------------------
      iid.new = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
      
      envir = environment(sys.call()[[1]]) 
            if(is.null(envir)) envir <- environment()
      if (!exists("my.cache", envir = envir, mode = "list")) {
        nn = n %/% m
        stopifnot (nn == as.integer(n/m))
        R = Diagonal(n, x = rep(1, n))
        A = matrix(0, nn, n)
        j = 1
        for(i in 1:nn) {
          A[i, j:(j+m-1)] = 1/m
          j = j + m
        }
        A = inla.as.sparse(A)
        D = Matrix::Diagonal(nn, x=1)
        assign("my.cache", list(R=R, A=A, D=D, nn=nn), envir = envir)
      } 
      
      interpret.theta = function() {
        return(list(kappa = exp(theta[1L])))
      }
      
      graph = function() {
        return (Q())
      }
      
      Q = function() {
        QQ = rbind(cbind(p$kappa * my.cache$R + tau * t(my.cache$A) %*% my.cache$A,
                         -tau * t(my.cache$A)),
                   cbind(-tau * my.cache$A, tau * my.cache$D))
        return(QQ)
      }
      
      mu = function() {
        return(numeric(0))
      }
      
      log.norm.const = function() {
        val = (n * (-0.5 * log(2 * pi) + 0.5 * log(p$kappa)) +
          (my.cache$nn * (-0.5 * log(2 * pi) + 0.5 * log(tau))))
        return(val)
      }
      
      log.prior = function() {
        val = dgamma(p$kappa, shape = shape0, rate = rate0, log = TRUE) + theta[1]
        return(val)
      }
      
      initial = function() {
        return(4)
      }
      
      quit = function() {
        return(invisible())
      }
      
      ## as some calls to this function does not define 'theta',  its convenient to have to
      ## defined still (like in the graph-function)
      if (is.null(theta))
        theta = initial()
      
      p = interpret.theta()
      val = do.call(match.arg(cmd), args = list())
      return(val)
    }  

    
    ## ---------------------------------------------------------
    ## New definition of the yearly + multi-year structured Q
    ## ---------------------------------------------------------
    st.new = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
    
    envir = environment(sys.call()[[1]]) 
          if(is.null(envir)) envir <- environment()
    inla.rw = utils::getFromNamespace("inla.rw", "INLA")
    inla.ginv = utils::getFromNamespace("inla.ginv", "INLA")
    # The new structure takes the following order
    # (x_11, ..., x_1T, ..., x_S1, ..., x_ST, xx_11, ..., xx_1t, ..., xx_S1, ..., xx_St)
    #  x_ij : random effect of region i, year j 
    # xx_ik : random effect of region i, period k

    if (!exists("my.cache", envir = envir, mode = "list")) {
      nn = n %/% m
      stopifnot (nn == as.integer(n/m))
      R1 = Diagonal(n, x = rep(1, n))
      R2 = inla.rw(n, order = order, scale.model=TRUE, sparse=TRUE)
      R3 = Diagonal(S, x = rep(1, S))
      R4 = Amat
      diag(R4) <- 0
      diag <- apply(R4, 1, sum)
      R4[R4 != 0] <- -1
      diag(R4) <- diag

      R4 <- INLA::inla.scale.model(R4, constr = list(A=matrix(1,1,dim(R4)[1]), e=0))
      # both independent
      if(type == 1){
          R <- R3 %x% R1
      # AR * independent    
      }else if(type == 2){
          R <- R3 %x% R2
      # independent * besag    
      }else if(type == 3){
          R <- R4 %x% R1
      # AR * besag
      }else if(type == 4){
          R <- R4 %x% R2
      }

      A = matrix(0, nn*S, n*S)
      j = 1
      for(i in 1:(nn*S)) {
        A[i, j:(j+m-1)] = 1/m
        j = j + m
      }
      A = inla.as.sparse(A)
      D = Diagonal(nn*S, x=1)
      assign("my.cache", list(R=INLA::inla.as.sparse(R), A=A, D=D, nn=nn), envir = envir)
    } 
    
    interpret.theta = function() {
      return(list(kappa = exp(theta[1L])))
    }
    
    graph = function() {
      return (Q())
    }
    
    Q = function() {
      QQ = rbind(cbind(p$kappa * my.cache$R + tau * t(my.cache$A) %*% my.cache$A,
                         -tau * t(my.cache$A)),
                   cbind(-tau * my.cache$A, tau * my.cache$D))
      return(QQ)
    }
    
    mu = function() {
      return(numeric(0))
    }
    ## Type I   : S * n
    ## Type II  : S * (n - order)
    ## Type III : (S-1) * n 
    ## Type IV  : (S-1) * (n - order)
    log.norm.const = function() {
      df <- S * n
      if(type == 2){
        df <- S * (n - order)
      }else if(type == 3){
        df <- (S-1) * n
      }else if(type == 4){
        df <- (S-1) * (n - order)
      }
      val = (df * (-0.5 * log(2 * pi) + 0.5 * log(p$kappa)) +
        (S * my.cache$nn * (-0.5 * log(2 * pi) + 0.5 * log(tau))))
      return(val)
    }
    
    log.prior = function() {
      val = dgamma(p$kappa, shape = shape0, rate = rate0, log = TRUE) + theta[1]
      return(val)
    }
    
    initial = function() {
      return(4)
    }
    
    quit = function() {
      return(invisible())
    }
    
    ## as some calls to this function does not define 'theta',  its convenient to have to
    ## defined still (like in the graph-function)
    if (is.null(theta))
      theta = initial()
    
    p = interpret.theta()
    val = do.call(match.arg(cmd), args = list())
    return(val)
  }  

    
    ## ---------------------------------------------------------
    ## Common Setup
    ## --------------------------------------------------------- 
    if(is.null(geo)){
      data <- data[which(data$region == "All"), ]
      if(length(data) == 0){
        stop("No geographics specified and no observation labeled 'All' either.")
      }
    } else{
      data <- data[which(data$region != "All"), ]
    }  
    #################################################################### Re-calculate hyper-priors
    # Todo: make it work with the new Q matrix!!
    
    if (redo.prior) {
      priors <- simhyper(R = 2, nsamp = 1e+05, nsamp.check = 5000, Amat = Amat, nperiod = length(year_names))
    }
    
    a.iid <- priors$a.iid
    b.iid <- priors$b.iid
    a.rw1 <- priors$a.iid
    b.rw1 <- priors$a.iid
    a.rw2 <- priors$a.iid
    b.rw2 <- priors$a.iid
    a.icar <- priors$a.iid
    b.icar <- priors$a.iid
    
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
    if(is.null(geo)){
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
    dat$region.struct <- dat$region.unstruct <- dat$region_number
    
    ################################################################### get the lsit of region and numeric index in one data frame
    if(is.yearly){
      n <- year_range[2] - year_range[1] + 1
      nn <- n %/% m
      N <- n + nn
      rw.model <- INLA::inla.rgeneric.define(model = rw.new,
                                       n = n, 
                                       m = m,
                                       order = rw,
                                       tau = exp(10),
                                       shape0 = a.rw2,
                                       rate0 = b.rw2) 
      iid.model.time <- INLA::inla.rgeneric.define(model = iid.new,
                                             n = n, 
                                             m = m,
                                             tau = exp(10),
                                             shape0 = a.iid,
                                             rate0 = b.iid)
      if(!is.null(geo)){
         st.model <- INLA::inla.rgeneric.define(model = st.new,
                                       n = n, 
                                       m = m,
                                       order = rw,
                                       S = region_count,
                                       Amat = Amat,
                                       type = type.st,
                                       tau = exp(10),
                                       shape0 = a.iid,
                                       rate0 = b.iid)
       }
      
      year_names_new <- c(as.character(c(year_range[1]:year_range[2])), year_names)
      time.index <- cbind.data.frame(idx = 1:N, Year = year_names_new)
      constr = list(A = matrix(c(rep(1, n), rep(0, nn)), 1, N), e = 0)
      
      if(type.st %in% c(2, 4)){
        tmp <- matrix(0, S, N * S)
        for(i in 1:S){
          tmp[i, ((i-1)*n + 1) : (i*n)] <- 1
        }
      }else{
        tmp <- NULL
      }
      
      # ICAR constraints
      if(type.st %in% c(3, 4)){
        tmp2 <- matrix(0, n, N*S)
        for(i in 1:n){
          tmp2[i , which((1:(n*S)) %% n == i-1)] <- 1
        }
      }else{
        tmp2 <- NULL
      }
      tmp <- rbind(tmp, tmp2)
      if(is.null(tmp)){
        constr.st <- NULL
      }else{
        constr.st <- list(A = tmp, e = rep(0, dim(tmp)[1]))
      }
      years <- data.frame(year = year_names_new[1:N], year_number = seq(1, N))
    }else{
      n <- 0
      N <- nn <- length(year_names)
      years <- data.frame(year = year_names, year_number = seq(1, N))      
    }
    
    # -- creating IDs for the temporal REs -- #
    if(is.yearly){
      dat$time.unstruct <- dat$time.struct <- years[match(dat$years, years[, 1]), 2]
    }else{
      dat$time.unstruct <- dat$time.struct <- years[match(dat$years, years[, 1]), 2]
    }
    
    ################################################################## get the number of surveys
    if(sum(!is.na(data$survey)) == 0){
      data$survey <- 1
      nosurvey <- TRUE
    }else{
      nosurvey <- FALSE
    }
    survey_count <- length(table(data$survey))
    ################################################################## -- these are the time X survey options -- #
    x <- expand.grid(1:nn, 1:survey_count)
    survey.time <- data.frame(time.unstruct = x[, 1], survey = x[, 2], survey.time = c(1:nrow(x)))
    
    # -- these are the area X survey options -- #
    x <- expand.grid(1:region_count, 1:survey_count)
    survey.area <- data.frame(region_number = x[, 1], survey = x[, 2], survey.area = c(1:nrow(x)))
    
    # -- these are the area X time options -- #
    # The new structure takes the following order
    # (x_11, ..., x_1T, ..., x_S1, ..., x_ST, xx_11, ..., xx_1t, ..., xx_S1, ..., xx_St)
    #  x_ij : random effect of region i, year j 
    # xx_ik : random effect of region i, period k
    if(is.yearly){
      x <- rbind(expand.grid(1:n, 1:region_count), 
                 expand.grid((n+1):N, 1:region_count))
    }else{
      x <- expand.grid(1:N, 1:region_count)
    }
    time.area <- data.frame(region_number = x[, 2], time.unstruct = x[, 1], time.area = c(1:nrow(x)))
    # fix for 0 instead of 1 when no geo file provided
    if(is.null(geo)){
      time.area$region_number <- 0
    }
    # -- these are the area X time X survey options -- #
    x <- expand.grid(1:region_count, 1:N, 1:survey_count)
    survey.time.area <- data.frame(region_number = x[, 1], time.unstruct = x[, 2], survey = x[, 3], survey.time.area = c(1:nrow(x)))
    
    # -- merge these all into the data sets -- #
    newdata <- dat
    if (!nosurvey) {
      newdata <- merge(newdata, survey.time, by = c("time.unstruct", "survey"))
      newdata <- merge(newdata, survey.area, by = c("region_number", "survey"))
      newdata <- merge(newdata, survey.time.area, by = c("region_number", "time.unstruct", "survey"))
    }
    if(!is.null(geo)){
      newdata <- merge(newdata, time.area, by = c("region_number", "time.unstruct"))
    }else{
      newdata$time.area <- NA
    }
    
    
    ########################## Model Selection ######
    
    # -- subset of not missing and not direct estimate of 0 -- #
    exdat <- newdata
    exdat <- exdat[!is.na(exdat$logit.est) && exdat$logit.est > (-20), ]
    
    
    ## ---------------------------------------------------------
    ## Setup yearly model
    ## ---------------------------------------------------------
    if(is.yearly && (!is.null(geo))){   
      if (is.null(formula)) {
        if(rw == 1){
          formula <- logit.est ~ f(time.struct, model = rw.model, diagonal = 1e-6, extraconstr = constr, values = 1:N) + f(region.unstruct,model="iid",param=c(a.iid,b.iid)) + f(region.struct, graph=Amat,model="besag",param=c(a.icar,b.icar), scale.model = TRUE) + f(time.unstruct,model=iid.model.time) + f(time.area,model=st.model, diagonal = 1e-6, extraconstr = constr.st, values = 1:(N*S))
        }else if(rw == 2 && type.st %in% c(2, 4)){
          formula <- logit.est ~ f(time.struct, model = rw.model, diagonal = 1e-6, extraconstr = constr, values = 1:N) + f(region.unstruct,model="iid",param=c(a.iid,b.iid)) + f(region.struct, graph=Amat,model="besag",param=c(a.icar,b.icar), scale.model = TRUE) + f(time.unstruct,model=iid.model.time) + f(time.area,model=st.model, diagonal = 1e-6, extraconstr = constr.st, values = 1:(N*S))
        }else if(rw == 2){
          formula <- logit.est ~ f(time.struct, model = rw.model, diagonal = 1e-6, extraconstr = constr, values = 1:N) + f(region.unstruct,model="iid",param=c(a.iid,b.iid)) + f(region.struct, graph=Amat,model="besag",param=c(a.icar,b.icar), scale.model = TRUE) + f(time.unstruct,model=iid.model.time) + f(time.area,model=st.model, diagonal = 1e-6, extraconstr = constr.st, values = 1:(N*S))
          
        }else{
          stop("Random walk order should be 1 or 2.")
        }
      }
      
      
      ## ---------------------------------------------------------
      ## Setup non-yearly model
      ## ---------------------------------------------------------
    }else if((!is.yearly) && (!is.null(geo))){
      if (is.null(formula)) {
        if(rw == 1){
          formula <- logit.est ~ f(region.unstruct,model="iid",param=c(a.iid,b.iid)) + f(region.struct, graph=Amat,model="besag",param=c(a.icar,b.icar), scale.model = TRUE) + f(time.struct,model="rw1",param=c(a.rw1,b.rw1))  + f(time.unstruct,model="iid",param=c(a.iid,b.iid)) + f(time.area,model="iid", param=c(a.iid,b.iid))
        }else if(rw == 2){
          formula <- logit.est ~ f(region.unstruct,model="iid",param=c(a.iid,b.iid)) + f(region.struct, graph=Amat,model="besag",param=c(a.icar,b.icar), scale.model = TRUE) + f(time.struct,model="rw2",param=c(a.rw2,b.rw2))  +
            f(time.unstruct,model="iid",param=c(a.iid,b.iid)) + f(time.area,model="iid", param=c(a.iid,b.iid))
        }else{
          stop("Random walk order should be 1 or 2.")
        }
      }
      ## ---------------------------------------------------------
      ## Setup yearly national model
      ## --------------------------------------------------------- 
    }else if(is.yearly && is.null(geo)){   
      if (is.null(formula)) {
        if(rw == 1){
          formula <- logit.est ~ f(time.struct, model = rw.model, diagonal = 1e-6, extraconstr = constr, values = 1:N) + f(time.unstruct,model=iid.model.time) 
        }else if(rw == 2){
          formula <- logit.est ~ f(time.struct, model = rw.model, diagonal = 1e-6, extraconstr = constr, values = 1:N) + f(time.unstruct,model=iid.model.time)
        }else{
          stop("Random walk order should be 1 or 2.")
        }
      }
      
      ## ---------------------------------------------------------
      ## Setup non-yearly national model
      ## ---------------------------------------------------------
    }else if((!is.yearly) && (is.null(geo))){
      if (is.null(formula)) {
        if(rw == 1){
          formula <- logit.est ~ f(time.struct,model="rw1",param=c(a.rw1,b.rw1))  + f(time.unstruct,model="iid",param=c(a.iid,b.iid), scale.model = TRUE) 
        }else if(rw == 2){
          formula <- logit.est ~ f(time.struct,model="rw2",param=c(a.rw2,b.rw2))  +
            f(time.unstruct,model="iid",param=c(a.iid,b.iid)) 
        }else{
          stop("Random walk order should be 1 or 2.")
        }
      }
    }  
    mod <- formula
    
    
    
    ## ---------------------------------------------------------
    ## Subnational lincomb for projection
    ## ---------------------------------------------------------
    if(!is.null(geo)){
      lincombs.info <- data.frame(Index = 1:(region_count*N), District = NA, Year = NA)
      index <- 0
      for(j in 1:region_count){
        for(i in 1:N){
          index <- index + 1    
          time <- rep(NA, N)
          # time.old <- rep(NA, m)
          area <- rep(NA, region_count)
          spacetime <- rep(NA, N*region_count) 
          
          space.time.id <- unique(time.area$time.area[time.area$time.unstruct == i & time.area$region_number == j])
          spacetime[space.time.id] <- 1
          time[i] <- 1
          area[j] <- 1
          time.unstruct <- time
          
          object.name <- paste("lc", index, sep = "")
          
          lincombs.info[index, c("District", "Year")] <- c(j,i)
          if(rw == 1){
            assign(object.name, INLA::inla.make.lincomb("(Intercept)" = 1,
                                                  time.area = spacetime,
                                                  time.struct= time ,
                                                  time.unstruct= time,
                                                  region.struct = area,
                                                  region.unstruct = area))          
          }else if(is.yearly && type.st %in% c(2, 4) && rw == 2){
            # the name of the third argument is changed later
            assign(object.name, INLA::inla.make.lincomb("(Intercept)" = 1,
                                                  time.area = spacetime,
                                                  time.struct= time ,
                                                  time.unstruct= time,
                                                  region.struct = area,
                                                  region.unstruct = area))
          }else{
            assign(object.name, INLA::inla.make.lincomb("(Intercept)" = 1,
                                                  time.area = spacetime,
                                                  time.struct= time ,
                                                  time.unstruct= time,
                                                  region.struct = area,
                                                  region.unstruct = area))
          }
          
          if(index == 1){
            lincombs.yearly <- get(object.name)
            names(lincombs.yearly)[index] <- object.name
          }else{
            tmp <- get(object.name)
            lincombs.yearly <- c(lincombs.yearly, tmp)
            names(lincombs.yearly)[index] <- object.name
          }
        }
      }
      
      ##------------------------------------------------------------##
      ## National model lincomb for projection
      ##------------------------------------------------------------##
    }else{
      lincombs.info <- data.frame(Index = 1:N, District = NA, Year = NA)
      index <- 0
      for(i in 1:N){
        index <- index + 1    
        time <- rep(NA, N)
        time[i] <- 1
        time.unstruct <- time
        
        object.name <- paste("lc", index, sep = "")
        
        lincombs.info[index, c("District", "Year")] <- c(0,i)
        if(rw == 1){
          assign(object.name, INLA::inla.make.lincomb("(Intercept)" = 1,
                                                time.struct= time ,
                                                time.unstruct= time))          
        }else{
          assign(object.name, INLA::inla.make.lincomb("(Intercept)" = 1,
                                                time.struct= time ,
                                                time.unstruct= time))
        }
        
        if(index == 1){
          lincombs.yearly <- get(object.name)
          names(lincombs.yearly)[index] <- object.name
        }else{
          lincombs.yearly <- c(lincombs.yearly, get(object.name))
          names(lincombs.yearly)[index] <- object.name
        }
      }
    }
    
    # if(is.yearly){
    # rbind yearly data with NA for the lincombs
    for(i in 1:N){
      tmp<-exdat[match(unique(exdat$region), exdat$region), ]
      tmp$time.unstruct<-tmp$time.struct<- i
      tmp$logit.est<-tmp$logit.prec<-tmp$survey<-NA
      tmp <- tmp[, colnames(tmp) != "time.area"]
      tmp <- merge(tmp, time.area, by = c("region_number", "time.unstruct"))
      tmp$years<-years[i, 1]
      tmp$u5m <- tmp$lower <- tmp$upper <- tmp$var.est <- NA
      if("u5m.nohiv" %in% colnames(data)){
        tmp$u5m.nohiv <- tmp$lower.nohiv <- tmp$upper.nohiv <- tmp$var.est.nohiv<- tmp$logit.prec.nohiv<- tmp$logit.est.nohiv <- NA          
      }
      exdat<-rbind(exdat,tmp)   
    }
    # }
    
    
    
    
    # -- fitting the model in INLA -- #

      inla11 <- INLA::inla(mod, family = "gaussian", control.compute = list(dic = T, mlik = T, cpo = T), data = exdat, control.predictor = list(compute = TRUE), control.family = list(hyper= list(prec = list(initial= log(1), fixed= TRUE ))), scale = exdat$logit.prec, 
                           lincomb = lincombs.yearly)
    
    return(list(model = mod, fit = inla11, Amat = Amat, newdata = exdat, time = seq(0, N - 1), area = seq(0, region_count - 
                                                                                                            1), survey.time = survey.time, survey.area = survey.area, time.area = time.area, survey.time.area = survey.time.area, 
                a.iid = a.iid, b.iid = b.iid, a.rw1 = a.rw1, b.rw1 = b.rw1, a.rw2 = a.rw2, b.rw2 = b.rw2, a.icar = a.icar, b.icar = b.icar, lincombs.info = lincombs.info))
    
  }
  
  
}