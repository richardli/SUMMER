#' Fit space-time smoothing models to mortality rates 
#' 
#' 
#' 
#'
#' @param data Combined dataset
#' @param geo Spatial polygon file, legacy parameter from previous versions of the package.
#' @param Amat Adjacency matrix for the regions
#' @param X Covariate matrix. It must contain either a column with name "region", or a column with name "years", or both. The covariates must not have missing values for all regions (if varying in space) and all time periods (if varying in time). The rest of the columns are treated as covariates in the mean model.
#' @param formula INLA formula. See vignette for example of using customized formula.
#' @param year_label string vector of year names
#' @param na.rm Logical indicator of whether to remove rows with NA values in the data. Default set to TRUE.
#' @param priors priors from \code{\link{simhyper}}
#' @param rw Take values 0, 1 or 2, indicating the order of random walk. If rw = 0, the autoregressive process is used instead of the random walk in the main trend. See the description of the argument ar for details.
#' @param ar Order of the autoregressive component. If ar is specified to be positive integer, the random walk components will be replaced by AR(p) terms in the interaction part. The main temporal trend remains to be random walk of order rw unless rw = 0.
#' @param is.yearly Logical indicator for fitting yearly or period model.
#' @param year_range Entire range of the years (inclusive) defined in year_label.
#' @param m Number of years in each period.
#' @param type.st type for space-time interaction
#' @param survey.effect logical indicator whether to include a survey iid random effect. If this is set to TRUE, there needs to be a column named 'survey' in the input data frame. In prediction, this random effect term will be set to 0. 
#' @param hyper which hyperpriors to use. Default to be using the PC prior ("pc"). 
#' @param pc.u hyperparameter U for the PC prior on precisions.
#' @param pc.alpha hyperparameter alpha for the PC prior on precisions.
#' @param pc.u.phi hyperparameter U for the PC prior on the mixture probability phi in BYM2 model.
#' @param pc.alpha.phi hyperparameter alpha for the PC prior on the mixture probability phi in BYM2 model.
#' @param pc.u.cor hyperparameter U for the PC prior on the autocorrelation parameter in the AR prior, i.e. Prob(cor > pc.u.cor) = pc.alpha.cor.
#' @param pc.alpha.cor hyperparameter alpha for the PC prior on the autocorrelation parameter in the AR prior.
#' @param pc.st.u hyperparameter U for the PC prior on precisions for the interaction term.
#' @param pc.st.alpha hyperparameter alpha for the PC prior on precisions for the interaction term.
#' @param a.iid hyperparameter for i.i.d random effects.
#' @param b.iid hyperparameter for i.i.d random effects.
#' @param a.rw hyperparameter for RW 1 or 2 random effects.
#' @param b.rw hyperparameter for RW 1 or 2random effects.
#' @param a.icar hyperparameter for ICAR random effects.
#' @param b.icar hyperparameter for ICAR random effects.
#' @param control.inla list of options to be passed to control.inla() in the inla() function. Default to the "adaptive" integration strategy.
#' @param options list of options to be passed to control.compute() in the inla() function.
#' @param verbose logical indicator to print out detailed inla() intermediate steps.
#' @seealso \code{\link{getDirect}}
#' @import Matrix
#' @importFrom stats dgamma
#' @importFrom Matrix Diagonal 
#' @return INLA model fit using the provided formula, country summary data, and geographic data
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
#' fit1 <- fitINLA(data = data, Amat = NULL, 
#'   year_label = years.all, year_range = c(1985, 2019), 
#'   rw = 2, is.yearly=FALSE, m = 5)
#' out1 <- getSmoothed(fit1)
#' plot(out1, is.subnational=FALSE)
#' 
#' #  subnational model
#' fit2 <- fitINLA(data = data, Amat = mat, 
#'   year_label = years.all, year_range = c(1985, 2019), 
#'   rw = 2, is.yearly=TRUE, m = 5, type.st = 4)
#' out2 <- getSmoothed(fit2, Amat = mat)
#' plot(out2, is.yearly=TRUE, is.subnational=TRUE)
#' 
#' 
#' }
#' @export
fitINLA <- function(data, Amat, X = NULL, formula = NULL, rw = 2, ar = 0, is.yearly = TRUE, year_label, year_range = c(1980, 2014), m = 5, na.rm = TRUE, priors = NULL, type.st = 1, survey.effect = FALSE, hyper = c("pc", "gamma")[1], pc.u = 1, pc.alpha = 0.01, pc.u.phi = 0.5, pc.alpha.phi = 2/3, pc.u.cor = 0.7, pc.alpha.cor = 0.9, a.iid = NULL, b.iid = NULL, a.rw = NULL, b.rw = NULL, a.icar = NULL, b.icar = NULL, options = list(dic = T, mlik = T, cpo = T, openmp.strategy = 'default'), control.inla = list(strategy = "adaptive", int.strategy = "auto"), verbose = FALSE, pc.st.u = NA, pc.st.alpha = NA, geo = NULL){

  if(!is.null(geo)){
    message("Argument geo is no longer necessary in the fitINLA2 function. Only Amat is needed.")
  }

  is.ar <- ar > 0 
  is.main.ar <- rw == 0
  if(is.ar) message("ar > 0: using the AR(p) process for the space-time interaction component.")
  if(rw %in% c(0, 1, 2) == FALSE) stop("Random walk only support rw = 1 or 2.")
  if(rw == 0) message("rw = 0: using the AR(p) process for the main temporal trend component.")

  if(m == 1){
    if(is.yearly) warning("Switched to period model because m = 1.")
    is.yearly = FALSE
  }
  if(is.ar && hyper=="gamma"){
    stop("AR1 model only implemented with PC priors for now.")
  }


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
  if (!is.element("Matrix", (.packages()))) {
    attachNamespace("Matrix")
  }
  # If INLA is installed, then attach the Namespace (so that all the relevant functions are available)
  if (isTRUE(requireNamespace("INLA", quietly = TRUE))) {
    if (!is.element("INLA", (.packages()))) {
      attachNamespace("INLA")
    }


    tau = exp(10)

   
    
    ## ---------------------------------------------------------
    ## Common Setup
    ## --------------------------------------------------------- 
    if(is.null(Amat)){
      data <- data[which(data$region == "All"), ]
      if(length(data) == 0){
        stop("Spatial adjacency matrix is not provided and no observation labeled 'All' either.")
      }
    } else{
      data <- data[which(data$region != "All"), ]
    }  


    #################################################################### Re-calculate hyper-priors
    
    if (is.null(priors)) {
      priors <- simhyper(R = 2, nsamp = 1e+05, nsamp.check = 5000, Amat = Amat, nperiod = length(year_label), only.iid = TRUE)
    }
    
    if(is.null(a.iid)) a.iid <- priors$a.iid
    if(is.null(b.iid)) b.iid <- priors$b.iid
    if(is.null(a.rw)) a.rw <- priors$a.iid
    if(is.null(b.rw)) b.rw <- priors$b.iid
    if(is.null(a.icar)) a.icar <- priors$a.iid
    if(is.null(b.icar)) b.icar <- priors$b.iid
    
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
    if(is.null(Amat)){
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
    
    ################################################################### get the lsit of region and numeric index in one data frame
    if(is.yearly){
      n <- year_range[2] - year_range[1] + 1
      if(n %% m  != 0) stop("The year_range specification is not a multiple of m. Please check that  year_range contains the first and last  year of the periods used. For example, if the last period is 2015-2019, the second element in year_range should be 2019.")
      nn <- n %/% m
      N <- n + nn
      rw.model <- INLA::inla.rgeneric.define(model = rw.new,
                                       n = n, 
                                       m = m,
                                       order = rw,
                                       tau = exp(10),
                                       shape0 = a.rw,
                                       rate0 = b.rw) 
      iid.model <- INLA::inla.rgeneric.define(model = iid.new,
                                             n = n, 
                                             m = m,
                                             tau = exp(10),
                                             shape0 = a.iid,
                                             rate0 = b.iid)

      rw.model.pc <- INLA::inla.rgeneric.define(model = rw.new.pc,
                                       n = n, 
                                       m = m,
                                       order = rw,
                                       tau = exp(10),
                                       u0 = pc.u,
                                       alpha0 = pc.alpha) 
      iid.model.pc <- INLA::inla.rgeneric.define(model = iid.new.pc,
                                             n = n, 
                                             m = m,
                                             tau = exp(10),
                                             u0 = pc.u,
                                             alpha0 = pc.alpha) 
      if(!is.null(Amat)){
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
         st.model.pc <- INLA::inla.rgeneric.define(model = st.new.pc,
                                       n = n, 
                                       m = m,
                                       order = rw,
                                       S = region_count,
                                       Amat = Amat,
                                       type = type.st,
                                       tau = exp(10),
                                       u0 = pc.u,
                                       alpha0 = pc.alpha) 
       }
      
      year_label_new <- c(as.character(c(year_range[1]:year_range[2])), year_label)
      time.index <- cbind.data.frame(idx = 1:N, Year = year_label_new)
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
      years <- data.frame(year = year_label_new[1:N], year_number = seq(1, N))
    }else{
      n <- 0
      N <- nn <- length(year_label)
      years <- data.frame(year = year_label, year_number = seq(1, N))      
    }
    
    # -- creating IDs for the temporal REs -- #
    if(is.yearly){
      dat$time.unstruct <- dat$time.struct <- dat$time.int <- years[match(dat$years, years[, 1]), 2]
    }else{
      dat$time.unstruct <- dat$time.struct <- dat$time.int <- years[match(dat$years, years[, 1]), 2]
    }
    

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
    if(is.null(Amat)){
      time.area$region_number <- 0
    }
    ################################################################## get the number of surveys
    if(survey.effect){
          if("survey" %in% colnames(data) == FALSE) stop("survey.effect is set to TRUE, but no survey column in the input data")
          survey_count <- length(table(data$survey))
          if(survey_count <= 1){
            warning("survey.effect is set to TRUE, but only one survey index is provided, the survey effect has been removed.")
            survey.effect <- FALSE
          }
    }
      
    # -- merge these all into the data sets -- #
    newdata <- dat
    if (survey.effect) {
      newdata$survey.id <- match(newdata$survey, unique(newdata$survey))
      survey.table <- data.frame(survey = unique(newdata$survey), survey.id = 1:survey_count)
    }else{
      newdata$survey.id <- NA
      survey.table <- NULL
    }
    if(!is.null(Amat)){
      newdata <- merge(newdata, time.area, by = c("region_number", "time.unstruct"))
    }else{
      newdata$time.area <- NA
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
    }


    
    ########################## Model Selection ######
    
    ## -- subset of not missing and not direct estimate of 0 -- #
    exdat <- newdata
    # exdat <- exdat[!is.na(exdat$logit.est) && exdat$logit.est > (-20), ]
    for(i in 1:N){
      tmp<-exdat[match(unique(exdat$region), exdat$region), ]
      tmp$time.unstruct<-tmp$time.struct<- tmp$time.int <- i
      tmp$logit.est<-tmp$logit.prec<-tmp$survey<-tmp$survey.id <- NA
      tmp <- tmp[, colnames(tmp) != "time.area"]
      tmp <- merge(tmp, time.area, by = c("region_number", "time.unstruct"))
      tmp$years<-years[i, 1]
      tmp$mean <- tmp$lower <- tmp$upper <- tmp$var.est <- NA
      if("mean.nohiv" %in% colnames(data)){
        tmp$mean.nohiv <- tmp$lower.nohiv <- tmp$upper.nohiv <- tmp$var.est.nohiv<- tmp$logit.prec.nohiv<- tmp$logit.est.nohiv <- NA          
      }
      exdat<-rbind(exdat,tmp)   
    }
   if(!is.null(X)){
    exdat <- exdat[, colnames(exdat) %in% covariate.names == FALSE]
    Xnew <- expand.grid(time.struct = 1:N, region.struct = 1:S)
    if(is.null(Amat)) Xnew$region.struct <- 0
    Xnew[, covariate.names] <- NA
    for(i in 1:dim(Xnew)[1]){
      tt <- ifelse(Xnew$time.struct[i] > n, Xnew$time.struct[i] - n, (Xnew$time.struct[i]-1) %/% m + 1)
      ii <- ifelse(is.null(Amat), 1, Xnew$region.struct[i])
      if("region" %in% by && (!"years" %in% by)){
        which <- which(X$region == region_names[ii])
      }else if((!"region" %in% by) && "years" %in% by){
        which <- which(X$years == year_label[tt]) 
      }else{
        which <- intersect(which(X$years == year_label[tt]), which(X$region == region_names[ii]))
      }
      Xnew[i, covariate.names] <- X[which, covariate.names]
    }
    exdat <- merge(exdat, Xnew, by = c("time.struct", "region.struct"))
   }
  
  if(is.ar){
    exdat$time.slope <- exdat$time.struct 
    center <- N/2 + 1e-5 # avoid exact zero in the lincomb creation
    exdat$time.slope <- (exdat$time.slope  - center) / (sd(1:N))
    exdat$region.slope <- exdat$region.struct
  }


  if(is.null(formula)){
        period.constr <- NULL
        # Tmax <- length(year_label)            
        # if(rw == 2) period.constr <- list(A = matrix(c(rep(1, Tmax)), 1, Tmax), e = 0)

     ## ---------------------------------------------------------
    ## Setup PC prior model
    ## ---------------------------------------------------------
    if(tolower(hyper) == "pc"){
        hyperpc1 <- list(prec = list(prior = "pc.prec", param = c(pc.u , pc.alpha)))
        pc.st.u <- ifelse(is.na(pc.st.u), pc.u, pc.st.u)
        pc.st.alpha <- ifelse(is.na(pc.st.alpha), pc.alpha, pc.st.alpha)
        hyperpc1.interact <- list(prec = list(prior = "pc.prec", param = c(pc.st.u , pc.st.alpha)))
        hyperpc2 <- list(prec = list(prior = "pc.prec", param = c(pc.u , pc.alpha)), 
                         phi = list(prior = 'pc', param = c(pc.u.phi , pc.alpha.phi)))
        hyperar1 = list(prec = list(prior = "pc.prec", param = c(pc.u , pc.alpha)), 
                        theta2 = list(prior = "pc.cor1", param = c(pc.u.cor, pc.alpha.cor)))
        hyperar2 = list(theta2 = list(prior = "pc.cor1", param = c(pc.u.cor, pc.alpha.cor)))
        ## -----------------------
        ## Period + National + PC
        ## ----------------------- 
        if(!is.yearly && is.null(Amat)){

              formula <- logit.est ~ 
                            f(time.struct,model=paste0("rw", rw), constr = TRUE,  extraconstr = period.constr, hyper = hyperpc1) + 
                            f(time.unstruct,model="iid", hyper = hyperpc1) 

        ## -----------------------
        ## Yearly + National + PC
        ## -----------------------
        }else if(is.yearly && is.null(Amat)){

          formula <- logit.est ~ 
              f(time.struct, model = rw.model.pc, diagonal = 1e-6, extraconstr = constr, values = 1:N) + 
              f(time.unstruct,model=iid.model.pc) 
            

        ## -------------------------
        ## Period + Subnational + PC
        ## ------------------------- 
        }else if(!is.yearly && (!is.null(Amat))){

            formula <- logit.est ~ 
                f(time.struct, model=paste0("rw", rw), hyper = hyperpc1, scale.model = TRUE, extraconstr = period.constr)  + 
                f(time.unstruct,model="iid", hyper = hyperpc1) + 
                f(region.struct, graph=Amat,model="bym2", hyper = hyperpc2, scale.model = TRUE, adjust.for.con.comp = TRUE)  

            if(type.st == 1){
                formula <- update(formula, ~. + 
                    f(time.area,model="iid", hyper = hyperpc1.interact))
            }else if(type.st == 2){
                if(!is.ar){
                  formula <- update(formula, ~. + 
                    f(region.int,model="iid", group=time.int,control.group=list(model=paste0("rw", rw), scale.model = TRUE), hyper = hyperpc1.interact))
                }else{
                    formula <- update(formula, ~. + 
                    f(region.int,model="iid", hyper = hyperpc1.interact, group=time.int,control.group=list(model="ar", order = ar, hyper = hyperar2)))
                }

            }else if(type.st == 3){
                formula <- update(formula, ~. + 
                    f(region.int, model="besag", graph = Amat, group=time.int,control.group=list(model="iid"), hyper = hyperpc1.interact, scale.model = TRUE, adjust.for.con.comp = TRUE))
            }else{
             
              if(is.ar){
                formula <- update(formula, ~. + 
                     f(region.int, model="besag", graph = Amat, group=time.int, control.group=list(model="ar", order = ar, hyper = hyperar2), hyper = hyperpc1.interact, scale.model = TRUE, adjust.for.con.comp = TRUE)) 
               }else{
                    # defines type IV explicitly with constraints
                    # Use time.area as index
                    # S blocks each with time 1:T (in this code, 1:N)
                    # UPDATE for connected components:
                    # nc2 sum-to-zero constraints for each of the connected components of size >= 2. Scaled so that the geometric mean of the marginal variances in each connected component of size >= 2 is 1, and modified so that singletons have a standard Normal distribution.
                   inla.rw = utils::getFromNamespace("inla.rw", "INLA")
                   inla.scale.model.bym = utils::getFromNamespace("inla.scale.model.bym", "INLA")
                   inla.bym.constr.internal = utils::getFromNamespace("inla.bym.constr.internal", "INLA")
                   R2 <- inla.rw(N, order = rw, scale.model=TRUE, sparse=TRUE)
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
                   tmp <- rbind(tmp, tmp2)
                   constr.st <- list(A = tmp, e = rep(0, dim(tmp)[1]))
                  formula <- update(formula, ~. + 
                      f(time.area,model="generic0", Cmatrix = R, extraconstr = constr.st, rankdef = N*S -(N - rw)*(S - Q1$rankdef), hyper = hyperpc1.interact))              
               }
            }
          
        ## ------------------------- 
        ## Yearly + Subnational + PC
        ## ------------------------- 
        }else{
              formula <- logit.est ~ 
                  f(time.struct, model = rw.model.pc, diagonal = 1e-6, extraconstr = constr, values = 1:N) +
                  f(time.unstruct,model=iid.model.pc) + 
                  f(region.struct, graph=Amat,model="bym2", hyper = hyperpc2, scale.model = TRUE, adjust.for.con.comp = TRUE) + 
                  f(time.area,model=st.model.pc, diagonal = 1e-6, extraconstr = constr.st, values = 1:(N*S))
        }
       
        if(survey.effect){
          formula <- update(formula, ~. + f(survey.id, model = "iid", hyper = hyperpc1))
        }
    ## ---------------------------------------------------------
    ## Setup Gamma prior model
    ## ---------------------------------------------------------
    }else if(tolower(hyper) == "gamma"){
        ## ------------------- 
        ## Period + National
        ## ------------------- 
        if(!is.yearly && is.null(Amat)){
            formula <- logit.est ~ 
              f(time.struct,model=paste0("rw", rw),param=c(a.rw,b.rw), constr = TRUE)  + 
              f(time.unstruct,model="iid",param=c(a.iid,b.iid)) 
            
        ## ------------------- 
        ## Yearly + National
        ## -------------------   
        }else if(is.yearly && is.null(Amat)){
           formula <- logit.est ~ 
                  f(time.struct, model = rw.model, diagonal = 1e-6, extraconstr = constr, values = 1:N) + 
                  f(time.unstruct,model=iid.model) 
            
        ## ------------------- 
        ## Period + Subnational
        ## ------------------- 
        }else if(!is.yearly && (!is.null(Amat))){
       
            formula <- logit.est ~ 
                  f(time.struct,model=paste0("rw", rw), param=c(a.rw,b.rw), scale.model = TRUE, extraconstr = period.constr)  + 
                  f(time.unstruct,model="iid",param=c(a.iid,b.iid)) + 
                  f(region.struct, graph=Amat,model="besag",param=c(a.icar,b.icar), scale.model = TRUE, adjust.for.con.comp = TRUE) + 
                  f(region.unstruct,model="iid",param=c(a.iid,b.iid)) 
                  
            if(type.st == 1){
                formula <- update(formula, ~. + f(time.area,model="iid", param=c(a.iid,b.iid)))
            }else if(type.st == 2){
                formula <- update(formula, ~. + f(region.int,model="iid", group=time.int,control.group=list(model="rw2", scale.model = TRUE), param=c(a.iid,b.iid)))
            }else if(type.st == 3){
                formula <- update(formula, ~. + f(region.int,model="besag", graph = Amat, group=time.int,control.group=list(model="iid"),param=c(a.iid,b.iid), scale.model = TRUE, adjust.for.con.comp = TRUE))
            }else{
                

              # defines type IV explicitly with constraints
              # Use time.area as index
              # S blocks each with time 1:T (in this code, 1:N)
              # UPDATE!

             inla.rw = utils::getFromNamespace("inla.rw", "INLA")
             inla.scale.model.bym = utils::getFromNamespace("inla.scale.model.bym", "INLA")
             inla.bym.constr.internal = utils::getFromNamespace("inla.bym.constr.internal", "INLA")
             R2 <- inla.rw(N, order = rw, scale.model=TRUE, sparse=TRUE)
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
             tmp <- rbind(tmp, tmp2)
             constr.st <- list(A = tmp, e = rep(0, dim(tmp)[1]))
             

             formula <- update(formula, ~. + 
                    f(time.area,model="generic0", Cmatrix = R, extraconstr = constr.st, rankdef = N*S -(N - rw)*(S - Q1$rankdef), param=c(a.iid,b.iid)))
            }
         
          
        ## ------------------- 
        ## Yearly + Subnational
        ## ------------------- 
        }else{
            formula <- logit.est ~ 
              f(time.struct, model = rw.model, diagonal = 1e-6, extraconstr = constr, values = 1:N) + 
              f(time.unstruct,model=iid.model) + 
              f(region.struct, graph=Amat,model="besag",param=c(a.icar,b.icar), scale.model = TRUE, adjust.for.con.comp = TRUE) + 
              f(region.unstruct,model="iid",param=c(a.iid,b.iid)) + 
              f(time.area,model=st.model, diagonal = 1e-6, extraconstr = constr.st, values = 1:(N*S)) 
        }

        if(survey.effect){
          formula <- update(formula, ~. + f(survey.id, model = "iid", param=c(a.iid,b.iid)))
        }

    }else{
      stop("hyper needs to be either pc or gamma.")
    }

    if(!is.null(X)){
      formula <- as.formula(paste("logit.est~", as.character(formula)[3], "+", paste(covariate.names, collapse = " + ")))
    }

  
}
if(is.main.ar){
     formula <- update(formula, ~. -  
          f(time.struct, model=paste0("rw", rw), hyper = hyperpc1, scale.model = TRUE, extraconstr = period.constr) + 
          f(time.struct, model="ar", order = ar, hyper = hyperar1, extraconstr = period.constr, constr = TRUE) + time.slope)
}



    mod <- formula
    

    
    # borrow the old function from INLA
    inla.uncbind0 <- function (A, name.prefix = "col"){
        if (!is.matrix(A)) 
            return(NULL)
        result = list()
        rownames(A) = NULL
        tmp.result = apply(A, 2, function(x) list(x))
        result = sapply(tmp.result, function(x) c(x))
        return(result)
    }

    
    ## ---------------------------------------------------------
    ## Subnational lincomb for projection
    ## ---------------------------------------------------------
    if(!is.null(Amat)){
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
          if(!is.null(X)){
             if("region" %in% by && (!"years" %in% by)){
                which <- which(exdat$region.struct == j)
              }else if((!"region" %in% by) && "years" %in% by){
                which <- which(exdat$time.struct == i) 
              }else{
                which <- intersect(which(exdat$time.struct == i), which(exdat$region.struct == j))
              }
            sub <- matrix(exdat[which[1], covariate.names], nrow = 1)
            colnames(sub) <- covariate.names
            XX <- inla.uncbind0(sub)
          }else{
            XX <- NULL
          }
          
          object.name <- paste("lc", index, sep = "")
          
          lincombs.info[index, c("District", "Year")] <- c(j,i)
          # BYM model under Gamma prior
          if(hyper == "gamma"){
            if(is.ar && type.st == 1 && (!is.main.ar)){
                tmplin <- list("(Intercept)" = 1,
                                time.area = spacetime,
                                time.struct= time ,
                                time.unstruct= time,
                                # time.slope = (i - center)/sd(1:N),      
                                region.struct = area)
                if(is.main.ar) tmplin <- c(tmplin, time.slope = (i - center)/sd(1:N))
                 assign(object.name, INLA::inla.make.lincomb(c(tmplin, XX)))
            }else if(is.ar && (!is.main.ar)){
                tmplin <- list("(Intercept)" = 1,
                                region.int = area,
                                time.struct= time ,
                                time.unstruct= time,
                                # time.slope = (i - center)/sd(1:N),
                                region.struct = area)
                if(is.main.ar) tmplin <- c(tmplin, time.slope = (i - center)/sd(1:N))
                 assign(object.name, INLA::inla.make.lincomb(c(tmplin, XX)))
            }else  if(is.yearly || type.st%in% c(1, 4)){
                 assign(object.name, INLA::inla.make.lincomb(c(list("(Intercept)" = 1,
                                                    time.area = spacetime,
                                                    time.struct= time ,
                                                    time.unstruct= time,
                                                    region.struct = area,
                                                    region.unstruct = area), 
                                                    XX)))    
            }else{
              newidx <- rep(0, j + (i - 1) * region_count)
              newidx[j + (i - 1) * region_count] <- 1
              assign(object.name, INLA::inla.make.lincomb(c(list("(Intercept)" = 1,
                                                    time.struct= time ,
                                                    time.unstruct= time,
                                                    region.struct = area,
                                                    region.unstruct = area,
                                                    region.int = newidx), 
                                                    XX)))
            }
          # BYM2 model under PC prior
          }else{
            if(is.ar && type.st == 1){
                tmplin <- list("(Intercept)" = 1,
                                time.area = spacetime,
                                time.struct= time ,
                                time.unstruct= time, 
                                # time.slope = (i - center)/sd(1:N),
                                region.struct = area)
                 if(is.main.ar) tmplin <- c(tmplin, time.slope = (i - center)/sd(1:N))
                 assign(object.name, INLA::inla.make.lincomb(c(tmplin, XX)))
            }else if(is.ar){
                tmplin <- list("(Intercept)" = 1,
                                region.int = area,
                                time.struct= time ,
                                time.unstruct= time,
                                # time.slope = (i - center)/sd(1:N),
                                region.struct = area)
                 if(is.main.ar) tmplin <- c(tmplin, time.slope = (i - center)/sd(1:N))
                assign(object.name, INLA::inla.make.lincomb(c(tmplin, XX)))
            }else if(is.yearly || type.st %in% c(1, 4)){
                 assign(object.name, INLA::inla.make.lincomb(c(list("(Intercept)" = 1,
                                      time.area = spacetime,
                                      time.struct= time ,
                                      time.unstruct= time,
                                      region.struct = area), 
                                      XX)))
            }else{
              newidx <- rep(0, j + (i - 1) * region_count)
              newidx[j + (i - 1) * region_count] <- 1
              assign(object.name, INLA::inla.make.lincomb(c(list("(Intercept)" = 1,
                                                    time.struct= time ,
                                                    time.unstruct= time,
                                                    region.struct = area,
                                                    region.int = newidx), 
                                                    XX)))
            }
          }
           
          if(index == 1){
            lincombs.fit <- get(object.name)
            names(lincombs.fit)[index] <- object.name
          }else{
            tmp <- get(object.name)
            lincombs.fit <- c(lincombs.fit, tmp)
            names(lincombs.fit)[index] <- object.name
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
        if(!is.null(X)){
            if("region" %in% by && (!"years" %in% by)){
              which <- which(exdat$region.struct == 0)
            }else if((!"region" %in% by) && "years" %in% by){
              which <- which(exdat$time.struct == i) 
            }else{
              which <- intersect(which(exdat$time.struct == i), which(exdat$region.struct == 0))
            }
          sub <- matrix(exdat[which[1], covariate.names], nrow = 1)
          colnames(sub) <- covariate.names
          XX <- inla.uncbind0(sub)
        }else{
          XX <- NULL
        }
        
        object.name <- paste("lc", index, sep = "")
        
        lincombs.info[index, c("District", "Year")] <- c(0,i)
        if(rw == 1){
          assign(object.name, INLA::inla.make.lincomb(c(list("(Intercept)" = 1,
                                                time.struct= time ,
                                                time.unstruct= time), 
                                                XX)))         
        }else{
          assign(object.name, INLA::inla.make.lincomb(c(list("(Intercept)" = 1,
                                                time.struct= time ,
                                                time.unstruct= time), 
                                                XX)))
        }
        
        if(index == 1){
          lincombs.fit <- get(object.name)
          names(lincombs.fit)[index] <- object.name
        }else{
          lincombs.fit <- c(lincombs.fit, get(object.name))
          names(lincombs.fit)[index] <- object.name
        }
      }
    }
    

    fit <- INLA::inla(mod, family = "gaussian", control.compute = options, data = exdat, control.predictor = list(compute = TRUE), control.family = list(hyper= list(prec = list(initial= log(1), fixed= TRUE ))), scale = exdat$logit.prec, lincomb = lincombs.fit, control.inla = control.inla, verbose = verbose)
    
    return(list(model = mod, fit = fit, Amat = Amat, newdata = exdat, time = seq(0, N - 1), area = seq(0, region_count - 1), time.area = time.area, survey.table = survey.table, a.iid = a.iid, b.iid = b.iid, a.rw = a.rw, b.rw = b.rw, a.rw = a.rw, b.rw = b.rw, a.icar = a.icar, b.icar = b.icar, lincombs.info = lincombs.info, is.yearly = is.yearly, type.st = type.st, year_range = year_range))
  }
}