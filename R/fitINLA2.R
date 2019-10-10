#' Fit cluster-level space-time smoothing models to mortality rates 
#' 
#' 
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
#' @param age.groups a character vector of age groups in increasing order.
#' @param age.n number of months in each age groups in the same order.
#' @param age.rw.group vector indicating grouping of the ages groups. For example, if each age group is assigned a different random walk component, then set age.rw.group to c(1:length(age.groups)); if all age groups share the same random walk component, then set age.rw.group to a rep(1, length(age.groups)). The default for 6 age groups is c(1,2,3,3,3,3), which assigns a separate random walk to the first two groups and a common random walk for the rest of the age groups. The vector should contain values starting from 1.
#' @param family family of the model. This can be either binomial (with logistic normal prior) or betabiniomial.
#' @param Amat Adjacency matrix for the regions
#' @param geo Geo file
#' @param bias.adj the ratio of unadjusted mortality rates or age-group-specific hazards to the true rates or hazards. It needs to be a data frame that can be merged to thee outcome, i.e., with the same column names for time periods (for national adjustment), or time periods and region (for subnational adjustment). The column specifying the adjustment ratio should be named "ratio".
#' @param bias.adj.by vector of the column names specifying how to merge the bias adjustment to the count data. For example, if bias adjustment factor is provided in bias.adj for each region and time, then bias.adj.by should be `c("region", "time")`.
#' @param formula INLA formula.  See vignette for example of using customized formula.
#' @param year_label string vector of year names
#' @param priors priors from \code{\link{simhyper}}
#' @param rw Take values 1 or 2, indicating the order of random walk.
#' @param type.st type for space-time interaction
#' @param hyper which hyperpriors to use. Default to be using the PC prior ("pc"). 
#' @param pc.u hyperparameter U for the PC prior on precisions.
#' @param pc.alpha hyperparameter alpha for the PC prior on precisions.
#' @param pc.u.phi hyperparameter U for the PC prior on the mixture probability phi in BYM2 model.
#' @param pc.alpha.phi hyperparameter alpha for the PC prior on the mixture probability phi in BYM2 model.
#' @param a.iid hyperparameter for i.i.d random effects.
#' @param b.iid hyperparameter for i.i.d random effects.
#' @param a.rw hyperparameter for RW 1 or 2 random effects.
#' @param b.rw hyperparameter for RW 1 or 2random effects.
#' @param a.icar hyperparameter for ICAR random effects.
#' @param b.icar hyperparameter for ICAR random effects.
#' @param options list of options to be passed to control.compute() in the inla() function.
#' @param verbose logical indicator to print out detailed inla() intermediate steps.
#' @seealso \code{\link{getDirect}}
#' @import Matrix
#' @importFrom stats dgamma
#' @importFrom Matrix Diagonal 
#' @return INLA model fit using the provided formula, country summary data, and geographic data
#' @examples
#' message("Please check the package vignette on binomial models.")
#' 
#' @export
#' 
#' 

fitINLA2 <- function(data, family = c("betabinomial", "binomial")[1], age.groups = c("0", "1-11", "12-23", "24-35", "36-47", "48-59"), age.n = c(1,11,12,12,12,12), age.rw.group = 1:6, Amat, geo, bias.adj = NULL, bias.adj.by = NULL, formula = NULL, rw = 2, year_label, priors = NULL, type.st = 1, hyper = c("pc", "gamma")[1], pc.u = 1, pc.alpha = 0.01, pc.u.phi = 0.5, pc.alpha.phi = 2/3, a.iid = NULL, b.iid = NULL, a.rw = NULL, b.rw = NULL, a.icar = NULL, b.icar = NULL, options = list(config = TRUE), verbose = FALSE){


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

  if(is.null(age.groups)){
    age.n <- 1
    age.rw.group <- 1
  }

    tau = exp(10)

   
    
    ## ---------------------------------------------------------
    ## Common Setup
    ## --------------------------------------------------------- 
    if(!is.null(geo)){
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
    # if (na.rm) {
    #   na.count <- apply(data, 1, function(x) {
    #     length(which(is.na(x)))
    #   })
    #   to_remove <- which(na.count == 6)
    #   if (length(to_remove) > 0) 
    #     data <- data[-to_remove, ]
    # }
    # #################################################################### get the list of region and numeric index in one data frame
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
    dat$region.struct <- dat$region.unstruct <- dat$region.int <- dat$region_number
    
    #################################################################
    ## get the list of region and numeric index in one data frame    
    n <- 0
    N <- nn <- length(year_label)
    years <- data.frame(year = year_label, year_number = seq(1, N))      
  
    
    # -- creating IDs for the temporal REs -- #
    dat$time.unstruct <- dat$time.struct <- dat$time.int <- years[match(dat$years, years[, 1]), 2]
    

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
    
    x <- expand.grid(1:N, 1:region_count)
    
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
      newdata <- merge(newdata, time.area, 
        by = c("region_number", "time.unstruct"))
    }else{
      newdata$time.area <- NA
    }
    
    
    ########################## Model Selection ######
    
    # -- subset of not missing and not direct estimate of 0 -- #
    exdat <- newdata
    clusters <- unique(exdat$cluster)
    exdat$cluster.id <- match(exdat$cluster, clusters)
    exdat$nugget.id <- exdat$cluster.id
    # cluster.time <- expand.grid(cluster = clusters, time = 1:N)
    # cluster.time$nugget.id <- 1:dim(cluster.time)[1]
    # exdat <- merge(exdat, cluster.time, by.x = c("cluster", "time.struct"), by.y = c("cluster", "time"))
    # # exdat$nugget.id <- 1:dim(exdat)[1]

  replicate.rw <- length(unique(age.rw.group)) > 1


  if(is.null(formula)){
        period.constr <- NULL
        # Tmax <- length(year_label)            
        # if(rw == 2) period.constr <- list(A = matrix(c(rep(1, Tmax)), 1, Tmax), e = 0)
        if(rw %in% c(1, 2) == FALSE) stop("Random walk only support rw = 1 or 2.")
   
     ## ---------------------------------------------------------
    ## Setup PC prior model
    ## ---------------------------------------------------------
    if(tolower(hyper) == "pc"){
        hyperpc1 <- list(prec = list(prior = "pc.prec", param = c(pc.u , pc.alpha)))
        hyperpc2 <- list(prec = list(prior = "pc.prec", param = c(pc.u , pc.alpha)), 
                         phi = list(prior = 'pc', param = c(pc.u.phi , pc.alpha.phi)))
        
        ## -----------------------
        ## Period + National + PC
        ## ----------------------- 
        if(is.null(geo)){
          if(replicate.rw){
              formula <- Y ~
                            f(time.struct,model=paste0("rw", rw), constr = TRUE,  extraconstr = period.constr, hyper = hyperpc1, replicate =  age.rep.idx) + 
                            f(time.unstruct,model="iid", hyper = hyperpc1) 
          }else{
              formula <- Y ~
                            f(time.struct,model=paste0("rw", rw), constr = TRUE,  extraconstr = period.constr, hyper = hyperpc1) + 
                            f(time.unstruct,model="iid", hyper = hyperpc1)             
          }
        ## -------------------------
        ## Period + Subnational + PC
        ## ------------------------- 
        }else if(!is.null(geo)){

             if(replicate.rw){
              formula <- Y ~ 
                f(time.struct, model=paste0("rw", rw), hyper = hyperpc1, scale.model = TRUE, extraconstr = period.constr, values = 1:N, replicate =  age.rep.idx) 
             }else{
               formula <- Y ~ 
                f(time.struct, model=paste0("rw", rw), hyper = hyperpc1, scale.model = TRUE, extraconstr = period.constr, values = 1:N) 
             }
             formula <- update(formula, ~. + 
                  f(time.unstruct,model="iid", hyper = hyperpc1, values = 1:N) + 
                  f(region.struct, graph=Amat,model="bym2", hyper = hyperpc2, scale.model = TRUE, adjust.for.con.comp = TRUE)) 

            if(type.st == 1){
                formula <- update(formula, ~. + 
                    f(time.area,model="iid", hyper = hyperpc1))
            }else if(type.st == 2){
                formula <- update(formula, ~. + 
                    f(region.int,model="iid", group=time.int,control.group=list(model=paste0("rw", rw), scale.model = TRUE), hyper = hyperpc1))
            }else if(type.st == 3){
                formula <- update(formula, ~. + 
                    f(region.int, model="besag", graph = Amat, group=time.int,control.group=list(model="iid"), hyper = hyperpc1, scale.model = TRUE, adjust.for.con.comp = TRUE))
            }else{
                

              # defines type IV explicitly with constraints
              # Use time.area as index
              # S blocks each with time 1:T (in this code, 1:N)
              # UPDATE!

             inla.rw = utils::getFromNamespace("inla.rw", "INLA")
             R2 <- inla.rw(N, order = rw, scale.model=TRUE, sparse=TRUE)
             R4 = Amat
             diag(R4) <- 0
             diag <- apply(R4, 1, sum)
             R4[R4 != 0] <- -1
             diag(R4) <- diag
             R4 <- INLA::inla.scale.model(R4, constr = list(A=matrix(1,1,dim(R4)[1]), e=0))
             R <- R4 %x% R2
             tmp <- matrix(0, S, N * S)
             for(i in 1:S){
               tmp[i, ((i-1)*N + 1) : (i*N)] <- 1
             }
             tmp2 <- matrix(0, N, N * S)
             for(i in 1:N){
                tmp2[i , which((1:(N*S)) %% N == i-1)] <- 1
             }
             tmp <- rbind(tmp, tmp2)
             constr.st <- list(A = tmp, e = rep(0, dim(tmp)[1]))
             

             formula <- update(formula, ~. + 
                    f(time.area,model="generic0", Cmatrix = R, extraconstr = constr.st, rankdef = N*S -(N - rw)*(S - 1), hyper = hyperpc1))
            }
          
        ## ------------------------- 
        ## Yearly + Subnational + PC
        ## ------------------------- 
        }

    ## ---------------------------------------------------------
    ## Setup Gamma prior model
    ## ---------------------------------------------------------
    }else if(tolower(hyper) == "gamma"){
        ## ------------------- 
        ## Period + National
        ## ------------------- 
        if(is.null(geo)){
            if(replicate.rw){
              formula <- Y ~
                f(time.struct,model=paste0("rw", rw),param=c(a.rw,b.rw), constr = TRUE, extraconstr = period.constr, hyper = hyperpc1, replicate =  age.rep.idx)  + 
                f(time.unstruct,model="iid",param=c(a.iid,b.iid)) 

            }else{
              formula <- Y ~
                f(time.struct,model=paste0("rw", rw),param=c(a.rw,b.rw), constr = TRUE)  + 
                f(time.unstruct,model="iid",param=c(a.iid,b.iid))               
            }
            
        }else{
          if(replicate.rw){
            formula <- Y ~
                  f(time.struct,model=paste0("rw", rw), param=c(a.rw,b.rw), scale.model = TRUE, extraconstr = period.constr)  
          }else{
            formula <- Y ~
                  f(time.struct,model=paste0("rw", rw), param=c(a.rw,b.rw), scale.model = TRUE, extraconstr = period.constr)  
          }
          formula <- update(formula,  ~. +
                  f(time.unstruct,model="iid",param=c(a.iid,b.iid)) + 
                  f(region.struct, graph=Amat,model="besag",param=c(a.icar,b.icar), scale.model = TRUE, adjust.for.con.comp = TRUE) + 
                  f(region.unstruct,model="iid",param=c(a.iid,b.iid))) 
                  
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
                 # UPDATE for connected components:
                 # nc2 sum-to-zero constraints for each of the connected components of size >= 2. Scaled so that the geometric mean of the marginal variances in each connected component of size >= 2 is 1, and modified so that singletons have a standard Normal distribution.

                inla.rw = utils::getFromNamespace("inla.rw", "INLA")
                R2 <- inla.rw(N, order = rw, scale.model=TRUE, sparse=TRUE)
                R4 = Amat
                diag(R4) <- 0
                diag <- apply(R4, 1, sum)
                R4[R4 != 0] <- -1
                diag(R4) <- diag
                R4 <- INLA::inla.scale.model(R4, constr = list(A=matrix(1,1,dim(R4)[1]), e=0))
                R <- R4 %x% R2
                tmp <- matrix(0, S, N * S)
                for(i in 1:S){
                  tmp[i, ((i-1)*N + 1) : (i*N)] <- 1
                }
                tmp2 <- matrix(0, N, N * S)
                for(i in 1:N){
                   tmp2[i , which((1:(N*S)) %% N == i-1)] <- 1
                }
                tmp <- rbind(tmp, tmp2)
                constr.st <- list(A = tmp, e = rep(0, dim(tmp)[1]))
                

                formula <- update(formula, ~. + 
                    f(time.area,model="generic0", Cmatrix = R, extraconstr = constr.st, rankdef = N*S -(N - rw)*(S - 1), param=c(a.iid,b.iid)))
            }
         
          
        ## ------------------- 
        ## Yearly + Subnational
        ## ------------------- 
        }
    }else{
      stop("hyper needs to be either pc or gamma.")
    }

    if(family == "binomial"){
      if(tolower(hyper) == "gamma"){
          formula <- update(formula, ~.+ f(nugget.id,model="iid",model="iid", param=c(a.iid,b.iid)))
      }else if(tolower(hyper) == "pc"){
          formula <- update(formula, ~.+ f(nugget.id,model="iid", hyper = hyperpc1))
      }else{
          stop("hyper needs to be either pc or gamma.")
      }
    }
    formula <- update(formula, ~. -1 + age + strata)
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


  ## add yearly observations with NA outcome and 1 trial, does not contribute to likelihood
  total <- NA
  exdat <- subset(exdat, total != 0)
  
  # tmp<-exdat[1:N, ]
  # tmp$time.unstruct<-tmp$time.struct<-tmp$time.int <- 1:N
  # tmp$years<-years[1:N, 1]
  # tmp$total <- 1
  # tmp$Y <- NA
  # exdat<-rbind(exdat,tmp)   
  for(i in 1:N){
      tmp<-exdat[match(unique(exdat$region), exdat$region), ]
      tmp <- tmp[]
      tmp$time.unstruct<-tmp$time.struct<- tmp$time.int <- i
      tmp <- tmp[, colnames(tmp) != "time.area"]
      tmp <- merge(tmp, time.area, by = c("region_number", "time.unstruct"))
      tmp$years<-years[i, 1]
      tmp$total <- 1
      tmp$Y <- NA
      exdat<-rbind(exdat,tmp)   
    }


  exdat$strata <- factor(exdat$strata)
  if(!is.null(age.groups)){
      exdat$age <- factor(exdat$age, levels = age.groups)    
  }else{
     formula <- update(formula, ~.-age)   
  }
  exdat$age.idx <- match(exdat$age, age.groups)
  exdat$age.rep.idx <- age.rw.group[exdat$age.idx]
     
  fit <- INLA::inla(formula, family = family, control.compute = options, data = exdat, control.predictor = list(compute = FALSE), Ntrials = exdat$total, lincomb = NULL, control.inla = list(int.strategy = "auto"), verbose = verbose)

 # find the name for baseline strata
 levels <- grep("strata", rownames(fit$summary.fixed))   
 levels <- gsub("strata", "", rownames(fit$summary.fixed)[levels])
 strata.all <- as.character(unique(exdat$strata))
 strata.base <- strata.all[strata.all%in%levels == FALSE]

  return(list(model = formula, fit = fit, family= family, Amat = Amat, newdata = exdat, time = seq(0, N - 1), area = seq(0, region_count - 1), survey.time = survey.time, survey.area = survey.area, time.area = time.area, survey.time.area = survey.time.area, a.iid = a.iid, b.iid = b.iid, a.rw = a.rw, b.rw = b.rw, a.rw = a.rw, b.rw = b.rw, a.icar = a.icar, b.icar = b.icar, is.yearly = FALSE, type.st = type.st, age.groups = age.groups, age.n = age.n, age.rw.group = age.rw.group, strata.base = strata.base))
    
  }
}