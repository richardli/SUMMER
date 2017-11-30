#' Function to obtain smoothed estimates from fitted INLA object.
#' 
#' 
#' 
#' @param fitted Fitted object from \code{\link{fitINLA}}.
#' @param proj.time String vector of the time periods to project into the future.
#' @param ntime Number of time periods with observation. 
#' @param which.area String of which area to project.
#' @param quantiles Posterior credible interval to return.
#' @param return_raw Logical indicator of whether raw posterior draws or summary statistics is returned. Default to FALSE.
#' @seealso \code{\link{countrySummary}}
#' @return Smoothed model estimates from the provided INLA fit
projINLA_multi <- function(fitted, proj.time, ntime, which.area, quantiles, return_raw = FALSE) {
    quantlabel <- paste0(quantiles, "quant")
    model <- fitted$model
    data <- fitted$newdata
    fit <- fitted$fit
    Amat <- fitted$Amat
    n.region <- dim(Amat)[1]
    a.iid <- fitted$a.iid
    b.iid <- fitted$b.iid
    a.rw1 <- fitted$a.rw1
    b.rw1 <- fitted$b.rw1
    a.rw2 <- fitted$a.rw2
    b.rw2 <- fitted$b.rw2
    a.icar <- fitted$a.icar
    b.icar <- fitted$b.icar

    n.time <- c(ntime:(ntime + length(proj.time)))

    exist <- which(data$time.unstruct %in% n.time)
    if (sum(exist) > 0) {
        data <- data[-exist, ]
    }

    
    # tmp<-data[data$survey==max(data$survey) & data$time.struct==1,] refine this line if the last survey does not cover all
    # regions...
    exdatproj <- data
    for (nextTime in n.time) {
        tmp <- data[match(unique(data$region), data$region), ]
        tmp$time.unstruct <- tmp$time.struct <- tmp$idII <- tmp$groupIII <- tmp$groupIV <- nextTime
        tmp$logit.est <- tmp$logit.prec <- tmp$survey <- tmp$survey.time <- tmp$survey.area <- tmp$survey.time.area <- NA
        tmp$time.area <- (n.region * (nextTime - 1) + 1):(n.region * nextTime)
        tmp$years <- proj.time[nextTime - 5]
        tmp$u5m <- tmp$lower <- tmp$upper <- tmp$var.est <- NA
        exdatproj <- rbind(exdatproj, tmp[, colnames(exdatproj)])
    }
    
    
    smoothed <- rep(NA, max(n.time))
    upper <- rep(NA, max(n.time))
    lower <- rep(NA, max(n.time))
    raw <- matrix(0, max(n.time), 1000)
    
    # - which space-time interaction IDs do we need - #
    spacetimenum <- unique(exdatproj[exdatproj$region == which.area, c("time.area")])
    # - which numeric value this area is coded - #
    spacenum <- unique(exdatproj[exdatproj$region == which.area, c("region_num")])
    
    for (i in 1:max(n.time)) {
        time <- rep(NA, max(n.time))
        time[i] <- 1
        area <- rep(NA, n.region)
        area[spacenum] <- 1
        space.time <- rep(NA, n.region * max(n.time))
        space.time[spacetimenum[i]] <- 1
        
        if (!isTRUE(requireNamespace("INLA", quietly = TRUE))) {
          stop("You need to install the packages 'INLA'. Please run in your R terminal:\n install.packages('INLA', repos='https://www.math.ntnu.no/inla/R/stable')")
        }
        # If INLA is installed, then attach the Namespace (so that all the relevant functions are available)
        if (isTRUE(requireNamespace("INLA", quietly = TRUE))) {
          if (!is.element("INLA", (.packages()))) {
            attachNamespace("INLA")
          }
          lc1 <- INLA::inla.make.lincomb(`(Intercept)` = 1, region.unstruct = area, region.struct = area, time.struct = time, time.unstruct = time, 
                                         time.area = space.time)
          
          mod <- INLA::inla(model, family = "gaussian", data = exdatproj, lincomb = lc1, control.predictor = list(compute = TRUE), 
                            control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))), scale = exdatproj$logit.prec, quantiles = quantiles)
          smoothed[i] <- expit(mod$summary.lincomb.derived[, quantlabel[2]])
          upper[i] <- expit(mod$summary.lincomb.derived[, quantlabel[3]])
          lower[i] <- expit(mod$summary.lincomb.derived[, quantlabel[1]])
          raw[i, ] <- INLA::inla.rmarginal(1000, mod$marginals.lincomb.derived[[1]])
        }
        #cat(".")
    }
    #cat("\n")
    if (return_raw) {
        return(raw)
    } else {
        return(data.frame(cbind(smoothed, lower, upper)))
    }
}
