#' Smooth via area level model
#' 
#' Generates small area estimates  by smoothing direct estimates using an area 
#' level model
#'
#' @param formula an object of class 'formula' describing the model to be fitted. 
#'   If direct.est is specified, the right hand side of the formula is not necessary.
#' @param domain formula specifying variable containing domain labels
#' @param design an object of class "svydesign" containing the data for the model
#' @param responseType of the response variable, currently supports 'binary' (default with logit link function) or 'gaussian'.
#' @param Amat adjacency matrix for the regions. If set to NULL, the IID spatial effect will be used.
#' @param direct.est data frame of direct estimates, with first column containing domain, second column containing direct estimate, and third column containing variance of direct estimate. 
#' @param X.area areal covariates data frame. One of the column names needs to match 'domain', in order to be linked to the data input. Currently only supporting time-invariant domain-level covariates.
#' @param domain.size domain size data frame. One of the column names needs to match 'domain' in order to be linked to the data input and there must be a size column containing domain sizes.
#' @param pc.u 	hyperparameter U for the PC prior on precisions.
#' @param pc.alpha hyperparameter alpha for the PC prior on precisions.
#' @param pc.u.phi hyperparameter U for the PC prior on the mixture probability phi in BYM2 model.
#' @param pc.alpha.phi hyperparameter alpha for the PC prior on the mixture probability phi in BYM2 model.
#' @param CI 	the desired posterior credible interval to calculate
#' @param n.sample number of draws from posterior used to compute summaries
#' @param var.tol tolerance parameter; if variance of an area's direct estimator is below this value, that direct estimator is dropped from model
#'
#' @return A list with elements
#' \item{direct.est}{direct estimates}
#' \item{s.dir.iid.fit}{fitted INLA object for iid domain effects model}
#' \item{s.dir.iid.est}{non-spatial smoothed estimates}
#' \item{s.dir.sp.fit}{fitted INLA object for spatial domain effects model}
#' \item{s.dir.sp.est}{spatially smoothed estimates (if adjacency matrix provided)}
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' library(survey)
#' data(DemoData2)
#' data(DemoMap2)
#' des0 <- svydesign(ids = ~clustid+id, strata = ~strata,
#'                   weights = ~weights, data = DemoData2, nest = T)
#' Xmat <- aggregate(age~region, data = DemoData2, FUN = mean)
#' 
#' EXAMPLE 1: Continuous response model
#' cts.res <- smoothArea(tobacco.use ~ 1, domain = ~region,
#'                       Amat = DemoMap2$Amat, design = des0,
#'                       pc.u = 1,
#'                       pc.alpha = 0.01,
#'                       pc.u.phi = 0.5,
#'                       pc.alpha.phi = 2/3)
#'                       
#' EXAMPLE 2: Including area level covariates
#' cts.cov.res <- smoothArea(tobacco.use ~ age, domain = ~region,
#'                           Amat = DemoMap2$Amat, design = des0,
#'                           X.area = Xmat,
#'                           pc.u = 1,
#'                           pc.alpha = 0.01,
#'                           pc.u.phi = 0.5,
#'                           pc.alpha.phi = 2/3)
#'                           
#' EXAMPLE 3: Binary response model
#' bin.res <- smoothArea(tobacco.use ~ 1, domain = ~region,
#'                       responseType = "binary",
#'                       Amat = DemoMap2$Amat, design = des0,
#'                       pc.u = 1,
#'                       pc.alpha = 0.01,
#'                       pc.u.phi = 0.5,
#'                       pc.alpha.phi = 2/3)
#'                       
#' EXAMPLE 4: Including area level covariates in binary response model
#' bin.cov.res <- smoothArea(tobacco.use ~ age, domain = ~region,
#'                           responseType = "binary",
#'                           Amat = DemoMap2$Amat, design = des0,
#'                           X.area = Xmat,
#'                           pc.u = 1,
#'                           pc.alpha = 0.01,
#'                           pc.u.phi = 0.5,
#'                           pc.alpha.phi = 2/3)
#'}
#'
smoothArea <- function(formula,
                       domain,
                       design = NULL,
                       responseType = c("gaussian", "binary")[1],
                       Amat = NULL,
                       direct.est = NULL,
                       X.area = NULL,
                       domain.size = NULL,
                       pc.u = 1,
                       pc.alpha = 0.01, 
                       pc.u.phi = 0.5,
                       pc.alpha.phi = 2/3,
                       CI = .95, 
                       n.sample = 250,
                       var.tol = 1e-10) {
  

  # SETUP
  domain.var <- all.vars(domain)[1]
  resp.frm <- as.formula(paste0("~", all.vars(update(formula, . ~ NULL))[1]))
  cov.frm <- update(formula, NULL ~ .)
  
  out <- list()

  if (is.null(direct.est)) {
    # calculate direct estimates
    if (!is.null(domain.size)) {
      dir.est <- survey::svyby(resp.frm, domain, design = design, survey::svytotal, na.rm = T)
      dir.est$domain.size <- domain.size$size[match(dir.est[, 1], domain.size[[domain.var]])]
      dir.est[, 2] = dir.est[, 2] / dir.est$domain.size
      dir.est[, 3] = (dir.est[, 3] / dir.est$domain.size) ^ 2
      dir.est <- dir.est[, 1:3]
    } else {
      dir.est <- survey::svyby(resp.frm, domain, design = design, survey::svymean, na.rm = T)
      dir.est[, 3] <- dir.est[, 3] ^ 2
    }
    
    
  } else {
    # assumes structure of direct.est!!!!
    dir.est <- direct.est
  }
  rownames(dir.est) <- NULL
  colnames(dir.est) <- c("domain", "dir.est", "dir.est.var")
  out$direct.est <- dir.est
  
  # remove any areas with zero sampling variances
  mod.dat <- dir.est
  mod.dat$dir.est <- ifelse(mod.dat$dir.est.var > var.tol, mod.dat$dir.est, NA)
  mod.dat$dir.est.var <- ifelse(mod.dat$dir.est.var > var.tol, mod.dat$dir.est.var, NA)
  mod.dat$dir.est.prec <- 1 / mod.dat$dir.est.var
  mod.dat$domain.id <- 1:nrow(mod.dat)
  if(!is.null(X.area)) {
    mod.X.area <- X.area
    mod.X.area$domain <- X.area[[domain.var]]
    mod.X.area <- mod.X.area[, names(mod.X.area) != domain.var]
    mod.dat <- merge(mod.dat, mod.X.area,  by = "domain")
  }
  mod.dat <- mod.dat[match(1:nrow(mod.dat), mod.dat$domain.id), ]
  mm.area <- stats::model.matrix(cov.frm, mod.dat)
  
  h <- function(x) x
  h.inv <- function(x) x
  if (responseType == "binary") {
    h <- function(x) logit(x)
    h.inv <- function(x) expit(x)
    mod.dat$dir.est.prec <- 
      (mod.dat$dir.est^2*(1-mod.dat$dir.est)^2) / mod.dat$dir.est.var
    mod.dat$dir.est <- h(mod.dat$dir.est)

  }
  mod.dat$scale <- mod.dat$dir.est.prec
  
  s.dir.ftxt <- 
    paste("dir.est ~ 1")
  if (length(all.vars(cov.frm)) > 0) {
    s.dir.ftxt <- 
      paste(s.dir.ftxt, paste(all.vars(cov.frm), collapse = " + "), sep = " + ")
  }
  s.dir.iid.ftxt <- 
    paste0(s.dir.ftxt, " + f(domain.id, model = 'iid', hyper = hyperpc.iid.int)")
  # set priors
  hyperpc.iid.int <- list(
    prec = list(prior = "pc.prec",
                param = c(pc.u , pc.alpha))
  )
  # SMOOTHED DIRECT w/ IID RE
  s.dir.iid.fit <-
    INLA::inla(as.formula(s.dir.iid.ftxt),
               family = "gaussian", data = mod.dat, 
               scale = mod.dat$dir.est.prec,
               control.family = 
                 list(hyper = list(prec = list(initial= log(1), fixed= TRUE))),
               control.predictor = list(compute = TRUE),
               control.compute=list(config = TRUE))
  sample.id <-c(list(domain.id = 1:nrow(mod.dat)),
                as.list(setNames(rep(1, ncol(mm.area)), colnames(mm.area))))
  s.dir.iid.sample <-
    INLA::inla.posterior.sample(n = n.sample, s.dir.iid.fit, sample.id)
  fe.idx <- grep(colnames(mm.area)[1], rownames(s.dir.iid.sample[[1]]$latent))
  fe.idx <- fe.idx:(fe.idx + ncol(mm.area) - 1)
  s.dir.iid.mat <-
    do.call(cbind, lapply(s.dir.iid.sample,
                          function(x) x$latent[1:nrow(mod.dat)] +
                            mm.area %*% x$latent[fe.idx]))

  s.dir.iid.mat <- h.inv(s.dir.iid.mat)
  out$s.dir.iid.fit <- s.dir.iid.fit

  out$s.dir.iid.est <- 
    data.frame(domain = mod.dat$domain,
               mean = rowMeans(s.dir.iid.mat),
               median = apply(s.dir.iid.mat, 1,
                              function(x) median(x, na.rm = T)),
               var = apply(s.dir.iid.mat, 1, var),
               lower = apply(s.dir.iid.mat, 1,
                             function(x) quantile(x, (1-CI)/2, na.rm = T)),
               upper = apply(s.dir.iid.mat, 1,
                             function(x) quantile(x, 1-(1-CI)/2, na.rm = T)),
               method = paste0("Smoothed Direct: IID"))
  
  # SMOOTHED DIRECT w/ BYM2 RE
  if (!is.null(Amat)) {
    mod.dat$domain.id <- match(mod.dat$domain, rownames(Amat))
    mod.dat <- mod.dat[match(1:nrow(mod.dat), mod.dat$domain.id), ]
    mm.area <- model.matrix(cov.frm, mod.dat)
    hyperpc.bym.int <- list(
      prec = list(prior = "pc.prec", param = c(pc.u , pc.alpha)),  
      phi = list(prior = 'pc', param = c(pc.u.phi , pc.alpha.phi))
    )
    s.dir.sp.ftxt <- 
      paste0(s.dir.ftxt,
             "+ f(domain.id, model = 'bym2', graph = Amat,", 
             "hyper = hyperpc.bym.int, scale.model = TRUE)")
    s.dir.sp.fit <-
      INLA::inla(as.formula(s.dir.sp.ftxt),
                 family = "gaussian", data = mod.dat, 
                 scale = mod.dat$dir.est.prec,
                 control.family = 
                   list(hyper = list(prec = list(initial= log(1), fixed= TRUE))),
                 control.predictor = list(compute = TRUE),
                 control.compute=list(config = TRUE))
    sample.id <-c(list(domain.id = 1:nrow(Amat)),
                  as.list(setNames(rep(1, ncol(mm.area)), colnames(mm.area))))
    s.dir.sp.sample <-
      INLA::inla.posterior.sample(n = n.sample, s.dir.sp.fit, sample.id)
    fe.idx <- grep(colnames(mm.area)[1], rownames(s.dir.sp.sample[[1]]$latent))
    fe.idx <- fe.idx:(fe.idx + ncol(mm.area) - 1)
    s.dir.sp.mat <-
      do.call(cbind, lapply(s.dir.sp.sample,
                            function(x) x$latent[1:nrow(Amat)] +
                              mm.area %*% x$latent[fe.idx]))
    s.dir.sp.mat <- h.inv(s.dir.sp.mat)
    out$s.dir.sp.fit <- s.dir.sp.fit
    out$s.dir.sp.est <- 
      data.frame(domain = mod.dat$domain,
                 mean = rowMeans(s.dir.sp.mat),
                 median = apply(s.dir.sp.mat, 1,
                                function(x) median(x, na.rm = T)),
                 var = apply(s.dir.sp.mat, 1, var),
                 lower = apply(s.dir.sp.mat, 1,
                               function(x) quantile(x, (1-CI)/2, na.rm = T)),
                 upper = apply(s.dir.sp.mat, 1,
                               function(x) quantile(x, 1-(1-CI)/2, na.rm = T)),
                 method = paste0("Smoothed Direct: BYM2"))
    out$s.dir.sp.est <- 
      out$s.dir.sp.est[match(out$direct.est$domain, out$s.dir.sp.est$domain),]
    
  }
  return(out)
}
