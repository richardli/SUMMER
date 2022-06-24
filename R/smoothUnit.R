#' Smooth via unit level model
#' 
#' Generates small area estimates by smoothing direct estimates using a unit 
#' level model
#'
#' @param formula an object of class "formula" describing the model to be fitted.
#' @param domain formula specifying variable containing domain labels
#' @param design an object of class "svydesign" containing the data for the model
#' @param family of the response variable, currently supports 'binomial' (default with logit link function) or 'gaussian'.
#' @param Amat Adjacency matrix for the regions. If set to NULL, the IID spatial effect will be used.
#' @param X.pop unit-level covariates data frame. One of the column name needs to match the domain specified, in order to be linked to the data input. Currently only supporting time-invariant domain-level covariates.
#' @param domain.size Domain size data frame. One of the column name needs to match the domain specified, in order to be linked to the data input and there must be a size column containing domain sizes.
#' @param pc.u 	hyperparameter U for the PC prior on precisions.
#' @param pc.alpha hyperparameter alpha for the PC prior on precisions.
#' @param pc.u.phi hyperparameter U for the PC prior on the mixture probability phi in BYM2 model.
#' @param pc.alpha.phi hyperparameter alpha for the PC prior on the mixture probability phi in BYM2 model.
#' @param CI 	the desired posterior credible interval to calculate
#' @param n.sample number of draws from posterior used to compute summaries
#'
#' @return A list with elements
#' \item{direct.est}{direct estimates}
#' \item{model.fit}{fitted INLA object for iid domain effects model}
#' \item{model.est}{smoothed estimates}
#' @importFrom stats model.frame
#' @importFrom stats model.matrix 
#' @importFrom stats model.response 
#' @export
#'
#' @examples
#' \dontrun{
#' library(survey)
#' data(DemoData2)
#' data(DemoMap2)
#' des0 <- svydesign(ids = ~clustid+id, strata = ~strata,
#'                  weights = ~weights, data = DemoData2, nest = T)
#'                  
#' EXAMPLE 1: Continuous response model
#' cts.res <- smoothUnit(formula = tobacco.use ~ 1,
#'                       domain = ~region,
#'                       design = des0, X.pop = DemoData2)
#'                       
#' EXAMPLE 2: Binary response model
#' bin.res <- smoothUnit(formula = tobacco.use ~ 1,
#'                       family = "binomial",
#'                       domain = ~region,
#'                       design = des0, X.pop = DemoData2)
#'}

smoothUnit <- function(formula,
                       domain, 
                       design,
                       family = c("gaussian", "binomial")[1],
                       Amat = NULL,
                       X.pop = NULL,
                       domain.size = NULL,
                       pc.u = 1,
                       pc.alpha = 0.01, 
                       pc.u.phi = 0.5,
                       pc.alpha.phi = 2/3,
                       CI = .95, n.sample = 250) {
  
  if (design$has.strata) {
    warning("This model does not account for stratification yet.")
  }
  if (ncol(design$cluster) > 2) {
    stop("This function is not ready for > 2 stages of sampling")
  }
  # RESPONSE FAMILY
  if (family != "binomial" & family != "gaussian" ) {
    stop(paste0("family = ", family, " not supported"))
  }
  
  # SETUP
  domain.var <- all.vars(domain)[1]
  resp.frm <- as.formula(paste0("~", all.vars(formula)[1]))
  cov.frm <- update(formula, NULL ~ .)
  out <- list()
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
  rownames(dir.est) <- NULL
  colnames(dir.est) <- c("domain", "dir.est", "dir.est.var")
  out$direct.est <- dir.est
  
  # MAKE DATA FRAME
  mf <- model.frame(formula, design$variables)
  resp <- model.response(mf, "numeric")
  
  
  #Xvars <- delete.response(terms(formula))
  
  pop.dat <- data.frame(
    domain = X.pop[[domain.var]]
  )
  mm.pop <- model.matrix(cov.frm, X.pop)
  mod.dat <- data.frame(
    resp = as.vector(resp),
    domain = design$variables[[domain.var]]
  )
  mod.dat <- cbind(mod.dat, model.matrix(cov.frm, design$variables))
  
  domain.table <- data.frame(domain = unique(as.character(pop.dat$domain)))
  if (!is.null(Amat)) {
    domain.table <- data.frame(domain = rownames(Amat))
  }
  domain.table$domain.struct <- seq_len(nrow(domain.table))
  mod.dat$domain.struct <- match(mod.dat$domain, domain.table$domain)
  pop.dat$domain.struct <- match(pop.dat$domain, domain.table$domain)
  # MODEL FORMULA
  ftxt <- paste("resp ~ 1")
  if (length(all.vars(cov.frm)) > 0) {
    ftxt <- 
      paste(ftxt, paste(all.vars(cov.frm), collapse = " + "), sep = " + ")
  }
  
  if (is.null(Amat)) {
    # set priors
    hyperpc.iid.int <- list(
      prec = list(prior = "pc.prec",
                  param = c(pc.u , pc.alpha))
    )
    warning("No spatial information provided, using iid domain effects")
    ftxt <- paste0(ftxt, " + f(domain.struct, model = 'iid', hyper = hyperpc.iid.int)")
  } else {
    hyperpc.bym.int <- list(
      prec = list(prior = "pc.prec", param = c(pc.u , pc.alpha)),  
      phi = list(prior = 'pc', param = c(pc.u.phi , pc.alpha.phi))
    )
    ftxt <- paste0(ftxt, " + f(domain.struct, model = 'bym2', graph=Amat, hyper = hyperpc.bym.int)")
  }
  
  # ADD IID CLUSTER EFFECT
  if(ncol(design$cluster) > 1) {
    mod.dat$cluster <- design$cluster[, 1]
    ftxt <- paste0(ftxt, " + f(cluster, model = 'iid', hyper = hyperpc.iid.int)")
  }
  # FIT MODEL
  mod.frm <- as.formula(ftxt)
  fit <- INLA::inla(mod.frm, family = family, data = mod.dat,
                    control.compute = list(dic = T, mlik = T,
                                           cpo = T, config = TRUE),
                    control.predictor = list(compute = TRUE),
                    lincomb = NULL, quantiles = c((1-CI)/2, 0.5, 1-(1-CI)/2))
  
  
  # GENERATE DRAWS (IF BINOMIAL?)
  sampAll <- INLA::inla.posterior.sample(n = n.sample, result = fit, intern = TRUE)
  
  proj <- cbind(domain.table$domain,
                data.frame(mean=NA, var=NA, median=NA, lower=NA, upper=NA,
                           logit.mean=NA, logit.var=NA, 
                           logit.median=NA, logit.lower=NA, logit.upper=NA))
  fe.idx <- grep(colnames(mm.pop)[1], rownames(sampAll[[1]]$latent))
  fe.idx <- fe.idx:(fe.idx + length(all.vars(cov.frm)))
  re.idx <- grep("domain.struct", x = rownames(sampAll[[1]]$latent))
  summary.sample <- function(x) {
    pop.unit.ests <-
      x$latent[re.idx][pop.dat$domain.struct] + mm.pop %*% x$latent[fe.idx] 
    # pop.unit.ests <- pop.unit.ests + 
    #   rnorm(nrow(mm.pop), sd = sqrt(1 / exp(x$hyperpar[1])))
    if (family == "binomial") {
      pop.unit.ests <- expit(pop.unit.ests)
    }
    area.ests <- aggregate(pop.unit.ests, list(domain = pop.dat$domain.struct), mean)
    return(area.ests[match(1:length(re.idx), area.ests[, 1]), 2])
  }
  est.mat <- do.call(cbind, lapply(sampAll, summary.sample))
  out$model.fit <- fit
  out$model.est <-
    data.frame(domain = unique(pop.dat$domain),
               mean = rowMeans(est.mat),
               median = apply(est.mat, 1,
                              function(x) median(x, na.rm = T)),
               var = apply(est.mat, 1, var),
               lower = apply(est.mat, 1,
                             function(x) quantile(x, (1-CI)/2, na.rm = T)),
               upper = apply(est.mat, 1,
                             function(x) quantile(x, 1-(1-CI)/2, na.rm = T)))
  return(out)
}



