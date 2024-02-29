#' Smooth via basic unit level model
#' 
#' Generates small area estimates by smoothing direct estimates using a basic 
#' unit level model
#'
#' @param formula An object of class 'formula' describing the model to be fitted. 
#' @param domain One-sided formula specifying factors containing domain labels
#' @param design An object of class "svydesign" containing the data for the model 
#' @param family of the response variable, currently supports 'binomial' (default with logit link function) or 'gaussian'.
#' @param X.pop Data frame of population unit-level covariates. One of the column name needs to match the domain specified, in order to be linked to the data input. Currently only supporting time-invariant covariates.
#' @param adj.mat Adjacency matrix with rownames matching the domain labels. If set to NULL, the IID spatial effect will be used.
#' @param domain.size Data frame of domain sizes. One of the column names needs to match the name of the domain variable, in order to be linked to the data input and there must be a column names 'size' containing domain sizes. The default option is no transformation, but logit and log are implemented.
#' @param pc.u 	Hyperparameter U for the PC prior on precisions. See the INLA documentation for more details on the parameterization.
#' @param pc.alpha Hyperparameter alpha for the PC prior on precisions.
#' @param pc.u.phi Hyperparameter U for the PC prior on the mixture probability phi in BYM2 model.
#' @param pc.alpha.phi Hyperparameter alpha for the PC prior on the mixture probability phi in BYM2 model.
#' @param level The specified level for the posterior credible intervals
#' @param n.sample Number of draws from posterior used to compute summaries
#' @param return.samples If TRUE, return matrix of posterior samples of area level quantities
#' @param X.pop.weights Optional vector of weights to use when aggregating unit level predictions
#'
#' @return A svysae object
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' data(DemoData2)
#' data(DemoMap2)
#' library(survey)
#' des0 <- svydesign(ids = ~clustid+id, strata = ~strata,
#'                  weights = ~weights, data = DemoData2, nest = TRUE)
#'                  
#' # EXAMPLE 1: Continuous response model
#' cts.res <- smoothUnit(formula = tobacco.use ~ 1,
#'                       domain = ~region,
#'                       design = des0, X.pop = DemoData2)
#'                       
#' # EXAMPLE 2: Binary response model
#' bin.res <- smoothUnit(formula = tobacco.use ~ 1,
#'                       family = "binomial",
#'                       domain = ~region,
#'                       design = des0, X.pop = DemoData2)
#'}
smoothUnit <- function(formula,
                       domain, 
                       design,
                       family = c("gaussian", "binomial")[1],
                       X.pop = NULL,
                       adj.mat = NULL,
                       domain.size = NULL,
                       pc.u = 1,
                       pc.alpha = 0.01, 
                       pc.u.phi = 0.5,
                       pc.alpha.phi = 2/3,
                       level = .95, n.sample = 250,
                       return.samples = F,
                       X.pop.weights = NULL) {
  
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
  
  # SETUP ----------------------------------------------------------------------
  
  # set up return value
  out <- list()
  attr(out, "inla.fitted") <- c()
  
  # get domain variable
  domain.var <- all.vars(domain)
  if (length(domain.var) != 1) {
    stop("Domain labels must be contained in a single column. Spatio-temporal not supported yet.")
  } 
  
  # set up formulas
  resp.frm <- as.formula(paste0("~", all.vars(formula)[1]))
  cov.frm <- update(formula, NULL ~ .)
  
  # DIRECT ESTIMATES -----------------------------------------------------------
  # compute direct estimates (Hajek estimates if domain size unknown)
  if (!is.null(domain.size)) {
    direct.est <- svyby(resp.frm, domain, design = design, svytotal, na.rm = TRUE)
    direct.est$domain.size <- 
      domain.size$size[match(direct.est[, 1], domain.size[[domain.var]])]
    direct.est[, 2] = direct.est[, 2] / direct.est$domain.size
    direct.est[, 3] = (direct.est[, 3] / direct.est$domain.size) ^ 2
    direct.est <- direct.est[, 1:3]
  } else {
    direct.est <- svyby(resp.frm, domain, design = design, svymean, na.rm = TRUE)
    direct.est[, 3] <- direct.est[, 3] ^ 2
  }
  rownames(direct.est) <- NULL
  colnames(direct.est) <- c("domain", "mean", "var")
  direct.est <- 
    data.frame(domain = direct.est$domain,
               mean = direct.est$mean,
               median = direct.est$mean,
               var = direct.est$var,
               lower = direct.est$mean + qnorm((1-level)/2) * sqrt(direct.est$var),
               upper = direct.est$mean + qnorm(1 - (1-level)/2) * sqrt(direct.est$var),
               method = paste0("Direct"))
  
  
  if (!is.null(adj.mat) & !setequal(X.pop[[domain.var]], rownames(adj.mat))) {
    stop("Domains in X.pop do not match domains in adj.mat.")
  }
  if (any(is.na(match(X.pop[[domain.var]], direct.est$domain)))) {
    warning(cat("There are domains in X.pop not in design/direct estimates.",
                "\nGenerating estimates for all domains in X.pop."))
  }
  # if no adjacency matrix matches, take domain names from X.pop
  domain.table <- data.frame(domain = unique(as.character(X.pop[[domain.var]])))
  # if adjacency matrix provided, take domain names from row names
  if (!is.null(adj.mat)) {
    domain.table <- data.frame(domain = rownames(adj.mat))
  }
  
  direct.est <- 
    merge(direct.est, data.frame(domain = domain.table$domain), 
          by = "domain", all.y = TRUE)
  direct.est$method = "Direct"
  
  out$direct.est <- direct.est
  attr(out, "domain.names") <- sort(direct.est$domain)
  attr(out, "method.names") <- c("direct.est")
  
  # UNIT LEVEL MODEL -----------------------------------------------------------
  mf <- model.frame(formula, design$variables)
  resp <- model.response(mf, "numeric")
  
  pop.dat <- data.frame(
    domain = X.pop[[domain.var]]
  )
  mm.pop <- model.matrix(cov.frm, X.pop)
  mod.dat <- data.frame(
    resp = as.vector(resp),
    domain = design$variables[[domain.var]]
  )
  mod.dat <- cbind(mod.dat, model.matrix(cov.frm, design$variables))
  
  
  
  # domain labels as indexes for use in INLA
  domain.table$domain.struct <- seq_len(nrow(domain.table))
  mod.dat$domain.struct <- match(mod.dat$domain, domain.table$domain)
  pop.dat$domain.struct <- match(pop.dat$domain, domain.table$domain)
  
  # model formula
  ftxt <- paste("resp ~ 1")
  if (length(all.vars(cov.frm)) > 0) {
    ftxt <- 
      paste(ftxt, as.character(cov.frm)[-1], sep = " + ")
  }
  
  if (is.null(adj.mat)) {
    # set priors
    hyperpc.iid.int <- list(
      prec = list(prior = "pc.prec",
                  param = c(pc.u , pc.alpha))
    )
    warning("No spatial information provided, using iid domain effects")
    model.method <- "iid.model"
    ftxt <- paste0(ftxt, " + f(domain.struct, model = 'iid', hyper = hyperpc.iid.int)")
  } else {
    hyperpc.bym.int <- list(
      prec = list(prior = "pc.prec", param = c(pc.u , pc.alpha)),  
      phi = list(prior = 'pc', param = c(pc.u.phi , pc.alpha.phi))
    )
    model.method <- "bym2.model"
    ftxt <- paste0(ftxt, " + f(domain.struct, model = 'bym2', graph=adj.mat, hyper = hyperpc.bym.int)")
  }
  # fit model
  mod.frm <- as.formula(ftxt)
  fit <- INLA::inla(mod.frm, family = family, data = mod.dat,
                    control.compute = list(dic = TRUE, mlik = TRUE,
                                           cpo = TRUE, config = TRUE),
                    control.predictor = list(compute = TRUE),
                    lincomb = NULL, quantiles = c((1-level)/2, 0.5, 1-(1-level)/2))
  
  
  # generate draws
  samp.all <- INLA::inla.posterior.sample(n = n.sample, result = fit, intern = TRUE)

  # identify indices of fixed effects and random effects
  fe.idx <- grep(colnames(mm.pop)[1], rownames(samp.all[[1]]$latent))
  fe.idx <- fe.idx:(fe.idx + ncol(mm.pop) - 1)
  re.idx <- grep("domain.struct", x = rownames(samp.all[[1]]$latent))
  
  # aggregate sample predictions
  summary.sample <- function(x) {
    pop.unit.ests <-
      x$latent[re.idx][pop.dat$domain.struct] + mm.pop %*% x$latent[fe.idx] 
    if (family == "binomial") {
      pop.unit.ests <- expit(pop.unit.ests)
    }
    if (!is.null(X.pop.weights)) {
      area.ests <- 
        aggregate(pop.unit.ests * X.pop.weights,  list(domain = pop.dat$domain.struct), sum)
    } else {
      area.ests <- 
        aggregate(pop.unit.ests, list(domain = pop.dat$domain.struct), mean)
    }
    
    return(area.ests[match(1:nrow(domain.table), area.ests[, 1]), 2])
  }
  est.mat <- do.call(cbind, lapply(samp.all, summary.sample))
  out[[paste0(model.method, ".fit")]] <- fit
  out[[paste0(model.method, ".est")]]  <-
    data.frame(domain = domain.table$domain,
               mean = rowMeans(est.mat),
               median = apply(est.mat, 1,
                              function(x) median(x, na.rm = TRUE)),
               var = apply(est.mat, 1, var),
               lower = apply(est.mat, 1,
                             function(x) quantile(x, (1-level)/2, na.rm = TRUE)),
               upper = apply(est.mat, 1,
                             function(x) quantile(x, 1-(1-level)/2, na.rm = TRUE)),
               method = paste0("Unit level model: ", 
                               ifelse(is.null(adj.mat), "IID", "BYM2")))
  if (return.samples) {
    out[[paste0(model.method, ".sample")]]  <- est.mat
  }
  out[[paste0(model.method, ".est")]] <- 
    out[[paste0(model.method, ".est")]][match(out$direct.est$domain, 
                                              out[[paste0(model.method, ".est")]]$domain),]
  attr(out, "method.names") <- c(attr(out, "method.names"), paste0(model.method, ".est"))
  attr(out, "inla.fitted") <- c(attr(out, "inla.fitted"), model.method)
  
  # finish building return value
  out$call <- match.call()
  class(out) <- c("svysae")
  return(out)
}



