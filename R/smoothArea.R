#' Small area estimation via basic area level model
#' 
#' Generates small area estimates  by smoothing direct estimates using an
#' area level model
#'
#' @param formula An object of class 'formula' describing the model to be fitted. 
#'   If direct.est is specified, the right hand side of the formula is not necessary.
#' @param domain One-sided formula specifying factors containing domain labels
#' @param design An object of class "svydesign" containing the data for the model 
#' @param adj.mat Adjacency matrix with rownames matching the domain labels. If set to NULL, the IID spatial effect will be used. 
#' @param X.domain Data frame of areal covariates. One of the column names needs to match the name of the domain variable, in order to be linked to the data input. Currently only supporting time-invariant covariates.
#' @param direct.est Data frame of direct estimates, with first column containing the domain variable, second column containing direct estimate, and third column containing the variance of direct estimate. 
#' @param domain.size Data frame of domain sizes. One of the column names needs to match the name of the domain variable, in order to be linked to the data input and there must be a column names 'size' containing domain sizes.
#' @param transform Optional transformation applied to the direct estimates before fitting area level model. The default option is no transformation, but logit and log are implemented.
#' @param pc.u 	Hyperparameter U for the PC prior on precisions. See the INLA documentation for more details on the parameterization.
#' @param pc.alpha Hyperparameter alpha for the PC prior on precisions.
#' @param pc.u.phi Hyperparameter U for the PC prior on the mixture probability phi in BYM2 model.
#' @param pc.alpha.phi Hyperparameter alpha for the PC prior on the mixture probability phi in BYM2 model.
#' @param level The specified level for the posterior credible intervals
#' @param n.sample Number of draws from posterior used to compute summaries
#' @param var.tol Tolerance parameter; if variance of an area's direct estimator is below this value, that direct estimator is dropped from model
#' @param return.samples If TRUE, return matrix of posterior samples of area level quantities
#'
#' @importFrom stats model.frame model.matrix model.response qnorm
#' @importFrom utils combn
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
#'                   weights = ~weights, data = DemoData2, nest = TRUE)
#' Xmat <- aggregate(age~region, data = DemoData2, FUN = mean)
#' 
#' # EXAMPLE 1: Continuous response model
#' cts.res <- smoothArea(tobacco.use ~ 1,
#'                       domain = ~region,
#'                       design = des0,
#'                       adj.mat = DemoMap2$Amat, 
#'                       pc.u = 1,
#'                       pc.alpha = 0.01,
#'                       pc.u.phi = 0.5,
#'                       pc.alpha.phi = 2/3)
#' 
#' # EXAMPLE 2: Including area level covariates
#' cts.cov.res <- smoothArea(tobacco.use ~ age, 
#'                           domain = ~region,
#'                           design = des0,
#'                           adj.mat = DemoMap2$Amat, 
#'                           X.domain = Xmat,
#'                           pc.u = 1,
#'                           pc.alpha = 0.01,
#'                           pc.u.phi = 0.5,
#'                           pc.alpha.phi = 2/3)
#' 
#' # EXAMPLE 3: Binary response model
#' bin.res <- smoothArea(tobacco.use ~ 1, 
#'                       domain = ~region,
#'                       design = des0,
#'                       adj.mat = DemoMap2$Amat, 
#'                       transform = "logit",
#'                       pc.u = 1,
#'                       pc.alpha = 0.01,
#'                       pc.u.phi = 0.5,
#'                       pc.alpha.phi = 2/3)
#' 
#' # EXAMPLE 4: Including area level covariates in binary response model
#' bin.cov.res <- smoothArea(tobacco.use ~ age, 
#'                           domain = ~region,
#'                           design = des0,
#'                           adj.mat = DemoMap2$Amat, 
#'                           transform = "logit",
#'                           X.domain = Xmat,
#'                           pc.u = 1,
#'                           pc.alpha = 0.01,
#'                           pc.u.phi = 0.5,
#'                           pc.alpha.phi = 2/3)
#'}
smoothArea <- function(formula,
                       domain,
                       design = NULL,
                       adj.mat = NULL,
                       X.domain = NULL,
                       direct.est = NULL,
                       domain.size = NULL,
                       transform = c("identity", "logit", "log"),
                       pc.u = 1,
                       pc.alpha = 0.01, 
                       pc.u.phi = 0.5,
                       pc.alpha.phi = 2/3,
                       level = .95, 
                       n.sample = 250,
                       var.tol = 1e-10,
                       return.samples = F) {
  
  # setup output
  out <- list()
  attr(out, "inla.fitted") <- c()
  
  # identify domain variable
  domain.var <- all.vars(domain)
  if (length(domain.var) != 1) {
    stop("Domain labels must be contained in a single column. Spatio-temporal not supported yet.")
  } 
  
  # response formula
  resp.frm <- as.formula(paste0("~", all.vars(update(formula, . ~ NULL))[1]))
  # covariate formula
  cov.frm <- update(formula, NULL ~ .)
  
  # calculate direct estimates via svyby
  if (is.null(direct.est)) {
    # known domain sizes
    if (!is.null(domain.size)) {
      direct.est <- 
        survey::svyby(resp.frm, domain, design = design, svytotal, na.rm = TRUE)
      direct.est$domain.size <-
        domain.size[match(direct.est[, 1], domain.size[, 1]), 2]
      direct.est[, 2] = direct.est[, 2] / direct.est$domain.size
      direct.est[, 3] = (direct.est[, 3] / direct.est$domain.size) ^ 2
      direct.est <- direct.est[, 1:3]
    } else {
      direct.est <-
        survey::svyby(resp.frm, domain, design = design, svymean, na.rm = TRUE)
      direct.est[, 3] <- direct.est[, 3] ^ 2
    }
  } else {
    # assumes structure of direct.est
    direct.est <- direct.est
  }
  colnames(direct.est) <- c("domain", "mean", "var")
  direct.est <- 
    data.frame(domain = direct.est$domain,
               mean = direct.est$mean,
               median = direct.est$mean,
               var = direct.est$var,
               lower = direct.est$mean + 
                 qnorm((1-level)/2) * sqrt(direct.est$var),
               upper = direct.est$mean +
                 qnorm(1 - (1-level)/2) * sqrt(direct.est$var),
               method = paste0("Direct"))

  # prepare data for modeling, removing any areas with zero/low sampling variances
  mod.dat <- direct.est
  mod.dat$mean <- 
    ifelse(mod.dat$var > var.tol, mod.dat$mean, NA)
  mod.dat$var <-
    ifelse(mod.dat$var > var.tol, mod.dat$var, NA)
  mod.dat$prec  <- 1 / mod.dat$var
  
  # identify domains for estimation (including estimates for regions in X.domain or adj.mat)
  if (!is.null(X.domain)) {
    if (!is.null(adj.mat) & !setequal(X.domain[[domain.var]], rownames(adj.mat))) {
      stop("Domains in X.domain do not match domains in adj.mat.")
    }
    if (any(is.na(match(X.domain[[domain.var]], mod.dat$domain)))) {
      warning(cat("There are domains in X.domain not in design/direct estimates.",
                  "\nGenerating estimates for all domains in X.domain."))
    }
    mod.X.domain <- X.domain
    mod.X.domain$domain <- X.domain[[domain.var]]
    mod.X.domain <- mod.X.domain[, names(mod.X.domain) != domain.var]
    mod.dat <- merge(mod.dat, mod.X.domain,  by = "domain", all.y = TRUE)
    direct.est <- 
      merge(direct.est, data.frame(domain = mod.X.domain$domain), 
            by = "domain", all.y = TRUE)
    direct.est$method = "Direct"
  } else if(!is.null(adj.mat)) {
    if (any(is.na(match(mod.dat$domain, rownames(adj.mat))))) {
      stop("Domains in adj.mat do not match domains in design/direct estimates.")
    }
    if (any(is.na(match(rownames(adj.mat), mod.dat$domain)))) {
      warning(cat("There are domains in adj.mat not in design/direct estimates.",
                  "\nGenerating estimates for all domains in adj.mat."))
    }
    mod.dat <- 
      merge(mod.dat, data.frame(domain = rownames(adj.mat)), 
            by = "domain", all.y = TRUE)
    direct.est <- 
      merge(direct.est, data.frame(domain = rownames(adj.mat)), 
            by = "domain", all.y = TRUE)
    direct.est$method = "Direct"
  }
  mod.dat$domain.id <- 1:nrow(mod.dat)
  mm.domain <- model.matrix(cov.frm, mod.dat)
  out$direct.est <- direct.est
  attr(out, "domain.names") <- sort(direct.est$domain)
  attr(out, "method.names") <- c("direct.est")

  transform <- match.arg(transform)
  attr(out, "transform") <- transform
  
  # apply transformation to direct estimates
  if (transform == "identity") {
    h <- function(x) x
    h.inv <- function(x) x
    if (all(mod.dat$mean < 1 & mod.dat$mean > 0)) {
      warning("Direct estimates appear to be proportions. You may want to consider using transform = 'logit'.")
    }
  } else if (transform == "logit") {
    h <- function(x) logit(x)
    h.inv <- function(x) expit(x)
    mod.dat$prec  <- 
      (mod.dat$mean^2*(1-mod.dat$mean)^2) / mod.dat$var
    mod.dat$mean <- h(mod.dat$mean)
  } else if (transform == "log") {
    h <- function(x) log(x)
    h.inv <- function(x) exp(x)
    mod.dat$prec  <- 
      mod.dat$mean^2 / mod.dat$var
    mod.dat$mean <- h(mod.dat$mean)
  }
  mod.dat$scale <- mod.dat$prec 
  
  # prepare formula for INLA
  s.dir.ftxt <- 
    paste("mean ~ 1")
  if (length(all.vars(cov.frm)) > 0) {
    s.dir.ftxt <- 
      paste(s.dir.ftxt, paste(all.vars(cov.frm), collapse = " + "), sep = " + ")
  }
  iid.model.ftxt <- 
    paste0(s.dir.ftxt, " + f(domain.id, model = 'iid', hyper = hyperpc.iid.int)")
  # set priors
  hyperpc.iid.int <- list(
    prec = list(prior = "pc.prec",
                param = c(pc.u , pc.alpha))
  )
  # SMOOTHED DIRECT w/ IID RE
  iid.model.fit <-
    INLA::inla(as.formula(iid.model.ftxt),
               family = "gaussian", data = mod.dat, 
               scale = mod.dat$prec ,
               control.family = 
                 list(hyper = list(prec = list(initial= log(1), fixed= TRUE))),
               control.predictor = list(compute = TRUE),
               control.compute=list(config = TRUE))
  sample.id <-c(list(domain.id = 1:nrow(mod.dat)),
                as.list(setNames(rep(1, ncol(mm.domain)), colnames(mm.domain))))
  iid.model.sample <-
    INLA::inla.posterior.sample(n = n.sample, iid.model.fit, sample.id)
  fe.idx <- grep(colnames(mm.domain)[1], rownames(iid.model.sample[[1]]$latent))
  fe.idx <- fe.idx:(fe.idx + ncol(mm.domain) - 1)
  
  # sample domain means from posterior
  iid.model.mat <-
    do.call(cbind, lapply(iid.model.sample,
                          function(x) x$latent[1:nrow(mod.dat)] +
                            mm.domain %*% x$latent[fe.idx]))
  
  iid.model.mat <- h.inv(iid.model.mat)
  out$iid.model.fit <- iid.model.fit
  
  out$iid.model.est <- 
    data.frame(domain = mod.dat$domain,
               mean = rowMeans(iid.model.mat),
               median = apply(iid.model.mat, 1,
                              function(x) median(x, na.rm = TRUE)),
               var = apply(iid.model.mat, 1, var),
               lower = apply(iid.model.mat, 1,
                             function(x) quantile(x, (1-level)/2, na.rm = TRUE)),
               upper = apply(iid.model.mat, 1,
                             function(x) quantile(x, 1-(1-level)/2, na.rm = TRUE)),
               method = paste0("Area level model: IID"))

  out$iid.model.est <- 
    out$iid.model.est[match(out$direct.est$domain, out$iid.model.est$domain),]
  rownames(out$iid.model.est) <- NULL
  if (return.samples) {
    out$iid.model.sample <- iid.model.mat
  } else {
    out$iid.model.sample <- NULL
  }
  attr(out, "method.names") <- c(attr(out, "method.names"), "iid.model.est")
  attr(out, "inla.fitted") <- c(attr(out, "inla.fitted"), "iid.model")
  
  # SMOOTHED DIRECT w/ BYM2 RE
  if (!is.null(adj.mat)) {

    mod.dat$domain.id <- match(mod.dat$domain, rownames(adj.mat))
    mod.dat <- mod.dat[match(1:nrow(mod.dat), mod.dat$domain.id), ]
    mm.domain <- model.matrix(cov.frm, mod.dat)
    hyperpc.bym.int <- list(
      prec = list(prior = "pc.prec", param = c(pc.u , pc.alpha)),  
      phi = list(prior = 'pc', param = c(pc.u.phi , pc.alpha.phi))
    )
    # prepare formula for INLA
    bym2.model.ftxt <- 
      paste0(s.dir.ftxt,
             "+ f(domain.id, model = 'bym2', graph = adj.mat,", 
             "hyper = hyperpc.bym.int, scale.model = TRUE)")
    bym2.model.fit <-
      INLA::inla(as.formula(bym2.model.ftxt),
                 family = "gaussian", data = mod.dat, 
                 scale = mod.dat$prec ,
                 control.family = 
                   list(hyper = list(prec = list(initial= log(1), fixed= TRUE))),
                 control.predictor = list(compute = TRUE),
                 control.compute=list(config = TRUE))
    sample.id <-c(list(domain.id = 1:nrow(adj.mat)),
                  as.list(setNames(rep(1, ncol(mm.domain)), colnames(mm.domain))))
    bym2.model.sample <-
      INLA::inla.posterior.sample(n = n.sample, bym2.model.fit, sample.id)
    fe.idx <- grep(colnames(mm.domain)[1], rownames(bym2.model.sample[[1]]$latent))
    fe.idx <- fe.idx:(fe.idx + ncol(mm.domain) - 1)
    
    # sample domain means from posterior
    bym2.model.mat <-
      do.call(cbind, lapply(bym2.model.sample,
                            function(x) x$latent[1:nrow(adj.mat)] +
                              mm.domain %*% x$latent[fe.idx]))
    bym2.model.mat <- h.inv(bym2.model.mat)
    out$bym2.model.fit <- bym2.model.fit
    out$bym2.model.est <- 
      data.frame(domain = mod.dat$domain,
                 mean = rowMeans(bym2.model.mat),
                 median = apply(bym2.model.mat, 1,
                                function(x) median(x, na.rm = TRUE)),
                 var = apply(bym2.model.mat, 1, var),
                 lower = apply(bym2.model.mat, 1,
                               function(x) quantile(x, (1-level)/2, na.rm = TRUE)),
                 upper = apply(bym2.model.mat, 1,
                               function(x) quantile(x, 1-(1-level)/2, na.rm = TRUE)),
                 method = paste0("Area level model: BYM2"))
    out$bym2.model.est <- 
      out$bym2.model.est[match(out$direct.est$domain, out$bym2.model.est$domain),]
    rownames(out$bym2.model.est) <- NULL
    if (return.samples) {
      out$bym2.model.sample <- bym2.model.mat
    } else {
      out$bym2.model.sample <- NULL
    }
    attr(out, "method.names") <- c(attr(out, "method.names"), "bym2.model.est")
    attr(out, "inla.fitted") <- c(attr(out, "inla.fitted"), "bym2.model")
  }
  out$call <- match.call()
  class(out) <- c("svysae")
  return(out)
}

#' @method print svysae
#' @export
print.svysae <- function(x, ...) {
  x_att <- attributes(x)
  
  # print out call
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  # print estimation methods
  cat("Methods used: ")
  cat(x_att$method.names, sep = ", ")
  cat("\n\n")
  
  # print transform, if used
  if (!is.null(x_att$transform)) {
    cat("Transform: ")
    cat(x_att$transform)
    cat("\n\n")
  }
  
  # print estimates
  for (i in x_att$method.names) {
    cat(i, "\n")
    print(x[[i]])
    cat("\n")
  }
}

#' @method plot svysae
#' @export
plot.svysae <- function(x, return_list = F, plot.factor = NULL, ...) {
  combined_est <- do.call(rbind, x[attr(x, "method.names")])
  
  # pull the first method (should be "direct.est")
  ref_method <- attr(x, "method.names")[1]
  m <- nrow(x[[ref_method]])
  sorted_levels <- x[[ref_method]]$domain[order(x[[ref_method]]$mean)]
  combined_est$domain <- factor(combined_est$domain, 
                                levels = sorted_levels)
  
  # split across multiple plots for many estimates
  if (m > 30 & is.null(plot.factor)) {
    plot_breaks <- c(seq(0, m, by = 30), m)
    plot_labels <- paste0(
      "Areas ", 
      (seq_along(plot_breaks)[-1] - 2) * 30 + 1,
      " to ",
      pmin((seq_along(plot_breaks)[-1] - 1) * 30, m)
    )
    combined_est$plot <- 
      cut(1:m, breaks = plot_breaks, labels = plot_labels)
  } else if (is.null(plot.factor)) {
    combined_est$plot <- factor(1, levels = 1)
  }
  plot_list <- lapply(
    split(combined_est, combined_est$plot),
    function(x) {
      ggplot2::ggplot(x, ggplot2::aes(x = domain, y = mean, color = method)) +
        ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.5)) + 
        ggplot2::geom_linerange(ggplot2::aes(x = domain, ymin = lower, ymax = upper), 
                       position = ggplot2::position_dodge(width = 0.5)) + 
        ggplot2::scale_color_discrete(name = "Method") + 
        ggplot2::ylab("Small Area Estimate") + ggplot2::xlab("Domain") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
    }
  )
  if (!return_list) {
    for (i in seq_along(plot_list)) {
      print(plot_list[[i]])
    }
  } else {
    return(plot_list)
  }
}

#' Plot heatmap comparing pairwise posterior exceedence probabilities for svysae object
#' 
#'
#' @param x svysae object. Plots are created for all models in this object.
#' @param posterior.sample Matrix of posteriors samples of area level quantities with one row for each area and one column for each sample. This argument may be specified to only provide a heatmap for the desired samples.
#'
#' @return ggplot containing heat map of pairwise comparisons
#' 
#' @export 
#'
#' @examples 
#' \dontrun{
#' data(DemoData2)
#' data(DemoMap2)
#' library(survey)
#' des0 <- svydesign(ids = ~clustid+id, strata = ~strata,
#'                   weights = ~weights, data = DemoData2, nest = TRUE)
#' Xmat <- aggregate(age~region, data = DemoData2, FUN = mean)
#' 
#' cts.res <- smoothArea(tobacco.use ~ 1,
#'                       domain = ~region,
#'                       design = des0,
#'                       adj.mat = DemoMap2$Amat, 
#'                       pc.u = 1,
#'                       pc.alpha = 0.01,
#'                       pc.u.phi = 0.5,
#'                       pc.alpha.phi = 2/3,
#'                       return.samples = TRUE)
#' compareEstimates(cts.res)
#' }
compareEstimates <- function(x,
                             posterior.sample = NULL) {
  
  
  x_att <- attributes(x)
  if (is.null(posterior.sample)) {
    sample_list <- x[paste0(x_att$inla.fitted, ".sample")]
  } else {
    sample_list <- list(posterior.sample)
  }
  for (i in seq_along(sample_list)) {
    current_samples <- sample_list[[i]]
    domain_medians <- apply(current_samples, 1, median)
    median_order <- order(domain_medians)
    
    current_samples <- current_samples[median_order, ]
    # get indices for all possible pairs of admin 1 areas
    domain_comb <- combn(length(x_att$domain.names), 2)
    
    plot_dat <- data.frame(
      Domain1 = x_att$domain.names[median_order][domain_comb[1,]],
      Domain2 = x_att$domain.names[median_order][domain_comb[2,]]
    )
    plot_dat$Domain1 <- factor(plot_dat$Domain1,
                               x_att$domain.names[median_order])
    plot_dat$Domain2 <- factor(plot_dat$Domain2,
                               x_att$domain.names[median_order])
    plot_dat$Prob = apply(
      domain_comb, 2, 
      function(x) mean(current_samples[x[1],] > current_samples[x[2],])
    )
    plot_dat$est = NA
    
    median_dat <- data.frame(
      Domain1 = x_att$domain.names[median_order], 
      Domain2 = "Median", 
      Prob = NA,
      est = domain_medians
      
    )
    extra_cols <- c("", "Median", "Interval")
    # combine into single tibble for plotting
    plot_dat <- rbind(plot_dat, median_dat)
    
    plot_dat$Domain2 <- 
      factor(plot_dat$Domain2, levels = c(levels(plot_dat$Domain1), extra_cols))
    
    g_heat <- ggplot2::ggplot(data = plot_dat,
                              ggplot2::aes(x = Domain2, y = Domain1, fill = Prob)) + 
      ggplot2::theme_minimal() + 
      ggplot2::theme(legend.position="bottom", 
                     axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1, size = 10), 
                     panel.grid.major = ggplot2::element_blank(), 
                     panel.grid.minor = ggplot2::element_blank(),
                     text = ggplot2::element_text(size = 10)) +
      ggplot2::geom_tile() + 
      ggplot2::geom_text(
        ggplot2::aes(
          label = ifelse(is.na(Prob), "", sprintf("%0.2f", round(Prob, digits = 2))),
            color = Prob > .5
          ), 
        size = 5
      ) + 
      ggplot2::geom_text(
        data = plot_dat,
        ggplot2::aes(
          label = ifelse(is.na(est), "", sprintf("%0.2f", round(est, digits = 2)))
          ),
        size = 5
      ) +
      ggplot2::coord_equal() + 
      ggplot2::scale_fill_viridis_c(name = "Probability of Domain 1 > Domain 2",
                                    na.value = "white") + 
      ggplot2::scale_color_manual(values = c("white", "black"),  guide = "none") +
      ggplot2::labs(title = x_att$inla.fitted[i],
                    x = "Domain 2",
                    y = "Domain 1")
    suppressWarnings(print(g_heat))
  }
}
#' Mapping estimates for svysae object
#'
#' 
#' @param x syvsae object
#' @param geo.data sf object containing polygon data for the small areas. One of the columns should be named domain and contain the domain labels.
#' @param variable The posterior summary variable to plot. May be one of "median", "mean", or "var".
#' @param viridis.option viridis color scheme
#'
#' @return ggplot containing map of small area posterior summary statistics
#' 
#' @export
#'
#' @examples 
#' \dontrun{
#' data(DemoData2)
#' data(DemoMap2)
#' library(survey)
#' des0 <- svydesign(ids = ~clustid+id, strata = ~strata,
#'                   weights = ~weights, data = DemoData2, nest = TRUE)
#' Xmat <- aggregate(age~region, data = DemoData2, FUN = mean)
#' geo.data <- sf::st_as_sf(DemoMap2$geo)
#' geo.data$domain <- geo.data$REGNAME
#' cts.res <- smoothArea(tobacco.use ~ 1,
#'                       domain = ~region,
#'                       design = des0,
#'                       adj.mat = DemoMap2$Amat, 
#'                       pc.u = 1,
#'                       pc.alpha = 0.01,
#'                       pc.u.phi = 0.5,
#'                       pc.alpha.phi = 2/3,
#'                       return.samples = TRUE)
#' mapEstimates(cts.res, geo.data = geo.data, variable = "median")
#' mapEstimates(cts.res, geo.data = geo.data, variable = "var")
#' }
mapEstimates <- function(x, geo.data, variable, viridis.option = "viridis") {
  combined_est <- do.call(rbind, x[attr(x, "method.names")])
  
  # join estimates and geo, by domain column
  plot_dat <- merge(geo.data, combined_est, all.x = TRUE)
  ggplot2::ggplot(plot_dat, ggplot2::aes(fill = .data[[variable]])) + 
    ggplot2::geom_sf() + 
    ggplot2::facet_wrap(~method) + 
    ggplot2::theme_bw() +
    ggplot2::scale_fill_viridis_c(name = variable, option = viridis.option) 
}


