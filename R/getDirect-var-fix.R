# Rule 2 from AddPhantom_U5MR.R (Qianyu Dong), applied to normalized PSU IDs.
.u5mr_rule2 <- function(modt, area.col, tol = 1e-12) {
    area_vec <- sort(unique(modt[[area.col]]))
    vapply(area_vec, function(a) {
        d <- modt[modt[[area.col]] == a, ]
        if (length(unique(d$cluster)) == 1) return(TRUE)
        ok <- TRUE
        for (b in unique(d$age)) {
            db <- d[d$age == b, ]
            E <- sum(db$weight * db$exposure)
            if (E <= 0) next
            p.area <- sum(db$weight * db$died) / E
            Ec <- tapply(db$weight * db$exposure, db$cluster, sum)
            Dc <- tapply(db$weight * db$died, db$cluster, sum)
            p.cl <- Dc[Ec > 0] / Ec[Ec > 0]
            if (max(abs(p.cl - p.area)) > tol) { ok <- FALSE; break }
        }
        ok
    }, logical(1))
}

.u5mr_fix_audit <- function() {
    data.frame(region = character(), years = character(), status = character(),
        reason = character(), trigger = character(), real.psus = integer(),
        phantom.psus = integer(), phantom.rows = integer(),
        reference.groups = character(), phantom.weighted.deaths = numeric(),
        phantom.weighted.exposure = numeric(), stringsAsFactors = FALSE)
}

.prepare_u5mr_fix <- function(births, cluster.formula, Ntrials, refStrataVar) {
    vars <- all.vars(cluster.formula)
    terms <- stats::terms(cluster.formula)
    named.terms <- vapply(attr(terms, "term.labels"), function(term) {
        expr <- parse(text = term)[[1L]]
        is.symbol(expr) && as.character(expr) %in% vars
    }, logical(1))
    if (!length(vars) || attr(terms, "response") != 0L || !all(named.terms))
        stop("var.fix requires a simple named cluster formula, such as ~psu + household.")
    if (!all(vars %in% names(births))) stop("Cluster variables are missing from births.")
    if (!is.null(refStrataVar) && (!is.character(refStrataVar) || length(refStrataVar) != 1L ||
        is.na(refStrataVar) || !refStrataVar %in% names(births)))
        stop("refStrataVar must name one reference-group column in births.")
    exposure <- if (is.null(Ntrials)) rep(1, nrow(births)) else births[[Ntrials]]
    if (!is.numeric(exposure) || length(exposure) != nrow(births) ||
        any(!is.finite(exposure) | exposure < 0))
        stop("var.fix requires finite, nonnegative exposure counts in Ntrials.")
    if (!is.numeric(births$died) || any(!is.finite(births$died) |
        births$died < 0 | births$died > exposure))
        stop("var.fix requires finite death counts between zero and exposure.")
    if (!is.numeric(births$weights0) || any(!is.finite(births$weights0) | births$weights0 < 0))
        stop("var.fix requires finite, nonnegative sampling weights.")
    required <- c(vars, "strata", "time0", "age0", refStrataVar)
    if (anyNA(births[, required, drop = FALSE]))
        stop("var.fix requires nonmissing cluster, stratum, time, age and reference values.")
    dat <- data.frame(region0 = as.character(births$region0),
        time0 = as.character(births$time0), age0 = as.character(births$age0),
        weights0 = births$weights0, strata = as.character(births$strata),
        reference = if (is.null(refStrataVar)) rep("All", nrow(births)) else as.character(births[[refStrataVar]]),
        n = exposure, died = births$died, stringsAsFactors = FALSE)
    stage.vars <- paste0("stage", seq_along(vars))
    for (i in seq_along(vars)) dat[[stage.vars[i]]] <- as.character(births[[vars[i]]])
    mapping <- unique(dat[, c("time0", "strata", "reference")])
    if (anyDuplicated(mapping[, c("time0", "strata")]))
        stop("Each design stratum must map to exactly one reference group within each period.")
    retain <- exposure > 0
    keys <- c(stage.vars, "region0", "time0", "age0", "weights0", "strata", "reference")
    dat <- as.data.frame(data.table::as.data.table(dat[retain, ])[,
        list(n = sum(get("n")), died = sum(get("died"))), by = keys])
    dat$psu <- paste(match(dat$strata, unique(dat$strata)),
                     match(dat$stage1, unique(dat$stage1)), sep = ":")
    # Cache compact-data row numbers by period and area; missing area labels
    # remain in each period's national rows.
    index <- lapply(split(seq_len(nrow(dat)), dat$time0), function(rows) {
        list(all = rows, areas = split(rows, dat$region0[rows]))
    })
    # Choose an unused ID prefix once for all target-specific phantom PSUs.
    prefix <- ".SUMMER.phantom."
    real.ids <- unlist(dat[stage.vars], use.names = FALSE)
    while (any(startsWith(real.ids, prefix))) prefix <- paste0(prefix, ".")
    list(data = dat, reference = .u5mr_reference(dat), stages = stage.vars,
        formula = stats::reformulate(stage.vars), retain = retain,
        regions = sort(unique(as.character(births$region0[!is.na(births$region0)]))),
        index = index, phantom.prefix = prefix)
}

.u5mr_reference <- function(dat) {
    if (!nrow(dat)) return(data.frame())
    dt <- data.table::as.data.table(dat)
    psus <- dt[, list(exposure = sum(get("n")),
        weighted.exposure = sum(get("weights0") * get("n"))),
        by = c("time0", "reference", "psu")]
    size <- psus[, list(K = .N,
        weight.ph = mean(get("weighted.exposure") / get("exposure"))),
        by = c("time0", "reference")]
    bands <- dt[, list(weighted.exposure = sum(get("weights0") * get("n")),
        weighted.deaths = sum(get("weights0") * get("died"))),
        by = c("time0", "reference", "age0")]
    ref <- as.data.frame(merge(bands, size, by = c("time0", "reference"), sort = FALSE))
    ref$exposure.ph <- ref$weighted.exposure / (ref$K * ref$weight.ph)
    # This is p_a E*_a, also defined as zero when a band has zero weighted exposure.
    ref$deaths.ph <- ref$weighted.deaths / (ref$K * ref$weight.ph)
    ref
}

.u5mr_target <- function(fix, area, time) {
    period <- match(as.character(time), names(fix$index))
    if (is.na(period)) return(fix$data[integer(), , drop = FALSE])
    index <- fix$index[[period]]
    if (area == "All") {
        rows <- index$all
    } else {
        region <- match(as.character(area), names(index$areas))
        rows <- if (is.na(region)) integer() else index$areas[[region]]
    }
    fix$data[rows, , drop = FALSE]
}

.u5mr_target_rule2 <- function(dat) {
    if (!nrow(dat)) return(FALSE)
    frame <- data.frame(area = "target", cluster = dat$psu, age = dat$age0,
        weight = dat$weights0, exposure = dat$n, died = dat$died)
    unname(.u5mr_rule2(frame, "area"))
}

.u5mr_augment <- function(fix, target, area, time) {
    strata <- unique(target[, c("strata", "reference")])
    prefix <- fix$phantom.prefix
    phantoms <- vector("list", nrow(strata))
    for (i in seq_len(nrow(strata))) {
        ref <- fix$reference[fix$reference$time0 == as.character(time) &
            fix$reference$reference == strata$reference[i], , drop = FALSE]
        ref <- ref[is.finite(ref$exposure.ph) & ref$exposure.ph > 0, , drop = FALSE]
        ph <- fix$data[rep(1L, nrow(ref)), , drop = FALSE]
        ph$region0 <- rep(as.character(area), nrow(ref))
        ph$time0 <- rep(as.character(time), nrow(ref))
        ph$age0 <- ref$age0
        ph$strata <- rep(strata$strata[i], nrow(ref))
        ph$reference <- rep(strata$reference[i], nrow(ref))
        ph$weights0 <- ref$weight.ph
        ph$n <- ref$exposure.ph
        ph$died <- ref$deaths.ph
        for (j in seq_along(fix$stages))
            ph[[fix$stages[j]]] <- rep(paste0(prefix, i, ".", j), nrow(ref))
        ph$psu <- rep(paste0(prefix, i), nrow(ref))
        phantoms[[i]] <- ph
    }
    ph <- do.call(rbind, phantoms)
    augmented <- rbind(fix$data, ph)
    augmented$died <- augmented$died / augmented$n
    audit <- data.frame(region = as.character(area), years = as.character(time),
        status = "unresolved", reason = "pending", trigger = if (length(unique(target$psu)) == 1L)
            "single_psu" else "equal_hazards", real.psus = length(unique(target$psu)),
        phantom.psus = length(unique(ph$psu)), phantom.rows = nrow(ph),
        reference.groups = paste(unique(strata$reference), collapse = ";"),
        phantom.weighted.deaths = sum(ph$weights0 * ph$died),
        phantom.weighted.exposure = sum(ph$weights0 * ph$n), stringsAsFactors = FALSE)
    list(data = augmented, audit = audit,
        degenerate = .u5mr_target_rule2(rbind(target, ph)))
}
