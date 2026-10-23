fix_fixture <- function() {
    x <- expand.grid(psu = 1:4, age = c("0", "1-1"), time = c("p1", "p2"),
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    x$region <- c("A", "B", "B", "C")[x$psu]
    x$strata <- c("s1", "s1", "s2", "s2")[x$psu]
    x$reference <- c("urban", "urban", "rural", "rural")[x$psu]
    x$house <- 1L
    x$weights <- c(1, 2, 3, 4)[x$psu]
    x$total <- 10
    x$died <- c(0, 2, 6, 1)[x$psu] + as.integer(x$age == "1-1" & x$psu == 2)
    x$died[x$time == "p2"] <- x$died[x$time == "p2"] + 1
    x
}

fix_direct <- function(x, ..., compact = TRUE) {
    suppressMessages(getDirect(x, years = c("p1", "p2", "empty"),
        clusterVar = ~psu + house, weightsVar = "weights",
        Ntrials = if (compact) "total" else NULL, ...))
}

fix_prepare <- function(x, ref = NULL) {
    x$region0 <- x$region; x$time0 <- x$time
    x$age0 <- x$age; x$weights0 <- x$weights
    .prepare_u5mr_fix(x, ~psu + house, "total", ref)
}

fix_expand <- function(x) {
    do.call(rbind, lapply(seq_len(nrow(x)), function(i) {
        z <- x[rep(i, x$total[i]), ]
        z$died <- c(rep(1, x$died[i]), rep(0, x$total[i] - x$died[i]))
        z
    }))
}

fix_table <- function(x) {
    attr(x, "fixed_areas") <- attr(x, "variance.fix") <- NULL
    x
}

test_that("indexed targets preserve row order, missing areas and empty domains", {
    x <- fix_fixture()
    x$region <- c("", "a.b", "a.b", NA_character_)[x$psu]
    x$time <- ifelse(x$time == "p1", "", "p.2")
    x <- x[rev(seq_len(nrow(x))), ]
    p <- fix_prepare(x)
    for (time in c("", "p.2", "absent")) {
        for (area in c("All", "", "a.b", "absent")) {
            expected <- p$data[p$data$time0 == time, , drop = FALSE]
            if (area != "All")
                expected <- expected[!is.na(expected$region0) & expected$region0 == area, , drop = FALSE]
            expect_identical(.u5mr_target(p, area, time), expected)
            rows <- .u5mr_target(p, area, time, rows.only = TRUE)
            expect_identical(p$data[rows, , drop = FALSE], expected)
        }
    }
    x$total <- x$died <- 0
    p <- fix_prepare(x)
    expect_identical(.u5mr_target(p, "All", ""), p$data)
})

test_that("Rule 2 uses stratum-qualified PSU IDs and every exposed band", {
    x <- fix_fixture()
    x$psu[x$psu == 3] <- 2
    p <- fix_prepare(x)
    b <- .u5mr_target(p, "B", "p1")
    expect_equal(length(unique(b$psu)), 2)
    expect_false(.u5mr_target_rule2(b))
    b$died[b$age0 == "0"] <- 0
    expect_false(.u5mr_target_rule2(b))
    b$died <- b$n * 0.2
    expect_true(.u5mr_target_rule2(b))
    expect_true(.u5mr_target_rule2(.u5mr_target(p, "A", "p1")))
    expect_false(.u5mr_target_rule2(.u5mr_target(p, "A", "empty")))
})

test_that("references pool all areas by default and separate periods", {
    x <- fix_fixture()
    x$region[x$psu == 4] <- NA_character_
    p <- fix_prepare(x)
    expect_identical(unique(p$reference$reference), "All")
    r <- p$reference[p$reference$time0 == "p1" & p$reference$age0 == "0", ]
    expect_equal(r$K, 4)
    expect_equal(r$weight.ph, 2.5)
    expect_equal(r$weighted.exposure, 100)
    expect_equal(r$weighted.deaths, 26)
    expect_equal(r$exposure.ph, 10)
    expect_equal(r$deaths.ph, 2.6)
    r2 <- p$reference[p$reference$time0 == "p2" & p$reference$age0 == "0", ]
    expect_equal(r2$weighted.deaths, 36)
    separate <- fix_prepare(x, "strata")$reference
    expect_setequal(unique(separate$reference), c("s1", "s2"))
    expect_false(isTRUE(all.equal(p$reference, separate)))
    national <- fix_direct(x, var.fix = TRUE, national.only = TRUE)
    expect_equal(fix_table(national), fix_direct(x, national.only = TRUE))
})

test_that("all-PSU K and exposure-weighted PSU weights preserve phantom totals", {
    x <- fix_fixture()
    x <- x[!(x$psu == 4 & x$age == "1-1"), ]
    x$weights[x$psu == 1 & x$age == "1-1"] <- 5
    x$total[x$psu == 1 & x$age == "1-1"] <- 20
    p <- fix_prepare(x)
    r <- p$reference[p$reference$time0 == "p1" & p$reference$age0 == "1-1", ]
    expect_equal(r$K, 4)
    expect_equal(r$weight.ph, ((10 + 100)/30 + 2 + 3 + 4)/4)
    expect_equal(r$weight.ph * r$exposure.ph, (100 + 20 + 30)/4)
    expect_equal(r$weight.ph * r$deaths.ph, (0 + 6 + 18)/4)
})

test_that("augmentation preserves strata, stages and unique synthetic IDs", {
    x <- fix_fixture()
    x$psu <- as.character(x$psu)
    x$psu[x$psu == "1"] <- ".SUMMER.phantom.1.1"
    x$house <- factor(x$house)
    p <- fix_prepare(x, "reference")
    target <- .u5mr_target(p, "B", "p1")
    augmented <- .u5mr_augment(p, target, "B", "p1")
    ph <- augmented$data
    expect_equal(length(unique(ph$psu)), 2)
    expect_equal(nrow(ph), 4)
    expect_setequal(ph$strata, c("s1", "s2"))
    expect_true(all(ph$region0 == "B" & ph$time0 == "p1"))
    expect_length(intersect(ph$stage1, p$data$stage1), 0)
    expect_length(intersect(ph$stage2, p$data$stage2), 0)
    expect_true(all(ph$n > 0 & ph$died >= 0 & ph$died <= ph$n))
    expect_equal(ph$weights0 * ph$died,
        fix_prepare(x, "reference")$reference$weighted.deaths[
            match(paste(ph$time0, ph$reference, ph$age0),
                  with(p$reference, paste(time0, reference, age0)))] / 2)
    c <- .u5mr_augment(p, .u5mr_target(p, "C", "p2"), "C", "p2")$data
    expect_setequal(c$stage1, ph$stage1[ph$strata == "s2"])
    expect_setequal(c$stage2, ph$stage2[ph$strata == "s2"])
    expect_setequal(c$psu, ph$psu[ph$strata == "s2"])
})

test_that("correction agrees for expanded and compact input without changing real domains", {
    x <- fix_fixture()
    before <- x
    expanded <- fix_expand(x)
    compact <- fix_direct(x, var.fix = TRUE)
    long <- fix_direct(expanded, var.fix = TRUE, compact = FALSE)
    expect_equal(compact, long, tolerance = 1e-8)
    expect_identical(x, before)
    legacy <- fix_direct(x)
    keep <- compact$region %in% c("B", "All")
    expect_equal(fix_table(compact)[keep, ], legacy[keep, ], tolerance = 0)
    a <- compact[compact$region == "A" & compact$years == "p1", ]
    expect_true(is.finite(a$logit.est) && is.finite(a$var.est) && a$var.est > 0)
    audit <- attr(compact, "variance.fix")
    expect_true(all(audit$status == "fixed"))
    expect_setequal(attr(compact, "fixed_areas"), c("A", "C"))
    expect_false(anyDuplicated(audit[, c("region", "years")]) > 0)
    expect_true(all(is.na(compact$mean[compact$years == "empty"])))
    expect_identical(fix_direct(x), fix_direct(x, var.fix = FALSE))
    expect_null(attr(legacy, "variance.fix"))
})

test_that("unrelated target triggers do not change a corrected target", {
    x <- fix_fixture()
    y <- x
    y$region[y$region == "B"] <- ifelse(y$psu[y$region == "B"] == 2, "B1", "B2")
    a <- fix_table(fix_direct(x, var.fix = TRUE))
    b <- fix_table(fix_direct(y, var.fix = TRUE))
    expect_equal(a[a$region == "A", ], b[b$region == "A", ], tolerance = 0)
})

test_that("one shared design preserves target-local rows, PSU counts and fits", {
    build.design <- survey::svydesign
    fit.model <- survey::svyglm
    designs <- list()
    calls <- 0L
    testthat::local_mocked_bindings(
        svydesign = function(...) {
            calls <<- calls + 1L
            do.call(build.design, list(...))
        },
        svyglm = function(formula, design, ...) {
            if ("stage1" %in% names(design$variables))
                designs[[length(designs) + 1L]] <<- design
            do.call(fit.model, c(list(formula = formula, design = design), list(...)))
        }, .package = "survey")
    for (setting in c("missing_phantom", "national")) {
        x <- fix_fixture()
        x$psu <- as.character(x$psu)
        x$psu[x$psu == "1"] <- ".SUMMER.phantom.1.1"
        other <- x; other$house <- 2L
        x <- rbind(x, other)
        x$total <- x$total / 2
        x$died <- x$died / 2
        if (setting == "missing_phantom") {
            x$weights[x$strata == "s1" & x$time == "p1"] <- 0
            x$reference[x$time == "p2"] <- "pooled"
            x$region[x$region == "A"] <- "B"
            x$died[x$region == "B" & x$time == "p2"] <-
                x$total[x$region == "B" & x$time == "p2"] * 0.2
        } else {
            x$died <- x$total * 0.2
        }
        designs <- list()
        calls <- 0L
        result <- suppressWarnings(fix_direct(x, var.fix = TRUE, refStrataVar = "reference"))
        expect_equal(calls, 2L)
        p <- fix_prepare(x, "reference")
        audit <- attr(result, "variance.fix")
        expect_length(designs, nrow(audit))
        for (i in seq_len(nrow(audit))) {
            target <- .u5mr_target(p, audit$region[i], audit$years[i])
            ph <- .u5mr_augment(p, target, audit$region[i], audit$years[i])$data
            counts <- rbind(p$data, ph)
            counts$died <- counts$died / counts$n
            # Rebuild the old full-survey design with only this target's
            # phantoms, then select real rows followed by its phantom rows.
            rebuilt <- build.design(ids = p$formula, strata = ~strata,
                nest = TRUE, weights = ~weights0, data = counts)
            rows <- c(.u5mr_target(p, audit$region[i], audit$years[i], rows.only = TRUE),
                nrow(p$data) + seq_len(nrow(ph)))
            rebuilt <- rebuilt[rows, ]
            shared <- designs[[i]]
            expect_equal(unname(shared$fpc$sampsize), unname(rebuilt$fpc$sampsize), tolerance = 0)
            expect_equal(unname(shared$prob), unname(rebuilt$prob), tolerance = 0)
            expect_equal(shared$variables$age0, rebuilt$variables$age0, tolerance = 0)
            expect_equal(shared$variables$died, rebuilt$variables$died, tolerance = 0)
            expected <- suppressWarnings(fit.model(died ~ -1 + factor(age0), design = rebuilt,
                family = quasibinomial(), weights = rebuilt$variables$n))
            actual <- suppressWarnings(fit.model(died ~ -1 + factor(age0), design = shared,
                family = quasibinomial(), weights = shared$variables$n))
            expect_equal(coef(actual), coef(expected), tolerance = 1e-12)
            expect_equal(vcov(actual), vcov(expected), tolerance = 1e-12)
        }
    }
    x <- fix_fixture()
    x$region <- "B"
    calls <- 0L
    fix_direct(x, var.fix = TRUE)
    expect_equal(calls, 1L)
})

test_that("missing area labels still permit national estimation", {
    x <- fix_fixture()
    reference <- fix_direct(x, var.fix = TRUE, national.only = TRUE)
    x$region <- NA_character_
    result <- fix_direct(x, var.fix = TRUE)
    expect_equal(result, reference, tolerance = 0)
    expect_true(all(result$region == "All"))
    expect_true(all(is.finite(result$mean[result$years != "empty"])))
})

test_that("arbitrary field names and factor or character PSU IDs are supported", {
    x <- fix_fixture()
    answer <- fix_direct(x, var.fix = TRUE)
    x$psu <- factor(paste0("cluster", x$psu))
    x$house <- as.character(x$house)
    names(x)[match(c("region", "time", "age", "weights", "psu", "house", "total"), names(x))] <-
        c("district", "period", "band", "weight", "cluster", "household", "exposure")
    renamed <- suppressMessages(getDirect(x, c("p1", "p2", "empty"),
        regionVar = "district", timeVar = "period", ageVar = "band", weightsVar = "weight",
        clusterVar = ~cluster + household, Ntrials = "exposure", var.fix = TRUE))
    expect_equal(answer, renamed, tolerance = 1e-8)
    names(x)[names(x) == "cluster"] <- "PSU ID"
    renamed <- suppressMessages(getDirect(x, c("p1", "p2", "empty"),
        regionVar = "district", timeVar = "period", ageVar = "band", weightsVar = "weight",
        clusterVar = ~`PSU ID` + household, Ntrials = "exposure", var.fix = TRUE))
    expect_equal(answer, renamed, tolerance = 1e-8)
})

test_that("compaction preserves household stages and fractional counts", {
    x <- fix_fixture()
    other <- x; other$house <- 2L
    x <- rbind(x, other)
    x$total <- x$total / 2
    x$died <- x$died / 2
    p <- fix_prepare(x)
    expect_equal(nrow(p$data), nrow(x))
    expect_setequal(p$data$stage2, c("1", "2"))
    expect_equal(sum(p$data$n), sum(x$total))
    expect_equal(sum(p$data$died), sum(x$died))
    expect_no_error(fix_direct(x, var.fix = TRUE))
    one <- x[x$time == "p1" & x$age == "0", ]
    fit <- suppressMessages(getDirect(one, "p1", clusterVar = ~psu + house,
        weightsVar = "weights", Ntrials = "total", var.fix = TRUE, national.only = TRUE))
    expect_equal(nrow(fit), 1)
    expect_true(is.finite(fit$var.est))
})

test_that("degenerate references are unresolved and empty audits have stable types", {
    x <- fix_fixture()
    for (fraction in c(0, 1, 0.2)) {
        x$died <- x$total * fraction
        fit <- suppressWarnings(fix_direct(x, var.fix = TRUE))
        audit <- attr(fit, "variance.fix")
        expect_true(all(audit$status == "unresolved"))
        expect_true(all(audit$reason == "degenerate_reference"))
        expect_true("All" %in% attr(fit, "fixed_areas"))
        if (fraction == 0) expect_true(all(is.na(fit$mean)))
    }
    x <- fix_fixture()
    x$region <- "B"
    fit <- fix_direct(x, var.fix = TRUE)
    expect_identical(attr(fit, "fixed_areas"), character())
    expect_identical(attr(fit, "variance.fix"), .u5mr_fix_audit())
    x$total[x$region == "B"] <- 0
    x$died <- 0
    fit <- fix_direct(x, var.fix = TRUE)
    expect_true(all(is.na(fit$mean)))
    expect_identical(attr(fit, "variance.fix"), .u5mr_fix_audit())
})

test_that("invalid correction inputs fail clearly without affecting the legacy path", {
    x <- fix_fixture()
    x$reference[1] <- "other"
    expect_error(fix_direct(x, var.fix = TRUE, refStrataVar = "reference"), "exactly one reference")
    expect_no_error(fix_direct(x, var.fix = FALSE, refStrataVar = "missing"))
    x$reference[1] <- NA_character_
    expect_error(fix_direct(x, var.fix = TRUE, refStrataVar = "reference"), "nonmissing")
    expect_error(fix_direct(x, var.fix = TRUE, refStrataVar = "missing"), "reference-group column")
    x <- fix_fixture()
    expect_error(.prepare_u5mr_fix(x, ~interaction(psu, house), "total", NULL), "simple named")
    x$total[1] <- -1
    expect_error(fix_direct(x, var.fix = TRUE), "nonnegative exposure")
    x <- fix_fixture(); x$died[1] <- 11
    expect_error(fix_direct(x, var.fix = TRUE), "death counts")
})
