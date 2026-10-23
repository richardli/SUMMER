# Small direct inputs for inspecting formula/data/lincombs without running INLA.
slopes_fixture <- function(){
  W <- matrix(c(0,1,0,1, 1,0,1,0, 0,1,0,1, 1,0,1,0), 4, 4)
  dimnames(W) <- list(LETTERS[1:4], LETTERS[1:4])
  labels <- c("00-02", "03-05", "06-08")
  dat <- expand.grid(region = LETTERS[1:4], years = labels)
  dat$logit.est <- -2 + rep(c(-0.2, 0, 0.2), each = 4)
  dat$logit.prec <- 25
  dat$mean <- plogis(dat$logit.est)
  dat$lower <- plogis(dat$logit.est - 0.4)
  dat$upper <- plogis(dat$logit.est + 0.4)
  dat$var.est <- 0.001
  dat$survey <- "survey"
  list(data = dat, Amat = W, year.label = labels, year.range = c(2000, 2008),
       is.yearly = TRUE, m = 3, time.model = "rw2", st.time.model = "rw1",
       type.st = 4)
}

slopes_mock_inla <- function(.env = parent.frame()){
  skip_if_not_installed("INLA")
  INLA::inla.setOption(num.threads = 1)
  local_mocked_bindings(inla = function(formula, data, lincomb, ...){
    list(.args = list(formula = formula, data = data, lincomb = lincomb))
  }, inla.pc.rprec = function(n, u, alpha) 1 / u^2,
  .package = "INLA", .env = .env)
  # Controlled draws make the copied calibration stop at U = 1 for alpha = .01.
  local_mocked_bindings(rnorm = function(n, mean = 0, sd = 1){
    c(rep(2, round(sd)), rep(0, 100 - round(sd)))
  }, .package = "stats", .env = .env)
  local_mocked_bindings(simhyper = function(...) list(a.iid = 1, b.iid = 1),
                        .package = "SUMMER", .env = .env)
}

test_that("slope calibration is copied from smoothCluster", {
  slopes_mock_inla()
  fit <- suppressMessages(do.call(smoothDirect, c(slopes_fixture(), list(
    pc.st.slope.u = 1, pc.st.slope.alpha = 0.01))))
  env <- environment(fit$model)
  expect_equal(get("pc.st.slope.prec.u", env), 1)
  expect_identical(get("hyperpc.slope", env),
                   list(prec = list(prior = "pc.prec", param = c(1, 0.01))))
  calibration <- function(fun){
    code <- deparse(body(fun))
    first <- grep("rescale_U <- function", code, fixed = TRUE)
    last <- grep("return(ustar)", code, fixed = TRUE)
    parse(text = code[first:(last + 1L)], keep.source = FALSE)
  }
  expect_identical(calibration(smoothDirect), calibration(smoothCluster))
})

test_that("enabled yearly slopes use annual values and period averages", {
  slopes_mock_inla()
  args <- slopes_fixture()
  fit <- suppressMessages(do.call(smoothDirect, c(args, list(
    pc.st.slope.u = 1, pc.st.slope.alpha = 0.01))))
  expected <- c(seq(-0.5, 0.5, length.out = 9), c(-0.375, 0, 0.375))
  expect_equal(fit$newdata$st.slope, expected[fit$newdata$time.struct])
  expect_equal(fit$newdata$st.slope.id, fit$newdata$region.struct)
  expect_match(paste(deparse(fit$model), collapse = " "), "f\\(st.slope.id, st.slope")
  expect_length(fit$fit$.args$lincomb, 48)
  for(k in seq_along(fit$fit$.args$lincomb)){
    terms <- unlist(fit$fit$.args$lincomb[[k]], recursive = FALSE)
    term <- terms$st.slope.id
    t <- fit$lincombs.info$Year[k]
    if(expected[t] == 0){
      expect_length(term$weight, 0)
    }else{
      expect_equal(term$idx, fit$lincombs.info$District[k])
      expect_equal(term$weight, expected[t])
    }
  }
})

test_that("default and excluded RW2 slopes leave model inputs unchanged", {
  slopes_mock_inla()
  local_mocked_bindings(inla.pc.rprec = function(...) stop("prior was computed"),
                        .package = "INLA")
  args <- slopes_fixture()
  baseline <- suppressMessages(do.call(smoothDirect, args))
  expect_false("st.slope" %in% names(baseline$newdata))
  expect_false(grepl("st.slope", paste(deparse(baseline$model), collapse = " ")))
  for(main in c("rw1", "rw2")){
    args$time.model <- main
    args$st.time.model <- "rw2"
    set.seed(123)
    base <- suppressMessages(do.call(smoothDirect, args))
    seed <- .Random.seed
    set.seed(123)
    messages <- character()
    excluded <- withCallingHandlers(do.call(smoothDirect, c(args, list(
      pc.st.slope.u = 1, pc.st.slope.alpha = 0.01))), message = function(m){
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
    expect_identical(.Random.seed, seed)
    expect_equal(sum(grepl("slopes excluded", messages)), 1)
    expect_identical(excluded$newdata, base$newdata)
    expect_identical(deparse(excluded$model), deparse(base$model))
    expect_identical(excluded$fit$.args$lincomb, base$fit$.args$lincomb)
    expect_identical(get("constr.st", environment(excluded$model)),
                     get("constr.st", environment(base$model)))
  }
  # The inherited effective interaction is RW2 as well.
  args$st.time.model <- NULL
  expect_message(do.call(smoothDirect, c(args, list(pc.st.slope.u = 1))),
                 "slopes excluded")
})

test_that("unsupported or incomplete enabled slopes fail clearly", {
  args <- slopes_fixture()
  for(prior in list(list(pc.st.slope.u = 1), list(pc.st.slope.alpha = 0.1),
                    list(pc.st.slope.u = 0, pc.st.slope.alpha = 0.1),
                    list(pc.st.slope.u = Inf, pc.st.slope.alpha = 0.1),
                    list(pc.st.slope.u = 1, pc.st.slope.alpha = 1),
                    list(pc.st.slope.u = c(1, 2), pc.st.slope.alpha = 0.1))){
    expect_error(suppressMessages(do.call(smoothDirect, c(args, prior))),
                 "Supply both pc.st.slope")
  }
  args["Amat"] <- list(NULL)
  expect_error(suppressMessages(do.call(smoothDirect, c(args,
    list(pc.st.slope.u = 1, pc.st.slope.alpha = 0.1)))), "spatial and temporal")
  args <- slopes_fixture()
  args["time.model"] <- list(NULL)
  expect_error(suppressMessages(do.call(smoothDirect, c(args,
    list(pc.st.slope.u = 1, pc.st.slope.alpha = 0.1)))), "spatial and temporal")
})

test_that("ordinary slopes use correct grouped and Type II prediction indices", {
  slopes_mock_inla()
  args <- slopes_fixture()
  args$is.yearly <- FALSE
  args$type.st <- 2
  for(setting in list(c("rw1", "pc"), c("ar1", "pc"), c("rw1", "gamma"))){
    interaction <- setting[1]
    args$st.time.model <- interaction
    args$hyper <- setting[2]
    fit <- suppressMessages(do.call(smoothDirect, c(args, list(
      pc.st.slope.u = 1, pc.st.slope.alpha = 0.01))))
    expect_equal(fit$newdata$st.slope, c(-0.5, 0, 0.5)[fit$newdata$time.struct])
    for(k in seq_along(fit$fit$.args$lincomb)){
      terms <- unlist(fit$fit$.args$lincomb[[k]], recursive = FALSE)
      j <- fit$lincombs.info$District[k]
      t <- fit$lincombs.info$Year[k]
      if(interaction == "ar1"){
        expect_equal(terms$region.int$idx, j + (t - 1) * 4)
      }else{
        expect_equal(terms$time.area$idx, t + (j - 1) * 3)
      }
    }
  }
})

test_that("FH joint summaries include signed area slopes and zero midpoints", {
  skip_if_not_installed("INLA")
  fields <- c("(Intercept):1", "region.struct:1", "region.struct:2",
              "st.slope.id:1", "st.slope.id:2", "time.slope:1")
  latent <- cbind(c(-2, 0.2, -0.1, 0.4, -0.6, 0.3),
                  c(-1.8, 0.1, -0.2, 0.2, -0.4, 0.1),
                  c(-2.1, 0.3, 0.1, 0.6, -0.2, -0.1))
  rownames(latent) <- fields
  draws <- lapply(1:3, function(i) list(latent = latent[, i, drop = FALSE]))
  lc <- baseline <- list()
  expected <- expected.base <- matrix(0, 3, 6)
  for(j in 1:2){
    for(t in 1:3){
      k <- 3 * (j - 1) + t
      w <- c(-0.5, 0, 0.5)[t]
      area <- slope <- c(NA, NA)
      area[j] <- 1
      slope[j] <- w
      terms <- list("(Intercept)" = 1, region.struct = area, time.slope = w)
      baseline[[k]] <- INLA::inla.make.lincomb(terms)[[1]]
      if(w != 0) terms$st.slope.id <- slope
      lc[[k]] <- INLA::inla.make.lincomb(terms)[[1]]
      eta <- latent[1, ] + latent[1 + j, ] + w * latent[6, ]
      expected.base[, j + (t - 1) * 2] <- plogis(eta)
      expected[, j + (t - 1) * 2] <- plogis(eta + w * latent[3 + j, ])
    }
  }
  W <- matrix(c(0, 1, 1, 0), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  fit <- list(fit = list(.args = list(lincomb = lc),
                         summary.random = list(region.struct = data.frame()),
                         misc = list(configs = list(contents = data.frame(tag = "region.struct")))),
              Amat = W, has.Amat = TRUE, is.yearly = FALSE, is.temporal = TRUE,
              year.range = c(2000, 2002), year.label = as.character(2000:2002),
              lincombs.info = data.frame(Index = 1:6, District = rep(1:2, each = 3), Year = rep(1:3, 2)))
  out <- suppressMessages(getSmoothed(fit, joint = TRUE, draws = draws))
  expect_equal(out$overall$mean, colMeans(expected))
  expect_equal(out$overall$variance, apply(expected, 2, var))
  expect_equal(do.call(cbind, lapply(out$draws.est.overall, `[[`, "draws")), expected)
  fit$fit$.args$lincomb <- baseline
  out <- suppressMessages(getSmoothed(fit, joint = TRUE, draws = draws))
  expect_equal(out$overall$mean, colMeans(expected.base))
  expect_equal(out$overall$variance, apply(expected.base, 2, var))
  expect_equal(do.call(cbind, lapply(out$draws.est.overall, `[[`, "draws")), expected.base)
})
