test_that("Type IV constraints are independent by connected component", {
  skip_if_not_installed("INLA")

  W <- matrix(0, 6, 6)
  W[cbind(c(1, 2, 2, 3, 4, 5), c(2, 1, 3, 2, 5, 4))] <- 1
  setup <- SUMMER:::.type4_inla_components(
    Amat = W, n.time = 7, rw.order = 2, constr = TRUE
  )

  # Two non-singleton components use one reference region apiece. The
  # singleton is not a spatial component constraint and retains its temporal
  # sum constraint.
  expect_equal(setup$reference, c(1L, 4L))
  expect_equal(setup$rankdef, 22)
  expect_equal(nrow(setup$extraconstr$A), 18)
  expect_equal(qr(setup$extraconstr$A)$rank, 18)
})

test_that("Type IV setup accepts row-standardized undirected adjacency", {
  skip_if_not_installed("INLA")

  W <- matrix(0, 5, 5)
  W[1, 2:5] <- 1
  W[2:5, 1] <- 1
  W <- W / rowSums(W)

  expect_silent(SUMMER:::.type4_inla_components(
    Amat = W, n.time = 7, rw.order = 1, constr = TRUE
  ))
})

test_that("Type II uses independent ANOVA constraints", {
  skip_if_not_installed("INLA")

  setup <- SUMMER:::.type2_inla_components(
    n.region = 4, n.time = 6, rw.order = 2, constr = TRUE
  )

  expect_equal(setup$rankdef, 8)
  expect_equal(nrow(setup$extraconstr$A), 9)
  expect_equal(qr(setup$extraconstr$A)$rank, 9)
  expect_equal(setup$reference, 1L)
})

test_that("smoothSurvey uses smoothDirect Type II and III labels", {
  skip_if_not_installed("INLA")

  regions <- letters[1:3]
  W <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), 3, 3,
              dimnames = list(regions, regions))
  direct <- expand.grid(region = regions, time = 1:3)
  direct$estimate <- seq_len(nrow(direct)) / 10
  direct$variance <- 0.1

  make_formula <- function(type) {
    fit <- smoothSurvey(
      data = NULL, direct.est = direct, Amat = W,
      regionVar = "region", timeVar = "time", responseVar = "estimate",
      direct.est.var = "variance", response.type = "gaussian",
      type.st = type, smooth = FALSE
    )
    paste(deparse(fit$formula), collapse = " ")
  }

  type2 <- make_formula(2)
  type3 <- make_formula(3)
  expect_match(
    type2,
    'f\\(time.area, model = "generic0".*Cmatrix = R.st'
  )
  expect_match(
    type3,
    'f\\(region.int, model = "besag".*control.group = list\\(model = "iid"'
  )
})
