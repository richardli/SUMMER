
# The following functions are tested in this file:
#     expit()
#     logit()

test_that("Expit works", {
  expect_equal(expit(2), exp(2)/(1 + exp(2)))
  expect_equal(expit(-2), exp(-2)/(1 + exp(-2)))
})

test_that("Logit works", {
  expect_equal(logit(0.2), log(0.2/(1-0.2)))
  expect_equal(logit(0.8), log(0.8/(1-0.8)))
})
