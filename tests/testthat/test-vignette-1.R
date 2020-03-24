
# Setup
library(survey)
data(BRFSS)
data(KingCounty)
BRFSS <- subset(BRFSS, !is.na(BRFSS$diab2))
BRFSS <- subset(BRFSS, !is.na(BRFSS$hracode))
mat <- getAmat(KingCounty, KingCounty$HRA2010v2_)

test_that("Adj matrix is of proper dimension", {
  expect_equal(dim(mat)[1], 48)
  expect_equal(dim(mat)[2], 48)
})

test_that("Adj matrix contains proper names", {
  expect_equal(rownames(mat)[3], "Ballard")
  expect_equal(colnames(mat)[3], "Ballard")
  expect_equal(rownames(mat)[18], "Downtown")
  expect_equal(colnames(mat)[18], "Downtown")
})

test_that("Regions are correctly adjacent/non-adjacent", {
  expect_equal(mat["Auburn-North","Auburn-North"], 0)
  expect_equal(mat["Auburn-North","Auburn-South"], 1)
  expect_equal(mat["Ballard","Downtown"], 0)
  expect_equal(mat["Ballard","NW Seattle"], 1)
  expect_equal(mat["Redmond","NE Seattle"], 0)
  expect_equal(mat["Redmond","Sammamish"], 1)
})

smoothed <- fitGeneric(
  data = BRFSS,
  geo = KingCounty,
  Amat = mat,
  responseType = "binary",
  responseVar = "diab2",
  strataVar = NULL,
  weightVar = NULL,
  regionVar = "hracode",
  clusterVar = NULL,
  CI = 0.95
)

test_that("Object created by fitGeneric is correct", {
  expect_equal(
    signif(as.numeric(smoothed$HT[1,1:8]),4),
    c(-1.813, 0.1727, 0.02983, 33.53,
      0.1403, 0.0004338, 278, 39)
  )
  expect_equal(
    signif(as.numeric(smoothed$smooth[1,3:12]),4),
    c(-1.855, 0.01863, -1.853, -2.125, -1.59,
      0.1362, 0.0002556, 0.1355, 0.1067, 0.1694)
  )
})
