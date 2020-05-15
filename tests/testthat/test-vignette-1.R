# This set of tests corresponds to the "Bayesian Small Area Estimation using Complex Survey Data" Vignette
# https://cran.r-project.org/web/packages/SUMMER/vignettes/summer-vignette.html

# The following functions are tested in this file:
#     getAmat()
#     fitGeneric()
#     getDirectList()
#     aggregateSurvey()
#     fitINLA()
#     getSmoothed()
#     rst()

# Define custom expectation function; expectation is satisfied if the percent
#    difference between `x` and `y` is no more than 10%
expect_10pd <- function(x, y) {
  percent_difference <- abs(x-y) / ((x+y)/2)
  return( expect_lte(percent_difference, 0.1) )
}

# The following tests correspond to the "Small Area Estimation with BRFSS Data"
#     vignette section

# Run code

library(survey)
data(BRFSS)
data(KingCounty)
BRFSS <- subset(BRFSS, !is.na(BRFSS$diab2))
BRFSS <- subset(BRFSS, !is.na(BRFSS$hracode))
mat <- getAmat(KingCounty, KingCounty$HRA2010v2_)
# !!!!! Subset data here
smoothed <- fitGeneric(data=BRFSS, geo=KingCounty, Amat=mat,
                       responseType="binary", responseVar="diab2",
                       strataVar=NULL, weightVar=NULL, regionVar="hracode",
                       clusterVar=NULL, CI=0.95)
svysmoothed <- fitGeneric(data=BRFSS, geo=KingCounty, Amat=mat,
                          responseType="binary", responseVar="diab2",
                          strataVar="strata", weightVar="rwt_llcp",
                          regionVar="hracode", clusterVar="~1", CI=0.95)
design <- svydesign(ids=~1, weights=~rwt_llcp, strata=~strata, data=BRFSS)
Education <- svyby(~educau, ~hracode, design, svymean)
Education <- Education[,1:5]
svysmoothed_w_cov <- fitGeneric(data=BRFSS, geo=KingCounty, Amat=mat,
                                X=Education, responseType="binary",
                                responseVar="diab2", strataVar="strata",
                                weightVar="rwt_llcp", regionVar="hracode",
                                clusterVar="~1", CI=0.95)
svysmoothed_w_cov$fit$summary.fixed

# Run tests

test_that("mat, created by getAmat(), is of proper dimension", {
  expect_equal(dim(mat)[1], c(48,48))
})

test_that("mat, created by getAmat(), contains correct names", {
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

test_that("smoothed$smooth, created by fitGeneric(), is correct", {
  expect_10pd( smoothed$smooth$mean[1], 1 )
  expect_10pd( smoothed$smooth$mean[2], 1 )
  expect_10pd( smoothed$smooth$variance[1], 1 )
  expect_10pd( smoothed$smooth$variance[2], 1 )
  expect_10pd( smoothed$smooth$median[1], 1 )
  expect_10pd( smoothed$smooth$median[2], 1 )
  expect_10pd( smoothed$smooth$lower[1], 1 )
  expect_10pd( smoothed$smooth$lower[2], 1 )
  expect_10pd( smoothed$smooth$upper[1], 1 )
  expect_10pd( smoothed$smooth$upper[2], 1 )
  expect_10pd( smoothed$smooth$mean.original[1], 1 )
  expect_10pd( smoothed$smooth$mean.original[2], 1 )
  expect_10pd( smoothed$smooth$variance.original[1], 1 )
  expect_10pd( smoothed$smooth$variance.original[2], 1 )
  expect_10pd( smoothed$smooth$median.original[1], 1 )
  expect_10pd( smoothed$smooth$median.original[2], 1 )
  expect_10pd( smoothed$smooth$lower.original[1], 1 )
  expect_10pd( smoothed$smooth$lower.original[2], 1 )
  expect_10pd( smoothed$smooth$upper.original[1], 1 )
  expect_10pd( smoothed$smooth$upper.original[2], 1 )
})

test_that("smoothed$HT, created by fitGeneric(), is correct", {
  expect_10pd( smoothed$HT$HT.est[1], 1 )
  expect_10pd( smoothed$HT$HT.est[2], 1 )
  expect_10pd( smoothed$HT$HT.sd[1], 1 )
  expect_10pd( smoothed$HT$HT.sd[2], 1 )
  expect_10pd( smoothed$HT$HT.variance[1], 1 )
  expect_10pd( smoothed$HT$HT.variance[2], 1 )
  expect_10pd( smoothed$HT$HT.prec[1], 1 )
  expect_10pd( smoothed$HT$HT.prec[2], 1 )
  expect_10pd( smoothed$HT$HT.est.original[1], 1 )
  expect_10pd( smoothed$HT$HT.est.original[2], 1 )
  expect_10pd( smoothed$HT$HT.variance.original[1], 1 )
  expect_10pd( smoothed$HT$HT.variance.original[2], 1 )
  expect_10pd( smoothed$HT$n[1], 1 )
  expect_10pd( smoothed$HT$n[2], 1 )
  expect_10pd( smoothed$HT$y[1], 1 )
  expect_10pd( smoothed$HT$y[2], 1 )
})

test_that("svysmoothed$smooth, created by fitGeneric(), is correct", {
  expect_10pd( svysmoothed$smooth$mean[1], 1 )
  expect_10pd( svysmoothed$smooth$mean[2], 1 )
  expect_10pd( svysmoothed$smooth$variance[1], 1 )
  expect_10pd( svysmoothed$smooth$variance[2], 1 )
  expect_10pd( svysmoothed$smooth$median[1], 1 )
  expect_10pd( svysmoothed$smooth$median[2], 1 )
  expect_10pd( svysmoothed$smooth$lower[1], 1 )
  expect_10pd( svysmoothed$smooth$lower[2], 1 )
  expect_10pd( svysmoothed$smooth$upper[1], 1 )
  expect_10pd( svysmoothed$smooth$upper[2], 1 )
  expect_10pd( svysmoothed$smooth$mean.original[1], 1 )
  expect_10pd( svysmoothed$smooth$mean.original[2], 1 )
  expect_10pd( svysmoothed$smooth$variance.original[1], 1 )
  expect_10pd( svysmoothed$smooth$variance.original[2], 1 )
  expect_10pd( svysmoothed$smooth$median.original[1], 1 )
  expect_10pd( svysmoothed$smooth$median.original[2], 1 )
  expect_10pd( svysmoothed$smooth$lower.original[1], 1 )
  expect_10pd( svysmoothed$smooth$lower.original[2], 1 )
  expect_10pd( svysmoothed$smooth$upper.original[1], 1 )
  expect_10pd( svysmoothed$smooth$upper.original[2], 1 )
})

test_that("svysmoothed$HT, created by fitGeneric(), is correct", {
  expect_10pd( svysmoothed$HT$HT.est[1], 1 )
  expect_10pd( svysmoothed$HT$HT.est[2], 1 )
  expect_10pd( svysmoothed$HT$HT.sd[1], 1 )
  expect_10pd( svysmoothed$HT$HT.sd[2], 1 )
  expect_10pd( svysmoothed$HT$HT.variance[1], 1 )
  expect_10pd( svysmoothed$HT$HT.variance[2], 1 )
  expect_10pd( svysmoothed$HT$HT.prec[1], 1 )
  expect_10pd( svysmoothed$HT$HT.prec[2], 1 )
  expect_10pd( svysmoothed$HT$HT.est.original[1], 1 )
  expect_10pd( svysmoothed$HT$HT.est.original[2], 1 )
  expect_10pd( svysmoothed$HT$HT.variance.original[1], 1 )
  expect_10pd( svysmoothed$HT$HT.variance.original[2], 1 )
  expect_10pd( svysmoothed$HT$n[1], 1 )
  expect_10pd( svysmoothed$HT$n[2], 1 )
  expect_10pd( svysmoothed$HT$y[1], 1 )
  expect_10pd( svysmoothed$HT$y[2], 1 )
})

test_that("svysmoothed_w_cov$fit, created by fitGeneric(), is correct", {
  expect_10pd( svysmoothed_w_cov$fit$summary.fixed[1,"mean"], 1 )
  expect_10pd( svysmoothed_w_cov$fit$summary.fixed[2,"mean"], 1 )
  expect_10pd( svysmoothed_w_cov$fit$summary.fixed[1,"sd"], 1 )
  expect_10pd( svysmoothed_w_cov$fit$summary.fixed[2,"sd"], 1 )
  expect_10pd( svysmoothed_w_cov$fit$summary.fixed[1,"0.025quant"], 1 )
  expect_10pd( svysmoothed_w_cov$fit$summary.fixed[2,"0.025quant"], 1 )
  expect_10pd( svysmoothed_w_cov$fit$summary.fixed[1,"0.5quant"], 1 )
  expect_10pd( svysmoothed_w_cov$fit$summary.fixed[2,"0.5quant"], 1 )
})

# The following tests correspond to the "U5MR Estimation in Space and Time"
#     vignette section

# Run code

data(DemoData)
data(DemoMap)
geo <- DemoMap$geo
mat <- DemoMap$Amat
years <- levels(DemoData[[1]]$time)
data_multi <- getDirectList(births=DemoData, years=years, regionVar="region",
                            timeVar="time", clusterVar="~clustid+id",
                            ageVar="age", weightsVar="weights", geo.recode=NULL)
data <- aggregateSurvey(data_multi)
years.all <- c(years, "15-19")
fit1 <- fitINLA(data=data, geo=NULL, Amat=NULL, year_label=years.all,
                year_range=c(1985, 2019), rw=2, is.yearly=FALSE, m=5)

fit2 <- fitINLA(data=data, geo=NULL, Amat=NULL, year_label=years.all,
                year_range=c(1985,2019), rw=2, is.yearly=TRUE, m=5, type.st=4)
out1 <- getSmoothed(fit1)
out2 <- getSmoothed(fit2)
fit3 <- fitINLA(data=data, geo=geo, Amat=mat, year_label=years.all,
                year_range=c(1985,2019), rw=2, is.yearly=FALSE, type.st=4)
out3 <- getSmoothed(fit3, Amat=mat)
fit4 <- fitINLA(data=data, geo=geo, Amat=mat, year_label=years.all,
                year_range=c(1985,2019), rw=2, is.yearly=TRUE, m=5, type.st=4)
out4 <- getSmoothed(fit4, Amat=mat)

# Run tests

test_that("data_multi, created by getDirectList(), is of proper dimension", {
  expect_equal(dim(data_multi), c(150,11))
})

test_that("data, created by aggregateSurvey(), is of proper dimension", {
  expect_equal(dim(data), c(30,10))
})

test_that("out1, created by getSmoothed(), is correct", {
  expect_equal( out1[1,"years"], "85-89" )
  expect_equal( out1[2,"years"], "90-94" )
  expect_10pd( out1[1,"logit.upper"], 1 )
  expect_10pd( out1[2,"logit.upper"], 1 )
  expect_10pd( out1[1,"logit.lower"], 1 )
  expect_10pd( out1[2,"logit.lower"], 1 )
  expect_10pd( out1[1,"logit.median"], 1 )
  expect_10pd( out1[2,"logit.median"], 1 )
  expect_10pd( out1[1,"upper"], 1 )
  expect_10pd( out1[2,"upper"], 1 )
  expect_10pd( out1[1,"lower"], 1 )
  expect_10pd( out1[2,"lower"], 1 )
  expect_10pd( out1[1,"median"], 1 )
  expect_10pd( out1[2,"median"], 1 )
})

test_that("out2, created by getSmoothed(), is correct", {
  expect_equal( out2[1,"years"], "1985" )
  expect_equal( out2[2,"years"], "1986" )
  expect_10pd( out2[1,"logit.upper"], 1 )
  expect_10pd( out2[2,"logit.upper"], 1 )
  expect_10pd( out2[1,"logit.lower"], 1 )
  expect_10pd( out2[2,"logit.lower"], 1 )
  expect_10pd( out2[1,"logit.median"], 1 )
  expect_10pd( out2[2,"logit.median"], 1 )
  expect_10pd( out2[1,"upper"], 1 )
  expect_10pd( out2[2,"upper"], 1 )
  expect_10pd( out2[1,"lower"], 1 )
  expect_10pd( out2[2,"lower"], 1 )
  expect_10pd( out2[1,"median"], 1 )
  expect_10pd( out2[2,"median"], 1 )
})

# !!!!! Continute

# !!!!! Same for out3, out4

# The following tests correspond to the "Simulate spatial(temporal) random
#     effects" vignette section

data(KingCounty)
u <- rst(n = 3, type = "s", Amat = getAmat(KingCounty, names = KingCounty$HRA2010v2_))
mapPlot(data.frame(r1 = u[1, ], r2 = u[2, ], r3 = u[3, ], region = colnames(u)), 
        geo = KingCounty, by.data = "region", by.geo = "HRA2010v2_", variables = c("r1", 
                                                                                   "r2", "r3"))

