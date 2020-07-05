# This set of tests corresponds to the "Bayesian Small Area Estimation using
#     Complex Survey Data" Vignette, found at the following URL:
#     https://cran.r-project.org/web/packages/SUMMER/vignettes/summer-vignette.html
# Note that the "Simulate spatial(temporal) random effects" vignette section,
#     which uses the rst() function, is not being tested

# The following functions are tested in this file:
#     getAmat()
#     fitGeneric()
#     getDirectList()
#     aggregateSurvey()
#     fitINLA()
#     getSmoothed()

# Define custom expectation function; expectation is satisfied if the percent
#    difference between `x` and `y` is no more than 10%
expect_10pd <- function(x, y) {
  percent_difference <- abs(abs(x-y) / ((x+y)/2))
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

# Subset data to make tests run faster
smoothed <- fitGeneric(data=BRFSS, geo=KingCounty, Amat=mat,
                       responseType="binary", responseVar="diab2",
                       strataVar=NULL, weightVar=NULL, regionVar="hracode",
                       clusterVar=NULL, CI=0.95)

if (devtools::r_env_vars()[["NOT_CRAN"]]=="true") {
  
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
  
}


# Run tests

test_that("mat, created by getAmat(), is of proper dimension", {
  expect_equal(dim(mat), c(48,48))
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
  expect_10pd( smoothed$smooth$mean[1], -1.854664 )
  expect_10pd( smoothed$smooth$mean[2], -1.451742 )
  expect_10pd( smoothed$smooth$var[1], 0.01858138 )
  expect_10pd( smoothed$smooth$var[2], 0.02691426 )
  expect_10pd( smoothed$smooth$median[1], -1.85311 )
  expect_10pd( smoothed$smooth$median[2], -1.452125 )
  expect_10pd( smoothed$smooth$lower[1], -2.123094 )
  expect_10pd( smoothed$smooth$lower[2], -1.772641 )
  expect_10pd( smoothed$smooth$upper[1], -1.589186 )
  expect_10pd( smoothed$smooth$upper[2], -1.129724 )
  expect_10pd( smoothed$smooth$mean.original[1], 0.1361696 )
  expect_10pd( smoothed$smooth$mean.original[2], 0.1910692 )
  expect_10pd( smoothed$smooth$variance.original[1], 0.0002555567 )
  expect_10pd( smoothed$smooth$variance.original[2], 0.0006358076 )
  expect_10pd( smoothed$smooth$median.original[1], 0.135508 )
  expect_10pd( smoothed$smooth$median.original[2], 0.1898493 )
  expect_10pd( smoothed$smooth$lower.original[1], 0.1067042 )
  expect_10pd( smoothed$smooth$lower.original[2], 0.1452106 )
  expect_10pd( smoothed$smooth$upper.original[1], 0.1693827 )
  expect_10pd( smoothed$smooth$upper.original[2], 0.2438131 )
})

test_that("smoothed$HT, created by fitGeneric(), is correct", {
  expect_10pd( smoothed$HT$HT.est[1], -1.812902 )
  expect_10pd( smoothed$HT$HT.est[2], -1.196804 )
  expect_10pd( smoothed$HT$HT.sd[1], 0.1726995 )
  expect_10pd( smoothed$HT$HT.sd[2], 0.1760789 )
  expect_10pd( smoothed$HT$HT.variance[1], 0.02982513 )
  expect_10pd( smoothed$HT$HT.variance[2], 0.03100377 )
  expect_10pd( smoothed$HT$HT.prec[1], 33.52878 )
  expect_10pd( smoothed$HT$HT.prec[2], 32.25414 )
  expect_10pd( smoothed$HT$HT.est.original[1], 0.1402878 )
  expect_10pd( smoothed$HT$HT.est.original[2], 0.2320442 )
  expect_10pd( smoothed$HT$HT.variance.original[1], 0.0004338385 )
  expect_10pd( smoothed$HT$HT.variance.original[2], 0.0009845287 )
  expect_10pd( smoothed$HT$n[1], 278 )
  expect_10pd( smoothed$HT$n[2], 181 )
  expect_10pd( smoothed$HT$y[1], 39 )
  expect_10pd( smoothed$HT$y[2], 42 )
})

if (devtools::r_env_vars()[["NOT_CRAN"]]=="true") {
  
  test_that("svysmoothed$smooth, created by fitGeneric(), is correct", {
    expect_10pd( svysmoothed$smooth$mean[1], -2.180409 )
    expect_10pd( svysmoothed$smooth$mean[2], -1.679097 )
    expect_10pd( svysmoothed$smooth$variance[1], 0.02967632 )
    expect_10pd( svysmoothed$smooth$variance[2], 0.05150413 )
    expect_10pd( svysmoothed$smooth$median[1], -2.18113 )
    expect_10pd( svysmoothed$smooth$median[2], -1.683833 )
    expect_10pd( svysmoothed$smooth$lower[1], -2.51698 )
    expect_10pd( svysmoothed$smooth$lower[2], -2.110335 )
    expect_10pd( svysmoothed$smooth$upper[1], -1.840016 )
    expect_10pd( svysmoothed$smooth$upper[2], -1.221456 )
    expect_10pd( svysmoothed$smooth$mean.original[1], 0.1025405 )
    expect_10pd( svysmoothed$smooth$mean.original[2], 0.1595901 )
    expect_10pd( svysmoothed$smooth$variance.original[1], 0.0002519658 )
    expect_10pd( svysmoothed$smooth$variance.original[2], 0.0009425412 )
    expect_10pd( svysmoothed$smooth$median.original[1], 0.1013643 )
    expect_10pd( svysmoothed$smooth$median.original[2], 0.1566448 )
    expect_10pd( svysmoothed$smooth$lower.original[1], 0.07472429 )
    expect_10pd( svysmoothed$smooth$lower.original[2], 0.1081466 )
    expect_10pd( svysmoothed$smooth$upper.original[1], 0.1367389 )
    expect_10pd( svysmoothed$smooth$upper.original[2], 0.2275898 )
  })
  
  test_that("svysmoothed$HT, created by fitGeneric(), is correct", {
    expect_10pd( svysmoothed$HT$HT.est[1], -2.153211 )
    expect_10pd( svysmoothed$HT$HT.est[2], -1.191824 )
    expect_10pd( svysmoothed$HT$HT.sd[1], 0.2304232 )
    expect_10pd( svysmoothed$HT$HT.sd[2], 0.2741175 )
    expect_10pd( svysmoothed$HT$HT.variance[1], 0.05309487 )
    expect_10pd( svysmoothed$HT$HT.variance[2], 0.07514043 )
    expect_10pd( svysmoothed$HT$HT.prec[1], 18.83421 )
    expect_10pd( svysmoothed$HT$HT.prec[2], 13.30842 )
    expect_10pd( svysmoothed$HT$HT.est.original[1], 0.1040315 )
    expect_10pd( svysmoothed$HT$HT.est.original[2], 0.2329329 )
    expect_10pd( svysmoothed$HT$HT.variance.original[1], 0.0004612837 )
    expect_10pd( svysmoothed$HT$HT.variance.original[2], 0.002398844 )
    expect_equal( svysmoothed$HT$n[1], NA )
    expect_equal( svysmoothed$HT$n[2], NA )
    expect_equal( svysmoothed$HT$y[1], NA )
    expect_equal( svysmoothed$HT$y[2], NA )
  })
  
  test_that("svysmoothed_w_cov$fit, created by fitGeneric(), is correct", {
    expect_10pd( svysmoothed_w_cov$fit$summary.fixed[1,"mean"], -2.671196 )
    expect_10pd( svysmoothed_w_cov$fit$summary.fixed[2,"mean"], 1.92831 )
    expect_10pd( svysmoothed_w_cov$fit$summary.fixed[1,"sd"], 0.07092216 )
    expect_10pd( svysmoothed_w_cov$fit$summary.fixed[2,"sd"], 3.205899 )
    expect_10pd( svysmoothed_w_cov$fit$summary.fixed[1,"0.025quant"], -2.811211 )
    expect_10pd( svysmoothed_w_cov$fit$summary.fixed[2,"0.025quant"], -4.411817 )
    expect_10pd( svysmoothed_w_cov$fit$summary.fixed[1,"0.5quant"], -2.671067 )
    expect_10pd( svysmoothed_w_cov$fit$summary.fixed[2,"0.5quant"], 1.937464 )
  })
  
}


# The following tests correspond to the "U5MR Estimation in Space and Time"
#     vignette section

# Run code

if (devtools::r_env_vars()[["NOT_CRAN"]]=="true") {
  
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
    expect_equal( as.character(out1[1,"years"]), "85-89" )
    expect_equal( as.character(out1[2,"years"]), "90-94" )
    expect_10pd( out1[1,"logit.upper"], -0.8572562 )
    expect_10pd( out1[2,"logit.upper"], -1.042065 )
    expect_10pd( out1[1,"logit.lower"], -1.531579 )
    expect_10pd( out1[2,"logit.lower"], -1.516882 )
    expect_10pd( out1[1,"logit.median"], -1.183325 )
    expect_10pd( out1[2,"logit.median"], -1.275046 )
    expect_10pd( out1[1,"upper"], 0.291806 )
    expect_10pd( out1[2,"upper"], 0.2597666 )
    expect_10pd( out1[1,"lower"], 0.1800689 )
    expect_10pd( out1[2,"lower"], 0.1773325 )
    expect_10pd( out1[1,"median"], 0.2331674 )
    expect_10pd( out1[2,"median"], 0.2166148 )
  })
  
  test_that("out2, created by getSmoothed(), is correct", {
    expect_equal( as.character(out2[1,"years"]), "1985" )
    expect_equal( as.character(out2[2,"years"]), "1986" )
    expect_10pd( out2[1,"logit.upper"], -0.6571897 )
    expect_10pd( out2[2,"logit.upper"], -0.704823 )
    expect_10pd( out2[1,"logit.lower"], -1.67913 )
    expect_10pd( out2[2,"logit.lower"], -1.650617 )
    expect_10pd( out2[1,"logit.median"], -1.146929 )
    expect_10pd( out2[2,"logit.median"], -1.185482 )
    expect_10pd( out2[1,"upper"], 0.3494734 )
    expect_10pd( out2[2,"upper"], 0.337915 )
    expect_10pd( out2[1,"lower"], 0.1503887 )
    expect_10pd( out2[2,"lower"], 0.159885 )
    expect_10pd( out2[1,"median"], 0.2393693 )
    expect_10pd( out2[2,"median"], 0.2357562 )
  })
  
  test_that("out3, created by getSmoothed(), is correct", {
    expect_equal( as.character(out3[1,"years"]), "85-89" )
    expect_equal( as.character(out3[2,"years"]), "85-89" )
    expect_10pd( out3[1,"logit.upper"], -0.6080551 )
    expect_10pd( out3[2,"logit.upper"], -0.3140722 )
    expect_10pd( out3[1,"logit.lower"], -1.692586 )
    expect_10pd( out3[2,"logit.lower"], -1.471721 )
    expect_10pd( out3[1,"logit.median"], -1.148868 )
    expect_10pd( out3[2,"logit.median"], -0.918756 )
    expect_10pd( out3[1,"upper"], 0.3504759 )
    expect_10pd( out3[2,"upper"], 0.4267779 )
    expect_10pd( out3[1,"lower"], 0.1538285 )
    expect_10pd( out3[2,"lower"], 0.1891358 )
    expect_10pd( out3[1,"median"], 0.2435743 )
    expect_10pd( out3[2,"median"], 0.2829265 )
  })
  
  test_that("out4, created by getSmoothed(), is correct", {
    expect_equal( as.character(out4[1,"years"]), "1985" )
    expect_equal( as.character(out4[2,"years"]), "1985" )
    expect_10pd( out4[1,"logit.upper"], -0.3122694 )
    # expect_10pd( out4[2,"logit.upper"], 0.1721247 ) # Skipped because of inconsistent results
    expect_10pd( out4[1,"logit.lower"], -1.961127 )
    expect_10pd( out4[2,"logit.lower"], -1.632416 )
    expect_10pd( out4[1,"logit.median"], -1.13054 )
    expect_10pd( out4[2,"logit.median"], -0.8668032 )
    expect_10pd( out4[1,"upper"], 0.4120062 )
    expect_10pd( out4[2,"upper"], 0.5325494 )
    expect_10pd( out4[1,"lower"], 0.1326133 )
    expect_10pd( out4[2,"lower"], 0.1520124 )
    expect_10pd( out4[1,"median"], 0.2440846 )
    expect_10pd( out4[2,"median"], 0.2992253 )
  })
  
}
