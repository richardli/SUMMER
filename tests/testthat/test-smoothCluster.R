##
## These tests cover all @examples of smoothCluster()
##

test_that("getCounts and smoothCluster works without covariates", {

 skip_on_cran()

 library(dplyr)
 library(INLA)
 library(sn)
 # make devtools::check() happy with single process
 inla.setOption( num.threads = 1 )

 data(DemoData)
 # Create dataset of counts
 counts.all <- NULL
 for(i in 1:length(DemoData)){
   counts <- getCounts(DemoData[[i]][, c("clustid", "time", "age", "died",
                                         "region", "strata")],
             variables = 'died', by = c("age", "clustid", "region", 
                                          "time", "strata"))
   counts <- counts %>% mutate(cluster = clustid, years = time, Y=died)
   counts$strata <- gsub(".*\\.","",counts$strata)
   counts$survey <- names(DemoData)[i] 
   counts.all <- rbind(counts.all, counts)
 }
 # check for counts dimension
 expect_equal(dim(counts.all), c(16908, 11))
 
 # fit cluster-level model on the periods
 periods <- levels(DemoData[[1]]$time)
 fit <- smoothCluster(data = counts.all, 
      Amat = DemoMap$Amat, 
      time.model = "rw2", 
      st.time.model = "rw1",
      strata.time.effect =  TRUE, 
      survey.effect = TRUE,
      family = "betabinomial",
      year.label = c(periods, "15-19"))

 # function is quite mature, check for finishing only
 expect_equal(class(fit), "SUMMERmodel")
 expect_equal(fit$family, "betabinomial")
 # checking fixed effect dimension: 6 age * 2 strata
 expect_equal(dim(fit$fit$summary.fixed)[1], 6*2)
 # checking time trends: 3 age groups (1, 2, next 4) * 2 strata * 7 time periods
 expect_equal(dim(fit$fit$summary.random$time.struct)[1], 3*2*7)
 expect_equal(fit$age.rw.group, c(1, 2, rep(3, 4), 4, 5, rep(6, 4)))

# check shared temporal trends across strata 
fit2 <- smoothCluster(data = counts.all, 
      Amat = DemoMap$Amat, 
      time.model = "rw2", 
      st.time.model = "rw1",
      strata.time.effect =  FALSE, 
      survey.effect = TRUE,
      family = "betabinomial",
      year.label = c(periods, "15-19"))

 # checking fixed effect dimension: 6 age * 2 strata + overall ur effect
 expect_equal(dim(fit2$fit$summary.fixed)[1], 6*2+1)
# checking time trends: 3 age groups (1, 2, next 4) * 7 time periods
 expect_equal(dim(fit2$fit$summary.random$time.struct)[1], 3*7)
 expect_equal(fit2$age.rw.group, c(1, 2, rep(3, 4)))

 
 est <- getSmoothed(fit, nsim = 100)
 # check dimensions
 expect_equal(dim(est$stratified), c(56, 12)) 
 # check strata variable
 expect_equal("strata" %in% colnames(est$stratified), TRUE) 
 
 # check dimensions
 expect_equal(dim(est$overall), c(28, 14))
 # check frame variable
 expect_equal("frame" %in% colnames(est$overall), TRUE) 
 # check overall is 0 since no weight returned
 expect_equal(est$overall$mean[1], 0)
 
 # check plot does not give error, not correctness
 g <- plot(est$stratified, plot.CI=TRUE) + ggplot2::facet_wrap(~strata) 
 expect_equal(class(g), c("gg", "ggplot"))

})



test_that("getCounts and smoothCluster works with covariates", {


  skip_on_cran()

 library(dplyr)
 data(DemoData)
 # Create dataset of counts
 counts.all <- NULL
 for(i in 1:length(DemoData)){
   counts <- getCounts(DemoData[[i]][, c("clustid", "time", "age", "died",
                                         "region", "strata")],
             variables = 'died', by = c("age", "clustid", "region", 
                                          "time", "strata"))
   counts <- counts %>% mutate(cluster = clustid, years = time, Y=died)
   counts$strata <- gsub(".*\\.","",counts$strata)
   counts$survey <- names(DemoData)[i] 
   counts.all <- rbind(counts.all, counts)
 }
 # fit cluster-level space-time model with covariate
 # notice without projected covariates, we use periods up to 10-14 only
 # construct a random covariate matrix for illustration
 periods <- levels(DemoData[[1]]$time)
 X <- expand.grid(years = periods, 
        region = unique(counts.all$region))
 X$X1 <- rnorm(dim(X)[1])
 X$X2 <- rnorm(dim(X)[1])
 fit.covariate <- smoothCluster(data = counts.all, 
    X = X,
      Amat = DemoMap$Amat, 
      time.model = "rw2", 
      st.time.model = "rw1",
      strata.time.effect =  TRUE, 
      survey.effect = TRUE,
      family = "betabinomial",
      year.label = c(periods))
 
 # function is quite mature, check for finishing only
 expect_equal(class(fit.covariate), "SUMMERmodel")
 expect_equal(fit.covariate$family, "betabinomial")
 # checking fixed effect dimension: 6 age * 2 strata + 2 covariates
 expect_equal(dim(fit.covariate$fit$summary.fixed)[1], 6*2+2)

 est <- getSmoothed(fit.covariate, nsim = 100)
 # check dimensions
 expect_equal(dim(est$stratified), c(48, 12)) 
 # check strata variable
 expect_equal("strata" %in% colnames(est$stratified), TRUE) 
 
 # check dimensions
 expect_equal(dim(est$overall), c(24, 14))
 # check frame variable
 expect_equal("frame" %in% colnames(est$overall), TRUE) 
 # check overall is 0 since no weight returned
 expect_equal(est$overall$mean[1], 0)


 
})



test_that("getCounts and smoothCluster works with space-only model", {


  skip_on_cran()

 library(dplyr)
 data(DemoData)
 # Create dataset of counts
 counts.all <- NULL
 for(i in 1:length(DemoData)){
   counts <- getCounts(DemoData[[i]][, c("clustid", "time", "age", "died",
                                         "region", "strata")],
             variables = 'died', by = c("age", "clustid", "region", 
                                          "time", "strata"))
   counts <- counts %>% mutate(cluster = clustid, years = time, Y=died)
   counts$strata <- gsub(".*\\.","",counts$strata)
   counts$survey <- names(DemoData)[i] 
   counts.all <- rbind(counts.all, counts)
 }

 # fit cluster-level model for one time point only
 # i.e., space-only model
 fit.sp <- smoothCluster(data = subset(counts.all, time == "10-14"), 
      Amat = DemoMap$Amat, 
      time.model = NULL, 
      survey.effect = TRUE,
      family = "betabinomial")

 # function is quite mature, check for finishing only
 expect_equal(class(fit.sp), "SUMMERmodel")
 expect_equal(fit.sp$family, "betabinomial")
 # checking fixed effect dimension: 6 age * 2 strata + 1 overall u/r
 expect_equal(dim(fit.sp$fit$summary.fixed)[1], 6*2+1)

 est <- getSmoothed(fit.sp, nsim = 100)
 # check dimensions
 expect_equal(dim(est$stratified), c(4*2, 12)) 
 # check strata variable
 expect_equal("strata" %in% colnames(est$stratified), TRUE) 
 
 # check dimensions
 expect_equal(dim(est$overall), c(4, 14))
 # check frame variable
 expect_equal("frame" %in% colnames(est$overall), TRUE) 
 # check overall is 0 since no weight returned
 expect_equal(est$overall$mean[1], 0)


 # fit cluster-level model for one time point and covariate
 # construct a random covariate matrix for illustration
 X <- data.frame(region = unique(counts.all$region),
       X1 = c(1, 2, 2, 1), 
       X2 = c(1, 1, 1, 2))
 fit.sp.covariate <- smoothCluster(data = subset(counts.all, time == "10-14"), 
      X = X, 
      Amat = DemoMap$Amat, 
      time.model = NULL, 
      survey.effect = TRUE,
      family = "betabinomial")
expect_equal(dim(fit.sp.covariate$fit$summary.fixed)[1], 6*2+1+2)

})
