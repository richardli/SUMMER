##
## These tests cover all @examples of smoothSurvey()
##

test_that("smoothSurvey: Area-level model works", {

   skip_on_cran()
   library(INLA)
   # make devtools::check() happy with single process
   inla.setOption( num.threads = 1 )

   data(DemoData2)
   data(DemoMap2)
   fit0 <- smoothSurvey(data=DemoData2,  
    Amat=DemoMap2$Amat, response.type="binary", 
    responseVar="tobacco.use", strataVar="strata", 
    weightVar="weights", regionVar="region", 
    clusterVar = "~clustid+id", CI = 0.95)
   expect_equal(class(fit0), "SUMMERmodel.svy")
   expect_equal(dim(fit0$direct), c(8, 6))
   expect_equal(dim(fit0$smooth), c(8, 11))

   # if only direct estimates without smoothing is of interest
   fit0.dir <- smoothSurvey(data=DemoData2,  
      Amat=DemoMap2$Amat, response.type="binary", 
      responseVar="tobacco.use", strataVar="strata", 
      weightVar="weights", regionVar="region", 
      clusterVar = "~clustid+id", CI = 0.95, smooth = FALSE)

   expect_equal(class(fit0.dir), "SUMMERmodel.svy")
   expect_equal(dim(fit0.dir$direct), c(8, 6))
   expect_equal(sum(is.na(fit0.dir$smooth$mean)), 8)


   # posterior draws can be returned with save.draws = TRUE
   fit0.draws <- smoothSurvey(data=DemoData2,  
      Amat=DemoMap2$Amat, response.type="binary", 
      responseVar="tobacco.use", strataVar="strata", 
      weightVar="weights", regionVar="region", 
      clusterVar = "~clustid+id", CI = 0.95, 
      save.draws = TRUE, nsim = 100)
   # notice the posterior draws are on the latent scale
   expect_equal(dim(fit0.draws$draws.est), c(8, 102))
   # check no NAs except in the time column
   expect_equal(sum(is.na(fit0.draws$draws.est[, -2])), 0)

   # Example with region-level covariates
   Xmat <- aggregate(age~region, data = DemoData2, 
                FUN = function(x) mean(x))
   fit1 <- smoothSurvey(data=DemoData2,  
      Amat=DemoMap2$Amat, response.type="binary", 
      X = Xmat,
      responseVar="tobacco.use", strataVar="strata", 
      weightVar="weights", regionVar="region", 
      clusterVar = "~clustid+id", CI = 0.95)
   expect_equal(class(fit1), "SUMMERmodel.svy")
   expect_equal(dim(fit1$direct), c(8, 6+1))
   expect_equal(dim(fit1$smooth), c(8, 11))
   expect_equal(sum(is.na(fit1$smooth$mean)), 0)




   # Example with using only direct estimates as input instead of the full data
   direct <- fit0$direct[, c("region", "direct.est", "direct.var")]
   fit2 <- smoothSurvey(data=NULL, direct.est = direct, 
                  Amat=DemoMap2$Amat, regionVar="region",
                  responseVar="direct.est", direct.est.var = "direct.var", 
                  response.type = "binary")

   expect_equal(class(fit2), "SUMMERmodel.svy")
   # check they are the same as using raw input
   expect_equal(fit2$smooth$mean, fit0$smooth$mean, tolerance = 0.01)
   expect_equal(fit2$smooth$var, fit0$smooth$var, tolerance = 0.01)

 

   # Example with using only direct estimates as input, 
   #   and after transformation into a Gaussian smoothing model
   # Notice: the output are on the same scale as the input 
   #   and in this case, the logit estimates.    
   direct.logit <- fit0$direct[, c("region", "direct.logit.est", "direct.logit.var")]
   fit3 <- smoothSurvey(data=NULL, direct.est = direct.logit, 
             Amat=DemoMap2$Amat, regionVar="region",
             responseVar="direct.logit.est", direct.est.var = "direct.logit.var",
             response.type = "gaussian")
   expect_equal(fit3$smooth$mean, fit0$smooth$logit.mean, tolerance = 0.01)
   expect_equal(fit3$smooth$var, fit0$smooth$logit.var, tolerance = 0.01)

   # Example with non-spatial smoothing using IID random effects
   fit4 <- smoothSurvey(data=DemoData2, response.type="binary", 
     responseVar="tobacco.use", strataVar="strata", 
     weightVar="weights", regionVar="region", 
     clusterVar = "~clustid+id", CI = 0.95)
   expect_equal(class(fit4), "SUMMERmodel.svy")
   expect_equal(dim(fit4$direct), c(8, 6))
   expect_equal(dim(fit4$smooth), c(8, 11))

})



test_that("smoothSurvey: Area-level model works with missing data", {

   skip_on_cran()
   library(INLA)
   # make devtools::check() happy with single process
   inla.setOption( num.threads = 1 )

   # Example with missing regions in the raw input
   # When no Amat is given, the missing area is completely gone
   DemoData2.sub <- subset(DemoData2, region != "central")
   fit.without.central <- smoothSurvey(data=DemoData2.sub,  
                       Amat=NULL, response.type="binary", 
                       responseVar="tobacco.use", strataVar="strata", 
                       weightVar="weights", regionVar="region", 
                       clusterVar = "~clustid+id", CI = 0.95)
   expect_equal(class(fit.without.central), "SUMMERmodel.svy")
   expect_equal(dim(fit.without.central$direct), c(7, 6))
   expect_equal(dim(fit.without.central$smooth), c(7, 11))


   # Example with missing regions in the raw input
   # When no Amat is given, but region.list is provided
   # the missing area is filled in
   fit.with.central <- smoothSurvey(data=DemoData2.sub,  
                       Amat=NULL, region.list = unique(DemoData2$region),
                       response.type="binary", 
                       responseVar="tobacco.use", strataVar="strata", 
                       weightVar="weights", regionVar="region", 
                       clusterVar = "~clustid+id", CI = 0.95)
   expect_equal(class(fit.with.central), "SUMMERmodel.svy")
   expect_equal(dim(fit.with.central$direct), c(8, 6))
   expect_equal(is.na(fit.with.central$direct$direct.est[fit.with.central$direct$region == "central"]), TRUE)
   expect_equal(dim(fit.with.central$smooth), c(8, 11))

   # Example with missing regions in the raw input
   # When Amat is given,
   fit.with.central.spatial <- smoothSurvey(data=DemoData2.sub,  
                       Amat=DemoMap2$Amat,
                       response.type="binary", 
                       responseVar="tobacco.use", strataVar="strata", 
                       weightVar="weights", regionVar="region", 
                       clusterVar = "~clustid+id", CI = 0.95)
   expect_equal(class(fit.with.central.spatial), "SUMMERmodel.svy")
   expect_equal(dim(fit.with.central.spatial$direct), c(8, 6))
   expect_equal(is.na(fit.with.central.spatial$direct$direct.est[fit.with.central.spatial$direct$region == "central"]), TRUE)
   expect_equal(dim(fit.with.central.spatial$smooth), c(8, 11))

})



test_that("smoothSurvey: Area-level model works with customized formula", {

   skip_on_cran()
   library(INLA)
   # make devtools::check() happy with single process
   inla.setOption( num.threads = 1 )


   # Using the formula argument, further customizations can be added to the 
   #  model fitted. For example, we can fit the Fay-Harriot model with 
   #  IID effect instead of the BYM2 random effect as follows.
   #  The "region.struct" and "hyperpc1" are picked to match internal object 
   #  names. Other object names can be inspected from the source of smoothSurvey.
   fit5 <- smoothSurvey(data=DemoData2,  
        Amat=DemoMap2$Amat, response.type="binary", 
        formula = "f(region.struct, model = 'iid', hyper = hyperpc1)",
        pc.u = 1, pc.alpha = 0.01,
        responseVar="tobacco.use", strataVar="strata", 
        weightVar="weights", regionVar="region", 
        clusterVar = "~clustid+id", CI = 0.95)
   # Check BYM2 is replaced by IID
   expect_equal(class(fit5), "SUMMERmodel.svy")
   expect_equal(fit5$fit$model.random, "IID model")



})

test_that("smoothSurvey: Cluster-level model works", {

   skip_on_cran()
   library(INLA)
   # make devtools::check() happy with single process
   inla.setOption( num.threads = 1 )


   # For unit-level models, we need to create stratification variable within regions
   data <- DemoData2
   data$urbanicity <- "rural"
   data$urbanicity[grep("urban", data$strata)] <- "urban"
   
   Xmat <- aggregate(age~region, data = DemoData2, 
                FUN = function(x) mean(x))
   # Beta-binomial likelihood is used in this model
   fit6 <- smoothSurvey(data=data, 
         Amat=DemoMap2$Amat, response.type="binary", 
         X = Xmat, is.unit.level = TRUE,
         responseVar="tobacco.use", strataVar.within = "urbanicity", 
         regionVar="region", clusterVar = "clustid", CI = 0.95)

   expect_equal(class(fit6), "SUMMERmodel.svy")
   expect_equal(dim(fit6$smooth), c(8*2, 13))


   # We may use aggregated PSU-level counts as input as well
   #    in the case of modeling a binary outcome 
   data.agg <- aggregate(tobacco.use~region + urbanicity + clustid, 
                    data = data, FUN = sum)
   data.agg.total <- aggregate(tobacco.use~region + urbanicity + clustid, 
                    data = data, FUN = length)
   colnames(data.agg.total)[4] <- "total"
   data.agg <- merge(data.agg, data.agg.total)

   fit7 <- smoothSurvey(data=data.agg, 
         Amat=DemoMap2$Amat, response.type="binary", 
         X = Xmat, is.unit.level = TRUE, is.agg = TRUE,
         responseVar = "tobacco.use", strataVar.within = "urbanicity", 
         totalVar = "total", regionVar="region", clusterVar = "clustid", CI = 0.95)

   expect_equal(fit6$smooth, fit7$smooth, tolerance = 0.01)




})

test_that("smoothSurvey: Cluster-level model works with Gaussian likelihood", {

   skip_on_cran()
   library(INLA)
   # make devtools::check() happy with single process
   inla.setOption( num.threads = 1 )


   # For unit-level models, we need to create stratification variable within regions
   data <- DemoData2
   data$urbanicity <- "rural"
   data$urbanicity[grep("urban", data$strata)] <- "urban"


   fit11 <- smoothSurvey(data = data, 
         Amat = DemoMap2$Amat, response.type = "gaussian", 
         is.unit.level = TRUE, responseVar="age", strataVar.within = "urbanicity",
         regionVar = "region", clusterVar = NULL, CI = 0.95)  
   expect_equal(class(fit11), "SUMMERmodel.svy")
   expect_equal(dim(fit11$smooth), c(8*2, 13))


   # Notice the usual output is now stratified within each region
   # The aggregated estimates require strata proportions for each region
   # For illustration, we set strata population proportions below
   prop <- data.frame(region = unique(data$region), 
                          urban = 0.3, 
                          rural = 0.7)
   fit12 <- smoothSurvey(data=data, 
               Amat=DemoMap2$Amat, response.type="gaussian", 
               is.unit.level = TRUE, responseVar="age", strataVar.within = "urbanicity",
               regionVar="region", clusterVar = NULL, CI = 0.95,
               weight.strata = prop)  
   expect_equal(class(fit12), "SUMMERmodel.svy")
   expect_equal(dim(fit12$smooth), c(8*2, 13))
   expect_equal(dim(fit12$smooth.overall), c(8, 13 - 1))


   # aggregated outcome
   head(fit12$smooth.overall)

   # Compare aggregated outcome with direct aggregating posterior means. 
   # There could be some differences due to the approximating nature, so only
   # check that the correlation is greater than 0.9.
   est.urb <- subset(fit11$smooth, strata == "urban")
   est.rural <- subset(fit11$smooth, strata == "rural")
   est.mean.post <- est.urb$mean * 0.3 + est.rural$mean * 0.7
   cor <- cor(fit12$smooth.overall$mean, est.mean.post)
   expect_equal(cor, 1, tolerance = 0.1)


   ##
   ## 6. Unit-level model with continuous response and unit-level covariate 
   ## 

   # For area-level prediction, area-level covariate mean needs to be  
   #   specified in X argument. And unit-level covariate names are specified
   #   in X.unit argument.

   set.seed(1)
   sim <- data.frame(region = rep(c(1, 2, 3, 4), 1000),
                 X1 = rnorm(4000), X2 = rnorm(4000))
   Xmean <- aggregate(.~region, data = sim, FUN = sum)
   sim$Y <- rnorm(4000, mean = sim$X1 + 0.3 * sim$X2 + sim$region)
   samp <- sim[sample(1:4000, 20), ]
   fit.sim <- smoothSurvey(data=samp , 
                X.unit = c("X1", "X2"),
                X = Xmean, Amat=NULL, response.type="gaussian", 
                is.unit.level = TRUE, responseVar="Y", regionVar = "region",  
                pc.u = 1, pc.alpha = 0.01, CI = 0.95) 
   expect_equal(class(fit.sim), "SUMMERmodel.svy")
   expect_equal(dim(fit.sim$smooth), c(4, 13))

})
