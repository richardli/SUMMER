##
## These tests cover all @examples of smoothDirect()
##

test_that("smoothDirect works for national model", {

   skip_on_cran()
   library(INLA)
   # make devtools::check() happy with single process
   inla.setOption( num.threads = 1 )

   data(DemoData)
   years <- levels(DemoData[[1]]$time)
   # obtain direct estimates
   data_multi <- getDirectList(births = DemoData, years = years,
   regionVar = "region",  timeVar = "time", clusterVar = "~clustid+id",
   ageVar = "age", weightsVar = "weights", geo.recode = NULL)
   data <- aggregateSurvey(data_multi)
   
   # Check aggregateSurvey works
   expect_equal(dim(data)[1], 30)


   #  national model
   years.all <- c(years, "15-19")
   fit1 <- smoothDirect(data = data, Amat = NULL, 
   year.label = years.all, year.range = c(1985, 2019), 
   time.model = 'rw2', m = 5, control.compute = list(config =TRUE))

   # function is quite mature, check for finishing only
   expect_equal(class(fit1), "SUMMERmodel")

   # check for smoothed output
   out1 <- getSmoothed(fit1)
   expect_equal(dim(out1)[1], 42)

})

test_that("smoothDirect works for subnational model", {

   skip_on_cran()
   library(INLA)
   # make devtools::check() happy with single process
   inla.setOption( num.threads = 1 )

   data(DemoData)
   years <- levels(DemoData[[1]]$time)
   years.all <- c(years, "15-19")

   # obtain direct estimates
   data_multi <- getDirectList(births = DemoData, years = years,
   regionVar = "region",  timeVar = "time", clusterVar = "~clustid+id",
   ageVar = "age", weightsVar = "weights", geo.recode = NULL)
   data <- aggregateSurvey(data_multi)


   fit1 <- smoothDirect(data = data, Amat = NULL, 
   year.label = years.all, year.range = c(1985, 2019), 
   time.model = 'rw2', m = 5, control.compute = list(config =TRUE))
    # check for smoothed output
   out1 <- getSmoothed(fit1)


   #  subnational model
   fit2 <- smoothDirect(data = data, Amat = DemoMap$Amat, 
   year.label = years.all, year.range = c(1985, 2019), 
   time.model = 'rw2', m = 5, type.st = 4)

  # function is quite mature, check for finishing only
   expect_equal(class(fit2), "SUMMERmodel")
 
  # check for smoothed output
   out2 <- getSmoothed(fit2)
   expect_equal(dim(out2)[1], 168)
 
   # check for consistency in columnnames
   expect_equal(colnames(out1), colnames(out2))

 })

test_that("smoothDirect works for subnational space-only model", {

   skip_on_cran()
   library(INLA)
   # make devtools::check() happy with single process
   inla.setOption( num.threads = 1 )

   data(DemoData)
   years <- levels(DemoData[[1]]$time)
   years.all <- c(years, "15-19")
   # obtain direct estimates
   data_multi <- getDirectList(births = DemoData, years = years,
   regionVar = "region",  timeVar = "time", clusterVar = "~clustid+id",
   ageVar = "age", weightsVar = "weights", geo.recode = NULL)
   data <- aggregateSurvey(data_multi)

   fit1 <- smoothDirect(data = data, Amat = NULL, 
   year.label = years.all, year.range = c(1985, 2019), 
   time.model = 'rw2', m = 5, control.compute = list(config =TRUE))
    # check for smoothed output
   out1 <- getSmoothed(fit1)


   #  subnational space-only model for one period
   fit3 <- smoothDirect(data = subset(data, years == "10-14"), 
            time.model = NULL, Amat = DemoMap$Amat)

   # function is quite mature, check for finishing only
   expect_equal(class(fit3), "SUMMERmodel")
 
   # check for smoothed output
   out3 <- getSmoothed(fit3)
   expect_equal(dim(out3)[1], 4)

   # check for consistency in columnnames
   expect_equal(colnames(out1), colnames(out3))

   # check plot does not give error, not correctness
    g <- plot(out3, plot.CI=TRUE)  
    expect_equal(class(g), c("gg", "ggplot"))

})
