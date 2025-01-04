test_that("getDiag works: period smoothDirect", {

   skip_on_cran()
   library(INLA)
   # make devtools::check() happy with single process
   inla.setOption( num.threads = 1 )

	data(DemoMap)
	years <- levels(DemoData[[1]]$time)

	# obtain direct estimates
	data <- getDirectList(births = DemoData, 
	years = years,  
	regionVar = "region", timeVar = "time", 
	clusterVar = "~clustid+id", 
	ageVar = "age", weightsVar = "weights", 
	geo.recode = NULL)
	# obtain direct estimates
	data_multi <- getDirectList(births = DemoData, years = years,
	  regionVar = "region",  timeVar = "time", clusterVar = "~clustid+id",
	  ageVar = "age", weightsVar = "weights", geo.recode = NULL)
	data <- aggregateSurvey(data_multi)

	#  national model
	years.all <- c(years, "15-19")
	fit1 <- smoothDirect(data = data, geo = DemoMap$geo, Amat = DemoMap$Amat, 
	  year.label = years.all, year.range = c(1985, 2019), 
	  rw = 2, is.yearly=FALSE, m = 5)
	random.time <- getDiag(fit1, field = "time")
	random.space <- getDiag(fit1, field = "space")
	random.spacetime <- getDiag(fit1, field = "spacetime")

   	expect_equal(nrow(random.time), 7*2)
   	expect_equal(nrow(random.space), 4*2)
   	expect_equal(nrow(random.spacetime), 4*7)
})

test_that("getDiag works: period-to-year smoothDirect", {

   skip_on_cran()
   skip_on_ci()
   library(INLA)
   # make devtools::check() happy with single process
   inla.setOption( num.threads = 1 )

	data(DemoMap)
	years <- levels(DemoData[[1]]$time)

	# obtain direct estimates
	data <- getDirectList(births = DemoData, 
	years = years,  
	regionVar = "region", timeVar = "time", 
	clusterVar = "~clustid+id", 
	ageVar = "age", weightsVar = "weights", 
	geo.recode = NULL)
	# obtain direct estimates
	data_multi <- getDirectList(births = DemoData, years = years,
	  regionVar = "region",  timeVar = "time", clusterVar = "~clustid+id",
	  ageVar = "age", weightsVar = "weights", geo.recode = NULL)
	data <- aggregateSurvey(data_multi)

	#  national model
	years.all <- c(years, "15-19")

   #  subnational model
   fit2 <- smoothDirect(data = data, Amat = DemoMap$Amat, 
	   year.label = years.all, year.range = c(1985, 2019), 
	   time.model = 'rw2', m = 5, type.st = 4)	

	random.time <- getDiag(fit2, field = "time")
	random.space <- getDiag(fit2, field = "space")
	random.spacetime <- getDiag(fit2, field = "spacetime")

	expect_equal(nrow(random.time), (7+35)*2)
	expect_equal(nrow(random.space), 4*2)
	expect_equal(nrow(random.spacetime), 4*(7+35))

})

test_that("getDiag works: smoothCluster", {

	skip_on_cran()
	library(INLA)
	library(dplyr)
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
	

	# fit cluster-level model on the periods
	periods <- levels(DemoData[[1]]$time)
	fit3 <- smoothCluster(data = counts.all, 
	  Amat = DemoMap$Amat, 
	  time.model = "rw2", 
	  st.time.model = "rw1",
	  strata.time.effect =  TRUE, 
	  survey.effect = TRUE,
	  family = "betabinomial",
	  year.label = c(periods, "15-19"))


	random.time <- getDiag(fit3, field = "time")
	random.space <- getDiag(fit3, field = "space")
	random.spacetime <- getDiag(fit3, field = "spacetime")

	expect_equal(nrow(random.time), 7*6*2+7)
	expect_equal(nrow(random.space), 4*2)
	expect_equal(nrow(random.spacetime), 4*7)

})



test_that("getDiag works: smoothCluster with linear trend", {

	skip_on_cran()
	skip_on_ci()
	library(INLA)
	library(dplyr)
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
	

	# fit cluster-level model on the periods
	periods <- levels(DemoData[[1]]$time)
	fit4 <- smoothCluster(data = counts.all, 
	  Amat = DemoMap$Amat, 
	  time.model = "ar1", 
	  st.time.model = "rw2",
	  strata.time.effect =  TRUE, 
	  linear.trend = TRUE,
	  survey.effect = TRUE,
	  family = "betabinomial",
	  year.label = c(periods, "15-19"))


	random.time <- getDiag(fit4, field = "time")
	random.space <- getDiag(fit4, field = "space")
	random.spacetime <- getDiag(fit4, field = "spacetime")

	expect_equal(nrow(random.time), 7*6*2+7)
	expect_equal(nrow(random.space), 4*2)
	expect_equal(nrow(random.spacetime), 4*7)

})


test_that("getDiag works: smoothCluster with random slopes", {

	skip_on_cran()
	library(INLA)
	library(dplyr)
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
	

	# fit cluster-level model on the periods
	periods <- levels(DemoData[[1]]$time)
	fit5 <- smoothCluster(data = counts.all, 
	  Amat = DemoMap$Amat, 
	  time.model = "rw2", 
	  st.time.model = "ar1",
	  pc.st.slope.u = 2, 
	  pc.st.slope.alpha = 0.1,
	  strata.time.effect =  TRUE, 
	  survey.effect = TRUE,
	  family = "betabinomial",
	  year.label = c(periods, "15-19"))


	random.time <- getDiag(fit5, field = "time")
	random.space <- getDiag(fit5, field = "space")
	random.spacetime <- getDiag(fit5, field = "spacetime")

	expect_equal(nrow(random.time), 7*6*2+7)
	expect_equal(nrow(random.space), 4*2)
	expect_equal(nrow(random.spacetime), 4*7)

})




test_that("getCounts and smoothCluster works with space-only model", {


	skip_on_cran()

	library(dplyr)
	library(INLA)
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

	# fit cluster-level model for one time point only
	# i.e., space-only model
	fit.sp <- smoothCluster(data = subset(counts.all, time == "10-14"), 
	  Amat = DemoMap$Amat, 
	  time.model = NULL, 
	  survey.effect = TRUE,
	  family = "betabinomial")


	random.space <- getDiag(fit.sp, field = "space")

	expect_equal(nrow(random.space), 4*2)

})

