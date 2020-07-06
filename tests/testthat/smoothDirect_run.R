rm(list=ls())

###########################################################
##
##  Test for smoothDirect
##
###########################################################

library(SUMMER)
data(DemoData)
years <- levels(DemoData[[1]]$time)
years.all <- c(years, "15-19")
data_multi <- getDirectList(births = DemoData, years = years,
regionVar = "region",  timeVar = "time", clusterVar = "~clustid+id",
ageVar = "age", weightsVar = "weights", geo.recode = NULL)
data <- aggregateSurvey(data_multi)

###########################################################

##
## national model, rw2
##  
fit <- smoothDirect(data = data, Amat = NULL, time.model = "rw2",
		year_label = years.all, year_range = c(1985, 2019), 
		is.yearly=TRUE, m = 5)
out1 <- getSmoothed(fit)


##
## subnational model, rw2
##  
fit <- smoothDirect(data = data, Amat = NULL, time.model = "rw2",
		year_label = years.all, year_range = c(1985, 2019), 
		is.yearly=TRUE, m = 5)
out2 <- getSmoothed(fit)


##
## subnational model, period, rw2 and ar1 interaction
##  
fit <- smoothDirect(data = data, Amat = DemoMap$Amat, time.model = "rw2", st.time.model = "ar1",
		year_label = years.all, year_range = c(1985, 2019), 
		is.yearly=FALSE, type.st = 4)
out3 <- getSmoothed(fit, Amat = DemoMap$Amat)

save(data, years.all, out1, out2, out3, file="smoothDirect_regress.RData")
