context("test-smoothDirect")

expect_5pd_prob <- function(x, y) {
  raw_difference <- abs(x-y)  
  return( expect_lte(max(raw_difference), 0.05) )
}

load("smoothDirect_regress.RData")

  
## run tests
test_that("Smooth Direct comparison: RW2 National", {
  
	
	## national model, rw2
	##  
	fit <- smoothDirect(data = data, Amat = NULL, time.model = "rw2",
			year_label = years.all, year_range = c(1985, 2019), 
			is.yearly=TRUE, m = 5)
	newout1 <- getSmoothed(fit)



  # regression test to see if the results match expectation
  expect_equal(newout1$region, out1$region)
  expect_equal(newout1$years, out1$years)
  expect_5pd_prob(newout1$upper, out1$upper)
  expect_5pd_prob(newout1$lower, out1$lower)
  expect_5pd_prob(newout1$median, out1$median)
  
})



test_that("Smooth Direct comparison: RW2 Subnational", {

	##
	## subnational model, rw2
	##  
	fit <- smoothDirect(data = data, Amat = NULL, time.model = "rw2",
			year_label = years.all, year_range = c(1985, 2019), 
			is.yearly=TRUE, m = 5)
	newout2 <- getSmoothed(fit)

	expect_equal(newout2$region, out2$region)
	expect_equal(newout2$years, out2$years)
	expect_5pd_prob(newout2$upper, out2$upper)
	expect_5pd_prob(newout2$lower, out2$lower)
	expect_5pd_prob(newout2$median, out2$median)
 
 }) 


test_that("Smooth Direct comparison: RW2 And AR1 Interaction", {

	##
	## subnational model, period, rw2 and ar1 interaction
	##  
	fit <- smoothDirect(data = data, Amat = DemoMap$Amat, time.model = "rw2", st.time.model = "ar1",
			year_label = years.all, year_range = c(1985, 2019), 
			is.yearly=FALSE, type.st = 4)
	newout3 <- getSmoothed(fit, Amat = DemoMap$Amat)


	expect_equal(newout3$region, out3$region)
	expect_equal(newout3$years, out3$years)
	expect_5pd_prob(newout3$upper, out3$upper)
	expect_5pd_prob(newout3$lower, out3$lower)
	expect_5pd_prob(newout3$median, out3$median)
 
 }) 
