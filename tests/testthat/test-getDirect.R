##
## These tests cover all @examples of getDirect() and getDirectList()
##

test_that("getDirect works", {
    
    skip_on_cran()

    data(DemoData)
    years <- c("85-89", "90-94", "95-99", "00-04", "05-09", "10-14")
    mean <- getDirect(births = DemoData[[1]],  years = years, 
    regionVar = "region", timeVar = "time", clusterVar = "~clustid+id", 
    ageVar = "age", weightsVar = "weights", geo.recode = NULL)

    # No dropping entries
    expect_equal(dim(mean)[1], 30)
    # No NA in mean
    expect_equal(sum(is.na(mean$mean)), 0)
    # No NA in logit precision
    expect_equal(sum(is.na(mean$logit.prec)), 0)
    # First region is "All"
    expect_equal(mean[1,1], "All")
})

test_that("getDirectList works", {
    
    skip_on_cran()

    data(DemoData)
    years <- c("85-89", "90-94", "95-99", "00-04", "05-09", "10-14")
    mean <- getDirectList(births = DemoData,  years = years, 
    regionVar = "region", timeVar = "time", clusterVar = "~clustid+id", 
    ageVar = "age", weightsVar = "weights", geo.recode = NULL)

    # No dropping entries
    expect_equal(dim(mean)[1], 150)
    # Two NA in mean
    expect_equal(sum(is.na(mean$mean)), 2)
    # Two NA in logit precision
    expect_equal(sum(is.na(mean$logit.prec)), 2) 
    # No NA in survey index
    expect_equal(sum(is.na(mean$survey)), 0)
    # First region is "All"
    expect_equal(mean[1,1], "All")
})

