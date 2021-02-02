#' The BRFSS dataset 
#' 
#' The Behavioral Risk Factor Surveillance System (BRFSS) is an annual telephone health survey conducted by the Centers for Disease Control and Prevention (CDC) that tracks health conditions and risk behaviors in the United States and its territories since 1984. This BRFSS dataset contains 16124 observations. The `diab2` variable is the binary indicator of Type II diabetes, `strata` is the strata indicator and `rwt_llcp` is the final design weight. Records with missing HRA code or diabetes status are removed from this dataset. See \url{https://www.cdc.gov/brfss/annual_data/2013/pdf/Weighting_Data.pdf} for more details of the weighting procedure.
#'
#' @format A data.frame of 26 variables. 
#' @references Washington State Department of Health, Center for Health Statistics. Behavioral Risk Factor Surveillance System, supported in part by the Centers for Disease Control and Prevention. Corporative Agreement U58/DP006066-01 (2015).  

#' @docType data
#'
#' @usage data(BRFSS)
"BRFSS"