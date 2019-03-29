#' Function to get Births file from DHS .dta files.
#' 
#'
#' @param filepath file path of raw .dta file from DHS. Only used when data frame is not provided in the function call.
#' @param data data frame of a DHS survey
#' @param surveyyear year of survey
#' @param variables vector of variables to be used in obtaining the person-month files. The variables correspond the the DHS recode manual VI. For early DHS data, the variable names may need to be changed.
#' @param strata vector of variable names used for strata. If a single variable is specified, then that variable will be used as strata indicator If multiple variables are specified, the interaction of these variables will be used as strata indicator. 
#' @param dob variable name for the date of birth.
#' @param alive variable name for the indicator of whether child was alive or dead at the time of interview.
#' @param age variable name for the age at death of the child in completed months.
#' @param date.interview variable name for the date of interview.
#' @param month.cut The cutoff of each bins of age group in the unit of months. Default values are 1, 12, 24, 36, 48, and 60, representing the age groups (0, 1), [1, 12), [12, 24), ..., [48, 60).
#' @param year.cut The cutoff of each bins of time periods, including both boundaries. Default values are 1980, 1985, ..., 2020, representing the time periods 80-84, 85-89, ..., 15-19.
#' 
#' @return This function returns a new data frame where each row indicate a person-month, with the additional variables specified in the function argument.
#' @examples 
#' \dontrun{
#' my_fp <- "/myExampleFilepath/surveyData.DTA"
#' DemoData <- getBirths(filepath = my_fp, surveyyear = 2015) 
#' }
#' 
#' @export
getBirths <- function(filepath = NULL, data = NULL, surveyyear, variables = c("caseid", "v001", "v002", "v004", "v005", "v021", "v022", "v023", "v024", "v025", "v139", "bidx"), strata=c("v024", "v025"), dob = "b3", alive = "b5", age = "b7", date.interview= "v008", month.cut = c(1,12,24,36,48,60), year.cut=seq(1980, 2020, by=5)) {
  if(is.null(data)){
      dat <- suppressWarnings(readstata13::read.dta13(filepath, generate.factors = TRUE))    
  }else{
      dat <- data.frame(data)
  }
  
  surveyyear <- surveyyear - 1900
  year.cut <- year.cut - 1900
  variables <- union(variables, strata)
  datnew <- dat[, variables] 
  
  
  datnew$dob <- dat[, dob]
  datnew$survey_year <- surveyyear
  datnew$obsStart <- dat[, dob]
  datnew$dod <- dat[, dob] + dat[, age]
  datnew$obsStop <- dat[, date.interview]
  datnew$obsStop[dat[, alive] == "no"] <- datnew$dod[dat[, alive] == "no"]
  datnew$died <- (dat[, alive] == "no")
  
  datnew$obsStop[datnew$obsStart == datnew$obsStop] <- datnew$obsStop[datnew$obsStart == datnew$obsStop] + 0.01
  
  datnew$id.new <- 1:nrow(datnew)
  
  # Surv(time = obsStart, time2=obsStop, event = died, origin=dob)~dob+survey_year+died+id + caseid + v001 + v002 + v004 + v005 + v021 + v022 + v023 + v024 + v025 + v139 + bidx
  formula <- as.formula(paste(c("Surv(time = obsStart, time2=obsStop, event = died, origin=dob)~dob+survey_year+died+id.new", union(variables, strata)),  collapse = "+"))

  Surv <- survival::Surv
  test <- survival::survSplit(formula,
                              data = datnew, cut = c(0:max(month.cut)), 
                              start = "agemonth", end = "tstop", event = "died")

  test$obsStart <- test$dob + test$agemonth
  test$obsStop <- test$dob + test$tstop
  test$obsmonth <- test$dob + test$agemonth
  test$year <- floor((test$obsmonth-1)/12)
  
  test <- test[test$agemonth<max(month.cut), ]
  test <- test[test$year>=year.cut[1], ]
  test <- test[test$year<test$survey_year, ]
  
  test$tstop <- NULL
  
  # test$ageGrp6 <- floor(test$agemonth/6)*6
  if(month.cut[1] == 0) month.cut <- month.cut[-1]
  bins <- rep(NA, length(month.cut))
  bins[1] <- paste("0", month.cut[1]-1, sep = "-")
  if(length(month.cut) > 1){
    for(i in 1:(length(month.cut)-1)){
      bins[i+1] <- paste(month.cut[i], month.cut[i+1]-1, sep = "-")
    }    
  }
  if(bins[1] == "0-0") bins[1] <- "0"

  test$age <- bins[1]
  for(i in 1:(length(month.cut)-1)){
    test$age[test$agemonth >= month.cut[i] & test$agemonth < month.cut[i+1]] <- bins[i+1]
  }
  # test$ageGrp12 <- floor(test$agemonth/12)*12
  # test$age <- test$ageGrp12
  # test$age[test$agemonth == 0] <- 0
  # test$age[test$agemonth >= 1 & test$agemonth < 12] <- 1
  # test$ageGrp6 <- factor(test$ageGrp6)
  # levels(test$ageGrp6) <- c("0-5", "6-11", "12-17", "18-23", "24-29", "30-35", "36-41", "42-47", "48-53", "54-59")
  # test$ageGrp12 <- factor(test$ageGrp12)
  # levels(test$ageGrp12) <- c("0-11", "12-23", "24-35", "36-47", "48-59")
  
  test$age <- factor(test$age)
  levels(test$age) <- bins #c("0","1-11","12-23","24-35","36-47","48-59")
  

  year.cut2 <- year.cut
  year.cut2[year.cut2 >= 100] <- as.character(year.cut2[year.cut2 >= 100] - 100)
  year.cut2[as.numeric(year.cut2) < 10] <- paste0("0", year.cut2[as.numeric(year.cut2) < 10])
  year.cut3 <- year.cut - 1
  year.cut3[year.cut3 >= 100] <- as.character(year.cut3[year.cut3 >= 100] - 100)
  year.cut3[as.numeric(year.cut3) < 10] <- paste0("0", year.cut3[as.numeric(year.cut3) < 10])

  test$time <- paste(year.cut[1], year.cut[2] - 1, sep = "-")
  year.bin <- paste(year.cut[1], year.cut[2] - 1, sep = "-")
  for(i in 2:(length(year.cut)-1)){
    test$time[test$year >= year.cut[i] & test$year < year.cut[i+1]] <- paste(year.cut2[i], year.cut3[i+1], sep = "-")
    year.bin <- c(year.bin, paste(year.cut2[i], year.cut3[i+1], sep = "-"))
  }
  # test$time[test$year > 89] <- "90-94"
  # test$time[test$year > 94] <- "95-99"
  # test$time[test$year > 99] <- "00-04"
  # test$time[test$year > 104] <- "05-09"
  # test$time[test$year > 109] <- "10-14"
  test$time <- factor(test$time, levels = year.bin)
  
  if(length(strata) == 0){
    test$strata <- NA
  }else if(length(strata) == 1){
    test$strata <- test[, strata]
  }else{
    test$strata <- do.call(paste, c(test[strata], sep="."))
  }
  test$survey_year <- test$survey_year + 1900
  return(test)
}