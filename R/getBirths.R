#' Function to get Births file from DHS .dta files.
#' 
#'
#' @param filepath file path of raw .dta file from DHS
#' @param surveyyear year of survey
#' @param variables vector of variables to be used in obtaining the person-month files. The variables correspond the the DHS recode manual VI. For early DHS data, the variable names may need to be changed.
#' @param strata vector of variable names used for strata. If a single variable is specified, then that variable will be used as strata indicator If multiple variables are specified, the interaction of these variables will be used as strata indicator. 
#' 
#' @return A list of birth-month data
#' 
#' @examples 
#' \dontrun{
#' my_fp <- "/myExampleFilepath/surveyData.DTA"
#' DemoData <- getBirths(filepath = my_fp, surveyyear = 2015) 
#' }
#' 
#' @export
getBirths <- function(filepath,surveyyear, variables = c("caseid", "v001", "v002", "v004", "v005", "v021", "v022", "v023", "v024", "v025", "v139", "bidx"), strata=c("v024", "v025")) {
  dat <- suppressWarnings(readstata13::read.dta13(filepath, generate.factors = TRUE))
  
  surveyyear <- surveyyear - 1900
  
  variables <- union(variables, strata)
  datnew <- dat[, variables] 
  
  
  datnew$dob <- dat$b3
  datnew$survey_year <- surveyyear
  datnew$obsStart <- dat$b3
  datnew$dod <- dat$b3 + dat$b7
  datnew$obsStop <- dat$v008
  datnew$obsStop[dat$b5 == "no"] <- datnew$dod[dat$b5 == "no"]
  datnew$died <- (dat$b5 == "no")
  
  datnew$obsStop[datnew$obsStart == datnew$obsStop] <- datnew$obsStop[datnew$obsStart == datnew$obsStop] + 0.01
  
  datnew$id <- 1:nrow(datnew)
  
  
  Surv <- survival::Surv
  test <- survival::survSplit(Surv(time = obsStart, time2=obsStop, event = died, origin=dob)~dob+survey_year+died+id + caseid + v001 + v002 + v004 + v005 + v021 + v022 + v023 + v024 + v025 + v139 + bidx,
                              data = datnew, cut = c(0:60), 
                              start = "agemonth", end = "tstop", event = "died")
  test$obsStart <- test$dob + test$agemonth
  test$obsStop <- test$dob + test$tstop
  
  test$obsmonth <- test$dob + test$agemonth
  
  test$year <- floor((test$obsmonth-1)/12)
  
  test <- test[test$agemonth<60,]
  
  test <- test[test$year>80,]
  test <- test[test$year<test$survey_year,]
  
  test$tstop <- NULL
  
  test$ageGrp6 <- floor(test$agemonth/6)*6
  test$ageGrp12 <- floor(test$agemonth/12)*12
  test$ageGrpD <- test$ageGrp12
  test$ageGrpD[test$agemonth == 0] <- 0
  test$ageGrpD[test$agemonth >= 1 & test$agemonth < 12] <- 1
  
  test$ageGrp6 <- factor(test$ageGrp6)
  levels(test$ageGrp6) <- c("0-5", "6-11", "12-17", "18-23", "24-29", "30-35", "36-41", "42-47", "48-53", "54-59")
  
  test$ageGrp12 <- factor(test$ageGrp12)
  levels(test$ageGrp12) <- c("0-11", "12-23", "24-35", "36-47", "48-59")
  
  test$ageGrpD <- factor(test$ageGrpD)
  levels(test$ageGrpD) <- c("0","1-11","12-23","24-35","36-47","48-59")
  
  
  test$per5 <- "80-84"
  test$per5[test$year > 84] <- "85-89"
  test$per5[test$year > 89] <- "90-94"
  test$per5[test$year > 94] <- "95-99"
  test$per5[test$year > 99] <- "00-04"
  test$per5[test$year > 104] <- "05-09"
  test$per5[test$year > 109] <- "10-14"
  test$per5 <- factor(test$per5, levels = c("80-84","85-89","90-94","95-99","00-04","05-09","10-14"))
  
  if(length(strata) == 0){
    test$strata <- NA
  }else if(length(strata) == 1){
    test$strata <- test[, strata]
  }else{
    test$strata <- do.call(paste, c(test[strata], sep="."))
  }
  test$survey_year <- test$survey_year + 1990
  return(test)
}