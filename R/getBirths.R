#' Reformat full birth records into person-month format
#' 
#'
#' @param filepath file path of raw .dta file from DHS. Only used when data frame is not provided in the function call.
#' @param data data frame of a DHS survey
#' @param surveyyear year of survey. Observations after this year will be excluded from the analysis.
#' @param variables vector of variables to be used in obtaining the person-month files. The variables correspond the the DHS recode manual VI. For early DHS data, the variable names may need to be changed.
#' @param strata vector of variable names used for strata. If a single variable is specified, then that variable will be used as strata indicator If multiple variables are specified, the interaction of these variables will be used as strata indicator. 
#' @param dob variable name for the date of birth.
#' @param alive variable name for the indicator of whether child was alive or dead at the time of interview.
#' @param age variable name for the age at death of the child in completed months.
#' @param age.truncate the smallest age in months where only full years are reported. The default value is 24, which corresponds to the DHS practice of recording only age in full years for children over 2 years old. That is, for children with age starting from 24 months old, we assume the age variable reported in multiples of 12 are truncated from its true value. For example, children between age 24 to 35 months are all recorded as 24. To account for the truncation of age, 5 months are added to all ages recorded in multiples of 12 starting from 24. To avoid this adjustment, set this argument to NA. 
#' @param date.interview variable name for the date of interview.
#' @param month.cut the cutoff of each bins of age group in the unit of months. Default values are 1, 12, 24, 36, 48, and 60, representing the age groups (0, 1), [1, 12), [12, 24), ..., [48, 60).
#' @param year.cut The cutoff of each bins of time periods, including both boundaries. Default values are 1980, 1985, ..., 2020, representing the time periods 80-84, 85-89, ..., 15-19. Notice that if each bin contains one year, the last year in the output is max(year.cut)-1. For example, if year.cut = 1980:2020, the last year in the output is 2019.
#' @param min.last.period The cutoff for how many years the last period must contain in order to be counted in the output. For example, if the last period is 2015-2019 and min.last.period = 3, person-months for the last period will only be returned if survey contains observations at least in 2017. This argument avoids the situation that estimates for the last period being based on only a small number of initial years, if applicable. Default to be 0. 
#' @param cmc.adjust number of months to add to the recorded month in the dataset. Some DHS surveys does not use Gregorian calendar (the calendar used in most of the world). For example, the Ethiopian calendar is 92 months behind the Gregorian calendar in general. Then we can set cmc.adjust to 92, which adds 92 months to all dates in the dataset, effectively transforming the Ethiopian calendar to the Gregorian calendar.  
#' @param compact logical indicator of whether the compact format is returned. In the compact output, person months are aggregated by cluster, age, and time. Total number of person months and deaths in each group are returned instead of the raw person-months.
#' @param compact.by vector of variables to summarize the compact form by. 
#' 
#' @return This function returns a new data frame where each row indicate a person-month, with the additional variables specified in the function argument.
#' @author Zehang Richard Li, Bryan Martin, Laina Mercer
#' @references Li, Z., Hsiao, Y., Godwin, J., Martin, B. D., Wakefield, J., Clark, S. J., & with support from the United Nations Inter-agency Group for Child Mortality Estimation and its technical advisory group. (2019). \emph{Changes in the spatial distribution of the under-five mortality rate: Small-area analysis of 122 DHS surveys in 262 subregions of 35 countries in Africa.} PloS one, 14(1), e0210645.
#' @references Mercer, L. D., Wakefield, J., Pantazis, A., Lutambi, A. M., Masanja, H., & Clark, S. (2015). \emph{Space-time smoothing of complex survey data: small area estimation for child mortality.} The annals of applied statistics, 9(4), 1889.
#' @examples 
#' \dontrun{
#' my_fp <- "/myExampleFilepath/surveyData.DTA"
#' DemoData <- getBirths(filepath = my_fp, surveyyear = 2015) 
#' }
#' 
#' @export
getBirths <- function(filepath = NULL, data = NULL, surveyyear = NA, variables = c("caseid", "v001", "v002", "v004", "v005", "v021", "v022", "v023", "v024", "v025", "v139", "bidx"), strata=c("v024", "v025"), dob = "b3", alive = "b5", age = "b7", age.truncate = 24, date.interview= "v008", month.cut = c(1,12,24,36,48,60), year.cut=seq(1980, 2020, by=5), min.last.period = 0, cmc.adjust = 0, compact = FALSE, compact.by = c('v001',"v024", "v025", "v005")) {
  if(is.null(data)){
      dat <- suppressWarnings(readstata13::read.dta13(filepath, generate.factors = TRUE))    
  }else{
      dat <- data.frame(data)
  }
  
  period.1yr <- year.cut[2] - year.cut[1] == 1
  surveyyear <- surveyyear - 1900
  year.cut <- year.cut - 1900
  variables <- union(variables, strata)
  
  # use mid point to impute age reported in years
  if(!is.na(age.truncate)){
      trunc <- intersect(which(dat[, age] %% 12 == 0), which(dat[, age] >= age.truncate))
      if(length(trunc) > 0){
          dat[trunc, age] <- dat[trunc, age] + 5
          message(paste0("Children with age at least ", age.truncate, " months are assumed to have recorded age truncated to full years. \nRecorded age + 5 months is used to adjust for the truncation for ages >= ", age.truncate , " and are multiples of 12." ))
      }else{
          message(paste0("Children with age at least ", age.truncate, " months are specified to have recorded age truncated to full years. But no records in the data needs the adjustment." ))
      }
  }else{
    message("No age truncation adjustment is performed.")
  }

  datnew <- dat[, variables] 
  dat[, alive] <- tolower(dat[, alive])
  
  datnew$dob <- dat[, dob] + cmc.adjust
  datnew$survey_year <- surveyyear
  datnew$obsStart <- dat[, dob] + cmc.adjust
  datnew$dod <- dat[, dob] + dat[, age] + cmc.adjust 
  datnew$obsStop <- dat[, date.interview] + cmc.adjust
  datnew$obsStop[dat[, alive] == "no"] <- datnew$dod[dat[, alive] == "no"]
  datnew$died <- (dat[, alive] == "no")
  
  datnew$obsStop[datnew$obsStart == datnew$obsStop] <- datnew$obsStop[datnew$obsStart == datnew$obsStop] + 0.01
  
  datnew$id.new <- 1:nrow(datnew)
  
  # Surv(time = obsStart, time2=obsStop, event = died, origin=dob)~dob+survey_year+died+id + caseid + v001 + v002 + v004 + v005 + v021 + v022 + v023 + v024 + v025 + v139 + bidx
  formula <- as.formula(paste(c("Surv(time = obsStart, time2=obsStop, event = died, origin=dob)~dob+survey_year+died+id.new", union(variables, strata)),  collapse = "+"))

  Surv <- survival::Surv
  # test <- survival::survSplit(formula,
  #                             data = datnew, cut = c(0:max(month.cut)), 
  #                             start = "agemonth", end = "tstop", event = "died")

  test <- survival::survSplit(formula,
                              data = datnew, cut = c(0.02, 1:max(month.cut)),
                              start = "agemonth", end = "tstop", event = "died")
  # the survival package splits the data into
  # (0, 0.02], (0.02, 1], (1, 2], ..., (59, 60], ...
  # which is equivalent to 
  # [0, 1),     [1, 2),    [2, 3),..., [60, 61), ...
  test$agemonth <- test$agemonth + 1
  test$agemonth[test$agemonth == 1] <- 0   
  test$agemonth[test$agemonth == 1.02] <- 1   

  test$obsStart <- test$dob + test$agemonth
  test$obsStop <- test$dob + test$tstop
  test$obsmonth <- test$dob + test$agemonth
  test$year <- floor((test$obsmonth-1)/12)
  
  test <- test[test$agemonth<max(month.cut), ]
  test <- test[test$year>=year.cut[1], ]
  test <- test[test$year<year.cut[length(year.cut)], ]
  if(!is.na(surveyyear)) test <- test[test$year<= surveyyear, ]
  
  # remove observations if last period has no more than min.last.period years.
  # e.g., if last period 15-19, min.last.period = 3
  #       determine if max year is at least 2017
  if(min.last.period > 0){
    # year.last.start <- year.cut[length(year.cut)-1] 
    year.last.start <- max(year.cut[year.cut <= max(test$year)])
    if(max(test$year) < year.last.start + min.last.period - 1){
        test <- test[test$year < year.last.start, ]
    }
  }
  
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
  
  # test$age <- factor(test$age)
  # levels(test$age) <- bins #c("0","1-11","12-23","24-35","36-47","48-59")
  test$age <- factor(as.character(test$age), levels = as.character(bins))
  
  if(period.1yr){
      test$time <- year.cut[1] + 1900 
      year.bin <-  year.cut[1] + 1900
      if(length(year.cut) > 2){
        for(i in 2:(length(year.cut)-1)){
          test$time[test$year >= year.cut[i] & test$year < year.cut[i+1]] <- year.cut[i] + 1900
          year.bin <- c(year.bin, year.cut[i] + 1900)
        }        
      }
  }else{
      year.cut2 <- year.cut
      year.cut2[year.cut2 >= 100] <- as.character(year.cut2[year.cut2 >= 100] - 100)
      year.cut2[as.numeric(year.cut2) < 10] <- paste0("0", year.cut2[as.numeric(year.cut2) < 10])
      year.cut3 <- year.cut - 1
      year.cut3[year.cut3 >= 100] <- as.character(year.cut3[year.cut3 >= 100] - 100)
      year.cut3[as.numeric(year.cut3) < 10] <- paste0("0", year.cut3[as.numeric(year.cut3) < 10])

      test$time <- paste(year.cut2[1], year.cut3[2], sep = "-")
      year.bin <- paste(year.cut2[1], year.cut3[2], sep = "-")
      if(length(year.cut) > 2){
        for(i in 2:(length(year.cut)-1)){
          test$time[test$year >= year.cut[i] & test$year < year.cut[i+1]] <- paste(year.cut2[i], year.cut3[i+1], sep = "-")
          year.bin <- c(year.bin, paste(year.cut2[i], year.cut3[i+1], sep = "-"))
        }        
      }
  }



  test$time <- factor(test$time, levels = year.bin)
  
  if(length(strata) == 0){
    test$strata <- NA
  }else if(length(strata) == 1){
    test$strata <- test[, strata]
  }else{
    test$strata <- do.call(paste, c(test[strata], sep="."))
  }
  test$survey_year <- test$survey_year + 1900

  if(compact){
      test.comp <- test[, c(compact.by, "age", "strata", "time", "died")]
      test.comp$total <- 1
      formula <- as.formula(paste0(".~age + time + strata + ", paste(compact.by, collapse = " + ")))
      test.comp <- stats::aggregate(formula, data = test.comp, FUN = sum, drop = TRUE)
     test <- test.comp
  }


  return(test)
}