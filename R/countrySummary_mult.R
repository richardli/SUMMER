#' Obtain the Horvitz-Thompson direct estimates and standard errors using delta method for multiple surveys.
#'
#' 
#' @param births A list of child-month data from multiple surveys from \code{\link{getBirths}}. The name of the list is used as the identifier in the output.
#' @param years String vector of the year intervals used
#' @param idVar Variable name for ID, typically 'v002'
#' @param regionVar Variable name for region, typically 'v024', for older surveys might be 'v101'
#' @param timeVar Variable name for time, typically 'per5'
#' @param ageVar Variable name for age group, default assumes the variable is called 'ageGrpD'
#' @param weightsVar Variable name for sampling weights, typically 'v005'
#' @param clusterVar Variable name for cluster, typically '~v001 + v002'
#' @param geo.recode The recode matrix to be used if region name is not consistent across different surveys. See \code{\link{ChangeRegion}}.
#'
#' @return a matrix of period-region summary of the Horvitz-Thompson direct estimates, the standard errors using delta method for a single survey, the 95\% confidence interval, the logit of the estimates, and the survey labels.
#' @seealso \code{\link{countrySummary}}
#' @examples
#' \dontrun{
#' data(Uganda)
#' years <- c("85-89", "90-94", "95-99", "00-04", "05-09", "10-14")
#' u5m <- countrySummary_mult(births = Uganda, years = years, idVar = "id", 
#' regionVar = "region", timeVar = "time", clusterVar = "~clustid+id", 
#' ageVar = "age", weightsVar = "weights", geo.recode = NULL)
#' }
#' @export
countrySummary_mult <- function(births, years, idVar = "v002", regionVar = "region", timeVar = "per5", clusterVar = "~v001+v002",
                                ageVar = "ageGrpD", weightsVar = "v005", geo.recode = NULL) {
    if (length(births) == 1) {
        stop("No multiple surveys detected. Use countrySummary.")
    }
    # vector of strings
    survey_years <- names(births)
    # number of surveys
    n_surv <- length(births)
    
    out_mult <- countrySummary(births = births[[1]], years = years, idVar = idVar, regionVar = regionVar, timeVar = timeVar, 
        ageVar = ageVar, weightsVar = weightsVar, clusterVar = clusterVar, geo.recode = geo.recode)
    out_mult$survey <- survey_years[1]
    for (i in 2:n_surv) {
        temp <- countrySummary(births = births[[i]], years = years, idVar = idVar, regionVar = regionVar, timeVar = timeVar, 
            ageVar = ageVar, weightsVar = weightsVar, clusterVar = clusterVar, geo.recode = geo.recode)
        temp$survey <- survey_years[i]
        out_mult <- rbind(out_mult, temp)
    }
    out_mult$surveyYears <- out_mult$survey
    out_mult$survey <- as.numeric(as.factor(out_mult$survey))
    return(out_mult)
}

