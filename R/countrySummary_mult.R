#' Obtain the Horvitz-Thompson direct estimates and standard errors using delta method for multiple surveys.
#'
#' 
#' @param births A list of child-month data from multiple surveys from \code{\link{getBirths}}. The name of the list is used as the identifier in the output.
#' @param years String vector of the year intervals used
#' @param regionVar Variable name for region, typically 'v024', for older surveys might be 'v101'
#' @param timeVar Variable name for the time period indicator in the input births data.
#' @param ageVar Variable name for age group. This variable need to be in the form of "a-b" where a and b are both ages in months. For example, "1-11" means age between 1 and 11 months, including both end points. An exception is age less than one month can be represented by "0" or "0-0".
#' @param weightsVar Variable name for sampling weights, typically 'v005'
#' @param clusterVar Variable name for the IDs in the second-stage cluster sampling, typically '~v001 + v002', i.e., the cluster number and household number. When no cluster sampling design exists, this variable usually is the household ID.
#' @param Ntrials Variable for the total number of person-months if the input data (births) is in the compact form.
#' @param geo.recode The recode matrix to be used if region name is not consistent across different surveys. See \code{\link{ChangeRegion}}.
#' @param national.only Logical indicator to obtain only the national estimates
#'
#' @return This is the extension to the \code{\link{getDirect}} function that returns estimates from multiple surveys. Additional columns in the output (survey and surveyYears) specify the estimates from different surveys.
#' @seealso \code{\link{getDirect}}
#' @author Zehang Richard Li, Bryan Martin, Laina Mercer
#' @references Li, Z., Hsiao, Y., Godwin, J., Martin, B. D., Wakefield, J., Clark, S. J., & with support from the United Nations Inter-agency Group for Child Mortality Estimation and its technical advisory group. (2019). \emph{Changes in the spatial distribution of the under-five mortality rate: Small-area analysis of 122 DHS surveys in 262 subregions of 35 countries in Africa.} PloS one, 14(1), e0210645.
#' @references Mercer, L. D., Wakefield, J., Pantazis, A., Lutambi, A. M., Masanja, H., & Clark, S. (2015). \emph{Space-time smoothing of complex survey data: small area estimation for child mortality.} The annals of applied statistics, 9(4), 1889.
#' @examples
#' \dontrun{
#' data(DemoData)
#' years <- c("85-89", "90-94", "95-99", "00-04", "05-09", "10-14")
#' mean <- getDirectList(births = DemoData, years = years, 
#' regionVar = "region", timeVar = "time", clusterVar = "~clustid+id", 
#' ageVar = "age", weightsVar = "weights", geo.recode = NULL)
#' }
#' @export
getDirectList <- function(births, years,  regionVar = "region", timeVar = "time", clusterVar = "~v001+v002", ageVar = "age", weightsVar = "v005", Ntrials = NULL, geo.recode = NULL, national.only = FALSE) {
    if (length(births) == 1) {
        stop("No multiple surveys detected. Use getDirect.")
    }
    # vector of strings
    survey_years <- names(births)
    # number of surveys
    n_surv <- length(births)
    
    out_mult <- getDirect(births = births[[1]], years = years, regionVar = regionVar, timeVar = timeVar, 
        ageVar = ageVar, weightsVar = weightsVar, clusterVar = clusterVar, Ntrials = Ntrials, geo.recode = geo.recode, national.only = national.only)
    out_mult$survey <- survey_years[1]
    for (i in 2:n_surv) {
        temp <- getDirect(births = births[[i]], years = years,  regionVar = regionVar, timeVar = timeVar, 
            ageVar = ageVar, weightsVar = weightsVar, clusterVar = clusterVar, Ntrials = Ntrials, geo.recode = geo.recode, national.only = national.only)
        temp$survey <- survey_years[i]
        out_mult <- rbind(out_mult, temp)
    }
    out_mult$surveyYears <- out_mult$survey
    out_mult$survey <- as.numeric(as.factor(out_mult$survey))
    return(out_mult)
}

