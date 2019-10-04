#' Auxiliary data for Kenya 2014 DHS. 
#' 
#' The list contains several data frames.  
#' 
#' 
#' \itemize{
#'  \item HIV2014, a data frame with three columns: time (in five year periods),  region (8 Admin-1 region groups), and the estimated bias of the reported U5MR due to HIV for each 5 year period from 1990-1994 to 2010-2014. The bias is represented as the ratio of the reported U5MR to the true U5MR.
#'  \item HIV2014.yearly, a data frame with three columns: time (in five year periods),  region (8 Admin-1 region groups), and the estimated bias of the reported U5MR due to HIV for each year from 1980 to 2014. The bias is represented as the ratio of the reported U5MR to the true U5MR.
#'  \item IGME2019. Yearly Estimates of national under-5 child mortality in Kenya from the 2019 UN-IGME estimates.
#' \item UrbanProp. Proportion of urban population by county and total population by county. Source: 2009 Kenya Population and Housing Census, and Table A2 of Kenya 2014 DHS report. 
#' }  
#' @references Neff Walker, Kenneth Hill, and Fengmin Zhao (2012) \emph{Child mortality estimation: methods used to adjust for bias due to aids in estimating trends in under-five mortality.}, \cr \emph{PLoS Medicine, 9(8):e1001298}.
#' 
#' @docType data
#'
#' @usage data(KenData)
#' 
"KenData"