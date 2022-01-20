#' Kenya 2009 Census Frame and Related Datasets
#' 
#' Datasets related to the 2009 census frame for Kenya based on the 2009 Kenya Population and Housing 
#' Census. General population totals are estimated for 2014. Based on 2014 population density 
#' estimates interpolated with exponential growth rate between 2010 and 2015 from WorldPop data.
#'
#' @format A number of data.frames with information about the 2009 Kenya Population and Housing Census and the population in Kenya at the time of the 2014 Demographic Health Survey. Some of the 
#' data.frames have been adjusted to contain information about neonatals born from 2010-2014 rather than general population in 2014. 
#' The dataset names are: easpaKenya, easpaKenyaNeonatal, 
#' poppaKenya, and poppsubKenya.
#' @references Kenya National Bureau of Statistics, Ministry of Health/Kenya, National AIDS Control Council/Kenya, Kenya Medical Research Institute, and National Council For Population And Development/Kenya, 2015. Kenya Demographic and Health Survey 2014. Rockville, Maryland, USA. URL: http://dhsprogram.com/pubs/pdf/FR308/FR308.pdf.
#' @references Stevens, F.R., Gaughan, A.E., Linard, C., Tatem, A.J., 2015. Disaggregat- ing census data for population mapping using random forests with remotely-sensed and ancillary data. PloS One 10, e0107042.
#' @references Tatem, A.J., 2017. WorldPop, open data for spatial demography. Scientific Data 4.
#' @source <http://dhsprogram.com/pubs/pdf/FR308/FR308.pdf>
#' @docType data
#' @name kenyaPopulationData
#' @usage data(kenyaPopulationData)
"easpaKenya"

#' @rdname kenyaPopulationData
"easpaKenyaNeonatal"

#' @rdname kenyaPopulationData
"poppaKenya"

#' @rdname kenyaPopulationData
"poppsubKenya"