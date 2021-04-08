#' Maps of Kenya Administrative Areas
#' 
#' Maps of Kenya Admin1 and Admin2 along with provinces (Admin1 areas from before the 2009 census)
#' 
#' @usage data(kenyaMaps)
#'
#' @format A number of SpatialPolygonsDataFrame objects and polygons with information about respective administrative areas for mapping. 
#' Also a spatial triangular mesh over Kenya east/north coordinates from [projKenya()] for the SPDE model. The 
#' dataset names are: adm2Kenya, adm1Kenya, kenyaMesh, kenyaPoly, provinceMapKenya.
#' 
#' @docType data
#' @name kenyaMaps
"adm1Kenya"

#' @rdname kenyaMaps
"adm2Kenya"

#' @rdname kenyaMaps
"kenyaPoly"

#' @rdname kenyaMaps
"provinceMapKenya"

#' @rdname kenyaMaps
"kenyaMesh"