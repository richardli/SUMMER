#' Maps of Kenya Administrative Areas
#' 
#' Maps of Kenya Admin1 and Admin2 along with provinces (Admin1 areas from before the 2009 census)
#' 
#' @usage data(kenyaMaps)
#'
#' @format A number of SpatialPolygonsDataFrame objects and polygons with information about respective administrative areas for mapping. 
#' Also a spatial triangular mesh over Kenya east/north coordinates from \code{\link{projKenya}} for the SPDE model. The 
#' dataset names are: adm2Kenya, adm1Kenya, kenyaMesh, kenyaPoly, provinceMapKenya.
#' 
#' @details In the Admin2 level shapefile we combine Admin2 areas in Kenya that have spatial area less than 50 km^2 with 
#' neighboring ones in the same Admin1 area in order to reduce numerical approximation errors caused 
#' by representing small areas using the pixel grid, resulting in 273 (possibly combined) Admin2 areas 
#' in total out of the original 290. The names of the areas that are combined are included in the new 
#' area name, separated by the "+" symbol.
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