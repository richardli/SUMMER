#' Functions for getting spatial polygon areas
#' 
#' Calculates the spatial area, in km^2, of every polygon in the spatialPolygonsDataFrame
#' 
#' @param mapDat SpatialPolygonsDataFrame object with map information
#' 
#' @return a vector of areas in km^2 of length \code{nrow(mapDat@data)}
#' 
#' @author John Paige
#' 
#' @examples 
#' \dontrun{
#' # download Kenya GADM shapefiles from SUMMERdata github repository
#' githubURL <- "https://github.com/paigejo/SUMMERdata/blob/main/data/kenyaMaps.rda?raw=true"
#' tempDirectory = "~/"
#' mapsFilename = paste0(tempDirectory, "/kenyaMaps.rda")
#' if(!file.exists(mapsFilename)) {
#'   download.file(githubURL,mapsFilename)
#' }
#' 
#' # load it in
#' out = load(mapsFilename)
#' out
#' adm1@data$NAME_1 = as.character(adm1@data$NAME_1)
#' adm1@data$NAME_1[adm1@data$NAME_1 == "Trans Nzoia"] = "Trans-Nzoia"
#' adm1@data$NAME_1[adm1@data$NAME_1 == "Elgeyo-Marakwet"] = "Elgeyo Marakwet"
#' }
#' 
#' @importFrom maptools readShapePoly
#' @importFrom sp CRS
#' @importFrom sp spTransform
#' 
#' @export
getSpatialAreas = function(mapDat) {
  if(!grepl("longlat", mapDat@proj4string)) {
    stop("getSpatialAreas assumes mapDat is in longitude latitude format, not another projection")
  }
  midLon = mean(mapDat@bbox[1,])
  midLat = mean(mapDat@bbox[2,])
  proj4string = paste0("+proj=laea +x_0=0 +y_0=0 +lon_0=", midLon, " +lat_0=", midLat, " +units=km")
  
  P4S.latlon <- sp::CRS("+proj=longlat +datum=WGS84")
  hrr.shp <- maptools::readShapePoly("HRR_Bdry", verbose=TRUE, proj4string=P4S.latlon)
  hrr.shp.2 <- sp::spTransform(hrr.shp, sp::CRS("+init=epsg:26978"))
  
}