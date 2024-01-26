#' Map projection for Kenya
#' 
#' Projection specifically chosen for Kenya. Project from lat/lon to northing/easting 
#' in kilometers.  Uses epsg=21097 with km units. May not work on all systems due to 
#' differences in the behavior between different PROJ and GDAL versions.
#' 
#' @param lon either longitude or, if inverse == TRUE, easting in km
#' @param lat either latitude or, if inverse == TRUE, northing in km
#' @param inverse if FALSE, projects from lon/lat to easting/northing.  Else from easting/northing to lon/lat
#' @return A 2 column matrix of easting/northing coordinates in km if inverse == FALSE. Otherwise, a 2 column matrix of longitude/latitude coordinates.
#' @author John Paige
#' @examples
#' eastLim = c(-110.6405, 832.4544)
#' northLim = c(-555.1739, 608.7130)
#' coordMatrixEN = cbind(eastLim, northLim)
#' coordMatrixLL = projKenya(coordMatrixEN, inverse=TRUE)
#' 
#' coordMatrixLL
#' # if the coordMatrixLL isn't the following, projKenya may not support 
#' # your installation of GDAL and/or PROJ:
#' #      east north
#' # [1,] 33.5  -5.0
#' # [2,] 42.0   5.5
#' 
#' projKenya(coordMatrixLL, inverse=FALSE)
#' # regardless of your PROJ/GDAL installations, the result of the 
#' # above line of could should be:
#' #            lon       lat
#' # [1,] -110.6405 -555.1739
#' # [2,]  832.4544  608.7130
#' 
#' @importFrom terra gdal
#' @importFrom sp SpatialPoints
#' @importFrom sp CRS
#' @export
projKenya = function(lon, lat=NULL, inverse=FALSE) {
  if(is.null(lat)) {
    lat = lon[,2]
    lon = lon[,1]
  }
  
  # determine version of PROJ
  ver = terra::gdal(lib="proj")
  PROJ6 <- as.numeric(substr(ver, 1, 1)) >= 6
  
  if(!inverse) {
    # from lon/lat coords to easting/northing
    if(!PROJ6) {
      lonLatCoords = sp::SpatialPoints(cbind(lon, lat), proj4string=sp::CRS("+proj=longlat"))
      coordsEN = sp::spTransform(lonLatCoords, sp::CRS("+init=epsg:21097 +units=m"))
    } else {
      lonLatCoords = sp::SpatialPoints(cbind(lon, lat), proj4string=sp::CRS(SRS_string="EPSG:4326"))
      coordsEN = sp::spTransform(lonLatCoords, sp::CRS(SRS_string="EPSG:21097"))
    }
    
    
    out = attr(coordsEN, "coords")
    colnames(out) = c("east", "north")
    
    # convert coordinates from m to km
    out = out/1000
  }
  else {
    # from easting/northing coords to lon/lat
    
    # first convert from km to m
    east = lon*1000
    north = lat*1000
    
    
    if(!PROJ6) {
      coordsEN = sp::SpatialPoints(cbind(east, north), proj4string=sp::CRS("+init=epsg:21097 +units=m"))
      lonLatCoords = sp::spTransform(coordsEN, sp::CRS("+proj=longlat"))
    } else {
      coordsEN = sp::SpatialPoints(cbind(east, north), proj4string=sp::CRS(SRS_string="EPSG:21097"))
      lonLatCoords = sp::spTransform(coordsEN, sp::CRS(SRS_string="EPSG:4326"))
    }
    
    out = attr(lonLatCoords, "coords")
    colnames(out) = c("lon", "lat")
  }
  
  out
}