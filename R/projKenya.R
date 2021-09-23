#' Map projection for Kenya
#' 
#' Projection specifically chosen for Kenya. Project from lat/lon to UTM northing/easting 
#' in kilometers.  Uses epsg=21097
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
#' @importFrom rgdal rgdal_extSoftVersion
#' @importFrom sp SpatialPoints
#' @importFrom sp CRS
#' @export
projKenya = function(lon, lat=NULL, inverse=FALSE) {
  if(is.null(lat)) {
    lat = lon[,2]
    lon = lon[,1]
  }
  
  # determine version of PROJ
  ver = rgdal::rgdal_extSoftVersion()
  theseNames = names(ver)
  thisI = which(grepl("PROJ", theseNames))
  PROJ6 <- as.numeric(substr(ver[thisI], 1, 1)) >= 6
  if(PROJ6) {
    warning("projKenya has not been tested with PROJ.4 version >= 6")
  }
  
  if(!inverse) {
    # from lon/lat coords to easting/northing
    if(!PROJ6) {
      lonLatCoords = sp::SpatialPoints(cbind(lon, lat), proj4string=sp::CRS("+proj=longlat"))
      coordsUTM = sp::spTransform(lonLatCoords, sp::CRS("+init=epsg:21097 +units=km"))
    } else {
      lonLatCoords = sp::SpatialPoints(cbind(lon, lat), proj4string=sp::CRS(SRS_string="EPSG:4326"))
      coordsUTM = sp::spTransform(lonLatCoords, sp::CRS("+init=epsg:21097 +units=km"))
    }
    out = attr(coordsUTM, "coords")
    colnames(out) = c("lon", "lat")
  }
  else {
    # from easting/northing coords to lon/lat
    east = lon
    north = lat
    if(!PROJ6) {
      coordsUTM = sp::SpatialPoints(cbind(east, north), proj4string=sp::CRS("+init=epsg:21097 +units=km"))
      lonLatCoords = sp::spTransform(coordsUTM, sp::CRS("+proj=longlat"))
    } else {
      coordsUTM = sp::SpatialPoints(cbind(east, north), proj4string=sp::CRS("+init=epsg:21097 +units=km"))
      lonLatCoords = sp::spTransform(coordsUTM, sp::CRS(SRS_string="EPSG:4326"))
    }
    
    out = attr(lonLatCoords, "coords")
  }
  
  out
}