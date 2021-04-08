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
#' @export
projKenya = function(lon, lat=NULL, inverse=FALSE) {
  if(is.null(lat)) {
    lat = lon[,2]
    lon = lon[,1]
  }
  
  if(!inverse) {
    # from lon/lat coords to easting/northing
    lonLatCoords = sp::SpatialPoints(cbind(lon, lat), proj4string=sp::CRS("+proj=longlat"))
    coordsUTM = sp::spTransform(lonLatCoords, sp::CRS("+init=epsg:21097 +units=km"))
    out = attr(coordsUTM, "coords")
    colnames(out) = c("lon", "lat")
  }
  else {
    # from easting/northing coords to lon/lat
    east = lon
    north = lat
    coordsUTM = sp::SpatialPoints(cbind(east, north), proj4string=sp::CRS("+init=epsg:21097 +units=km"))
    lonLatCoords = sp::spTransform(coordsUTM, sp::CRS("+proj=longlat"))
    out = attr(lonLatCoords, "coords")
  }
  
  out
}

#' Utility functions for Kenya administrative areas
#' 
#' @param pts 2 column matrix of lon/lat (if project == FALSE) or east/north (if project == TRUE) coordinates
#' @param project project to longitude/latitude coordinates
#' @param delta argument passed to fields.rdist.near in fields package
#' @param mean.neighbor argument passed to fields.rdist.near in fields package
#' @param constituencyNames character vector of Kenya constituency (Admin1) 
#' names to convert to county (Admin2) names using adm2Kenya
#' 
#' @author John Paige
#' 
#' @seealso \code{\link{projKenya}}
#' @name kenyaMapUtilities
NULL

#' @describeIn kenyaMapUtilities Calculates Admin2 areas of the input 
#' spatial coordinates using adm2Kenya variable from kenyaMaps dataset. 
#' NOTE: this calls \code{data(kenyaMaps)}
#' 
#' @export
getAdmin2Kenya = function(pts, project=FALSE, delta=.05, mean.neighbor=50) {
  # require(fields)
  
  # project pts to lon/lat coordinate system if user specifies
  if(project)
    pts = projKenya(pts, inverse=TRUE)
  
  # data(kenyaMaps, envir = environment())
  constituencyNames = SUMMER::adm2Kenya@data$CONSTITUEN
  
  # get constituency map polygons and set helper function for testing if pts are in the constituencies
  polys = SUMMER::adm2Kenya@polygons
  inRegion = function(i) {
    countyPolys = polys[[i]]@Polygons
    inside = sapply(1:length(countyPolys), function(x) {fields::in.poly(pts, countyPolys[[x]]@coords, inflation=0)})
    insideAny = apply(inside, 1, any)
    
    return(insideAny*i)
  }
  out = sapply(1:length(polys), inRegion)
  multipleRegs = apply(out, 1, function(vals) {sum(vals != 0) > 1})
  constituencyID = apply(out, 1, function(vals) {match(1, vals != 0)})
  constituencyNameVec = constituencyNames[constituencyID]
  
  # for all points not in a constituency polygon, determine the nearest constituency
  insideAny = apply(out, 1, function(x) {any(x != 0)})
  if(any(!insideAny)) {
    problemPointsI = which(!insideAny)
    
    # get nearby points (points within .2 lon/lat units), remove self matches
    nearbyPoints = fields::fields.rdist.near(pts[problemPointsI,], pts, delta=delta, mean.neighbor=mean.neighbor)
    selfI = nearbyPoints$ra == 0
    nearbyPoints$ind = nearbyPoints$ind[!selfI,]
    nearbyPoints$ra = nearbyPoints$ra[!selfI]
    nearbyI = lapply(sort(unique(nearbyPoints$ind[,1])), function(x) {nearbyPoints$ind[nearbyPoints$ind[,1] == x,2]})
    
    # get nearby constituencies, counties, and distances
    nearbyConstituencies = lapply(nearbyI, function(x) {constituencyNameVec[x]})
    nearbyLengths = sapply(nearbyI, function(x) {length(x)})
    nearbyDistances = c()
    # nearbyCounties = c()
    startI = 1
    for(i in 1:length(nearbyI)) {
      endI = startI + nearbyLengths[i] - 1
      nearbyDistances = c(nearbyDistances, list(nearbyPoints$ra[startI:endI]))
      startI = endI + 1
    }
    
    # sort nearby constituencies and indices by distance
    for(i in 1:length(nearbyI)) {
      thisDistances = nearbyDistances[[i]]
      sortI = sort(thisDistances, index.return=TRUE)$ix
      nearbyDistances[[i]] = nearbyDistances[[i]][sortI]
      nearbyConstituencies[[i]] = nearbyConstituencies[[i]][sortI]
      nearbyI[[i]] = nearbyI[[i]][sortI]
    }
    
    # get nearest non-NA constituency and assign it
    closestConstituency = sapply(nearbyConstituencies, function(x) {x[match(TRUE, !is.na(x))]})
    
    constituencyNameVec[problemPointsI] = closestConstituency
  }
  
  list(constituencyID=constituencyID, constituencyNames=constituencyNameVec, multipleRegs=multipleRegs)
}

#' @describeIn kenyaMapUtilities Computes what administrative regions the given points are in
#' @export
getAdmin1Kenya = function(pts, project=FALSE, delta=.05, mean.neighbor=50) {
  out = getAdmin2Kenya(pts, project, delta, mean.neighbor)
  constituencies = out$constituencyNames
  admin2ToAdmin1Kenya(constituencies)
}

#' @describeIn kenyaMapUtilities Gives Kenya Admin1 area containing Admin2 area
#' @export
admin2ToAdmin1Kenya = function(constituencyNames) {
  constituencies = SUMMER::adm2Kenya@data$CONSTITUEN
  counties = SUMMER::adm2Kenya@data$COUNTY_NAM
  matchI = match(constituencyNames, constituencies)
  counties[matchI]
}