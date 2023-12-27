#' Determines which administrative areas contain the given points
#' 
#' For any points not in an area, they are assigned the nearest area using 
#' fields::fields.rdist.near or fields::rdist depending on the number of points 
#' and the maximum memory in bytes with a warning.
#' 
#' @param pts 2 column matrix of lon/lat coordinates
#' @param shapefile A SpatialPolygonsDataFrame object
#' @param areaNameVar The column name in \code{slot(shapefile, "data")} 
#' corresponding to the area level of interest
#' @param delta Argument passed to fields::fields.rdist.near in fields package
#' @param mean.neighbor Argument passed to fields::fields.rdist.near in fields 
#' package
#' @param maxBytes Maximum allowed memory in bytes (default is 3Gb). Determines 
#' whether to call fields::fields.rdist.near which saves memory but requires 
#' delta and mean.neighbor inputs to be specified for fields::fields.rdist.near
#' 
#' @details delta and mean.neighbor arguments only used when some points 
#' are not in areas, perhaps due to inconsistencies in shapefiles.
#' 
#' @author John Paige
#' 
#' @seealso \code{\link{projKenya}}, \code{\link[fields]{fields.rdist.near}}
#' 
#' @examples 
#' \dontrun{
#' # download Kenya GADM shapefiles from SUMMERdata github repository
#' githubURL <- "https://github.com/paigejo/SUMMERdata/blob/main/data/kenyaMaps.rda?raw=true"
#' download.file(githubURL,"kenyaMaps.rda")
#' 
#' # load it in
#' load("kenyaMaps.rda")
#' 
#' # use the shapefile data to see what Admin1 and 2 areas the 
#' # points (0, 37) and (0.5, 38) are in
#' # (these are longitude/latitude coordinates)
#' pts = cbind(c(37, 38), c(0, .5))
#' head(slot(adm1, "data"))
#' admin1Areas = getAreaName(pts, adm1, "NAME_1")
#' admin2Areas = getAreaName(pts, adm2, "NAME_2")
#' }
#' 
#' @return A list of area IDs, area names, whether or not 
#' points are in multiple areas, and whether or not points 
#' are in no areas and assigned to the nearest one.
#' 
#' @importFrom fields fields.rdist.near
#' @importFrom fields rdist
#' @importFrom fields in.poly
#' @export
getAreaName = function(pts, shapefile, areaNameVar="NAME_1", 
                       delta=.05, mean.neighbor=50, maxBytes=3*2^30) {
  # require(fields)
  
  if(nrow(pts) == 1) {
    # must work around the fact that fields::in.poly() doesn't work for matrix with 1 row
    pts = rbind(pts, pts)
    oneRow = TRUE
  } else {
    oneRow = FALSE
  }
  
  # data(kenyaMaps, envir = environment())
  areaNames = shapefile@data[[areaNameVar]]
  
  # get area map polygons and set helper function for testing if pts are in the constituencies
  polys = shapefile@polygons
  inRegion = function(i) {
    countyPolys = polys[[i]]@Polygons
    inside = sapply(1:length(countyPolys), function(x) {
      fields::in.poly(pts, countyPolys[[x]]@coords, inflation=0)
      })
    insideAny = apply(inside, 1, any)
    
    return(insideAny*i)
  }
  out = sapply(1:length(polys), inRegion)
  if(oneRow) {
    out = matrix(out[1,], nrow=1)
  }
  multipleRegs = apply(out, 1, function(vals) {sum(vals != 0) > 1})
  areaID = apply(out, 1, function(vals) {match(1, vals != 0)})
  areaNameVec = areaNames[areaID]
  
  # for all points not in a area polygon, determine the nearest area
  insideAny = apply(out, 1, function(x) {any(x != 0)})
  if(any(!insideAny)) {
    warning("Some points not inside any areas. Assigning them to nearest area")
    problemPointsI = which(!insideAny)
    
    bytesUsed = length(problemPointsI) * nrow(pts)
    if(bytesUsed > maxBytes) {
      # get nearby points (points within delta lon/lat units), remove self matches
      nearbyPoints = fields::fields.rdist.near(matrix(pts[problemPointsI,], ncol=2), pts, 
                                               delta=delta, mean.neighbor=mean.neighbor)
      selfI = nearbyPoints$ra == 0
      nearbyPoints$ind = nearbyPoints$ind[!selfI,]
      nearbyPoints$ra = nearbyPoints$ra[!selfI]
      nearbyI = lapply(sort(unique(nearbyPoints$ind[,1])), function(x) {
        nearbyPoints$ind[nearbyPoints$ind[,1] == x,2]
        })
    } else {
      # get all points, remove self matches
      dists = fields::rdist(matrix(pts[problemPointsI,], ncol=2), pts)
      dists[cbind(1:length(problemPointsI), problemPointsI)] = Inf
      nearbyI = apply(dists, 1, which.min)
    }
    
    # get nearby constituencies, counties, and distances
    nearbyAreas = lapply(nearbyI, function(x) {areaNameVec[x]})
    nearbyLengths = sapply(nearbyI, function(x) {length(x)})
    
    if(bytesUsed > maxBytes) {
      nearbyDistances = c()
      startI = 1
      for(i in 1:length(nearbyI)) {
        endI = startI + nearbyLengths[i] - 1
        nearbyDistances = c(nearbyDistances, list(nearbyPoints$ra[startI:endI]))
        startI = endI + 1
      }
    } else {
      nearbyDistances = dists[cbind(1:nrow(dists), nearbyI)]
    }
    
    
    # sort nearby constituencies and indices by distance
    for(i in 1:length(nearbyI)) {
      thisDistances = nearbyDistances[[i]]
      sortI = sort(thisDistances, index.return=TRUE)$ix
      nearbyDistances[[i]] = nearbyDistances[[i]][sortI]
      nearbyAreas[[i]] = nearbyAreas[[i]][sortI]
      nearbyI[[i]] = nearbyI[[i]][sortI]
    }
    
    # get nearest non-NA area and assign it
    closestArea = sapply(nearbyAreas, function(x) {x[match(TRUE, !is.na(x))]})
    areaNameVec[problemPointsI] = closestArea
  }
  
  list(areaID=areaID, areaNames=areaNameVec, 
       multipleRegs=multipleRegs, notInAnyAreas=!insideAny)
}