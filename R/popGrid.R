#' Functions for pixellated population and urbanicity grids
#' 
#' Generate pixellated grid of coordinates (both longitude/latitude and east/north) 
#' over spatial domain of the given resolution
#' @param popMat pixellated grid data frame with variables `area` and `pop`
#' @param poppa data.frame of population per area separated by urban/rural with variables:
#' \describe{
#'   \item{area}{name of area}
#'   \item{popUrb}{total urban (general) population of area}
#'   \item{popRur}{total rural (general) population of area}
#'   \item{popTotal}{total (general) population of area}
#'   \item{pctUrb}{percentage of population in the area that is urban (between 0 and 100)}
#' }
#' @param poppsub data.frame of population per subarea separated by 
#' urban/rural using for population density grid normalization or urbanicity 
#' classification. Often based on extra fine scale population density grid. Has variables:
#' \describe{
#'   \item{subarea}{name of subarea}
#'   \item{area}{name of area}
#'   \item{popUrb}{total urban (general) population of subarea}
#'   \item{popRur}{total rural (general) population of subarea}
#'   \item{popTotal}{total (general) population of subarea}
#'   \item{pctUrb}{percentage of population in the subarea that is urban (between 0 and 100)}
#' }
#' @param kmRes the resolution of the pixelated grid in km
#' @param domainPoly a polygon representing the full spatial domain (e.g. country)
#' @param eastLim range in km easting over the spatial domain under the input projection
#' @param northLim range in km northing over the spatial domain under the input projection
#' @param mapProjection a projection function taking longitude and latitude and returning easting and 
#'                northing in km. Or the inverse if inverse is set to TRUE. For example, 
#'                \code{\link{projKenya}}
#' @param pop population density raster
#' @param pixelGrid a pixel grid data.frame such as that created by \code{\link{generatePixelGrid}}
#' @param areas character vector of length nPixels giving area names associated with each pixel
#' @param subareas character vector of length nPixels giving subarea names associated with each pixel
#' @param regions character vector of length nPixels giving a custom set of regions for which to generate 
#'  a population frame using population density
#' @param stratifyByUrban Whether to stratify the pixellated grid by urban/rural. If TRUE, sets 
#'                  threshold on population density when classifying urban/rural for each 
#'                  area to obtain urbanicity classification based on proportion urban 
#'                  population for that area
#' @param poppaTarget target population per area stratified by urban rural. Same format as poppa
#' @param adjustBy whether to adjust population density by the `area` or `subarea` level
#' @author John Paige
#' @examples 
#' \dontrun{
#' # create 5km pixellated grid over Kenya
#' data(kenyaMaps)
#' eastLim = c(-110.6405, 832.4544)
#' northLim = c(-555.1739, 608.7130)
#' pixelGrid = generatePixelGrid(kmRes=5, kenyaPoly, eastLim, northLim, 
#'   mapProjection=projKenya)
#' 
#' # add population density to the grid, assign urbanicity by thresholding 
#' # population density based on estimated proportion population urban/rural. 
#' # Note that popKenya is some general population density raster object for Kenya
#' data(kenyaPopulationData)
#' require(raster)
#' popMat = makeInterpPopMat(popKenya, pixelGrid, 
#'   areas=popMatKenya$area, subareas=popMatKenya$subarea, 
#'   poppa=poppaKenya, poppsub=poppsubKenya, stratifyByUrban=TRUE)
#' 
#' # adjust popMat to be target rather than general population density. First
#' # created target population frame (these numbers are based on IPUMS microcensus data)
#' mothersPerHouseholdUrb = 0.3497151
#' childrenPerMotherUrb = 1.295917
#' mothersPerHouseholdRur = 0.4787696
#' childrenPerMotherRur = 1.455222
#' targetPopPerStratumUrban = easpaKenya$HHUrb * mothersPerHouseholdUrb * childrenPerMotherUrb
#' targetPopPerStratumRural = easpaKenya$HHRur * mothersPerHouseholdRur * childrenPerMotherRur
#' easpaKenyaNeonatal = easpaKenya
#' easpaKenyaNeonatal$popUrb = targetPopPerStratumUrban
#' easpaKenyaNeonatal$popRur = targetPopPerStratumRural
#' easpaKenyaNeonatal$popTotal = easpaKenyaNeonatal$popUrb + easpaKenyaNeonatal$popRur
#' easpaKenyaNeonatal$pctUrb = 100 * easpaKenyaNeonatal$popUrb / easpaKenyaNeonatal$popTotal
#' easpaKenyaNeonatal$pctTotal = 
#'   100 * easpaKenyaNeonatal$popTotal / sum(easpaKenyaNeonatal$popTotal)
#' 
#' # generate the target population density by scaling the current population density grid 
#' # at the Admin1 x urban/rural level
#' popMatKenyaNeonatal = adjustPopMat(popMatKenya, easpaKenyaNeonatal)
#' 
#' # generate population table from the population density grid
#' poppsubNeonatal = poppRegionFromPopMat(popMatKenyaNeonatal, popMatKenyaNeonatal$subareas)
#' poppsubNeonatal = cbind(subarea=poppsubNeonatal$region, 
#'   area=admin2ToAdmin1Kenya(poppsubNeonatal$region), poppsubNeonatal[,-1])
#'   print(head(poppsubNeonatal))
#' }
#' @export
generatePixelGrid = function(kmRes=5, domainPoly, eastLim, northLim, mapProjection) {
  
  # get a rectangular grid
  eastGrid = seq(eastLim[1], eastLim[2], by=kmRes)
  northGrid = seq(northLim[1], northLim[2], by=kmRes)
  utmGrid = fields::make.surface.grid(list(east=eastGrid, north=northGrid))
  
  # project coordinates into lat/lon
  lonLatGrid = mapProjection(utmGrid, inverse=TRUE)
  
  # subset grid so it's in spatial domain
  inDomain = fields::in.poly(lonLatGrid, domainPoly)
  utmGrid = utmGrid[inDomain,]
  lonLatGrid = lonLatGrid[inDomain,]
  
  data.frame(lon=lonLatGrid[,1], lat=lonLatGrid[,2], east=utmGrid[,1], north=utmGrid[,2])
}

#' @describeIn generatePixelGrid Set thresholds of population density for urbanicity classifications within each area 
#' based on that area's percent population urban. Intended as a helper function of \code{\link{makeInterpPopMat}}
#' @export
setThresholdsArea = function(popMat, poppa) {
  
  getAreaThresh = function(areaName) {
    # do the setup
    thisArea = as.character(popMat$area) == areaName
    thisPop = popMat$pop[thisArea]
    thisTot = sum(thisPop)
    pctUrb = poppa$pctUrb[poppa$area == areaName]/100
    pctRural = 1 - pctUrb
    
    if(pctUrb == 1) {
      return(-Inf)
    } else if(pctUrb == 0) {
      return(Inf)
    }
    
    # calculate threshold by integrating ecdf via sorted value cumulative sum
    sortedPop = sort(thisPop)
    cumsumPop = cumsum(sortedPop)
    threshI = match(1, cumsumPop >= thisTot*pctRural)
    if((threshI != 1) && (threshI != length(thisPop))) {
      thresh = sortedPop[threshI]
    } else {
      # make sure not all pixels are urban or all are rural
      if(threshI == 1) {
        thresh = mean(sortedPop[1], sortedPop[2])
      } else {
        thresh = mean(sortedPop[length(thisPop)], sortedPop[length(thisPop)-1])
      }
    }
    
    thresh
  }
  
  # compute threshold for each county
  areas = poppa$area
  threshes = sapply(areas, getAreaThresh)
  
  list(areas=areas, threshes=threshes)
}

#' @describeIn generatePixelGrid Set thresholds of population density for urbanicity classifications within each subarea 
#' based on that subarea's percent population urban. Intended as a helper function of \code{\link{makeInterpPopMat}}
#' @export
setThresholdsSubarea = function(popMat, poppsub) {
  
  getSubareaThresh = function(subareaName) {
    # do the setup
    thisSubarea = as.character(popMat$subarea) == subareaName
    thisPop = popMat$pop[thisSubarea]
    thisTot = sum(thisPop)
    pctUrb = poppsub$popUrb[poppsub$subarea == subareaName]/poppsub$popTotal[poppsub$subarea == subareaName]
    pctRural = 1 - pctUrb
    
    if(pctUrb == 1) {
      return(-Inf)
    } else if(pctUrb == 0) {
      return(Inf)
    }
    
    # calculate threshold by integrating ecdf via sorted value cumulative sum
    sortedPop = sort(thisPop)
    cumsumPop = cumsum(sortedPop)
    threshI = match(1, cumsumPop >= thisTot*pctRural)
    if((threshI != 1) && (threshI != length(thisPop))) {
      thresh = sortedPop[threshI]
    } else {
      # make sure not all pixels are urban or all are rural
      if(threshI == 1) {
        thresh = mean(c(sortedPop[1], sortedPop[2]))
      } else {
        thresh = mean(c(sortedPop[length(thisPop)], sortedPop[length(thisPop)-1]))
      }
    }
    
    thresh
  }
  
  # compute threshold for each county
  subareas = poppsub$subarea
  threshes = sapply(subareas, getSubareaThresh)
  
  list(subareas=subareas, threshes=threshes)
}

#' @describeIn generatePixelGrid Generate the population density surface along with urbanicity classifications
#' @importFrom raster extract
#' @export
makeInterpPopMat = function(pop, pixelGrid, areas, subareas, poppa, poppsub=NULL, 
                            stratifyByUrban=TRUE) {
  
  # extract information about pixellated grid
  lonLatGrid = cbind(pixelGrid$lon, pixelGrid$lat)
  
  # get population density at those coordinates
  interpPopVals = raster::extract(pop, sp::SpatialPoints(lonLatGrid),method="bilinear")
  
  # determine which points are urban
  newPop = data.frame(list(lon=lonLatGrid[,1], lat=lonLatGrid[,2], pop=interpPopVals, area=areas, subareas=subareas))
  if(is.null(poppsub)) {
    # the set thresholds at the area level
    threshes = setThresholdsArea(newPop, poppa)
    popThreshes = sapply(1:nrow(newPop), function(i) {threshes$threshes[threshes$areas == newPop$area[i]]})
    urban = newPop$pop >= unlist(popThreshes)
  } else {
    # the set thresholds at the subarea level
    threshes = setThresholdsSubarea(newPop, poppsub)
    popThreshes = sapply(1:nrow(newPop), function(i) {threshes$threshes[threshes$subareas == newPop$subarea[i]]})
    urban = newPop$pop >= unlist(popThreshes)
  }
  
  newPop$urban = urban
  
  newPop$east = pixelGrid$east
  newPop$north = pixelGrid$north
  
  # if necessary, renormalize population values within subareas crossed with 
  # urban/rural to be the correct value
  if(!is.null(poppsub)) {
    for(i in 1:nrow(poppsub)) {
      thisSub = poppsub$subarea[i]
      substratumUrban = (subareas == thisSub) & urban
      substratumRural = (subareas == thisSub) & !urban
      factorUrban = poppsub$popUrb[i] / sum(newPop$pop[substratumUrban])
      factorRural = poppsub$popRur[i] / sum(newPop$pop[substratumRural])
      newPop$pop[substratumUrban] = newPop$pop[substratumUrban] * factorUrban
      newPop$pop[substratumRural] = newPop$pop[substratumRural] * factorRural
    }
  } else {
    # same as above but at the area level
    for(i in 1:nrow(poppa)) {
      thisArea = poppa$area[i]
      substratumUrban = (areas == thisArea) & urban
      substratumRural = (areas == thisArea) & !urban
      factorUrban = poppa$popUrb[i] / sum(newPop$pop[substratumUrban])
      factorRural = poppa$popRur[i] / sum(newPop$pop[substratumRural])
      newPop$pop[substratumUrban] = newPop$pop[substratumUrban] * factorUrban
      newPop$pop[substratumRural] = newPop$pop[substratumRural] * factorRural
    }
  }
  
  newPop
}

#' @describeIn generatePixelGrid adjust population densities in grid based on population frame
#' @export
adjustPopMat = function(popMat, poppaTarget=NULL, adjustBy=c("area", "subarea")) {
  adjustBy = match.arg(adjustBy)
  
  # sort get population per stratum from poppaTarget
  if(adjustBy == "area") {
    areas=sort(unique(poppaTarget$area))
  } else {
    areas=sort(unique(poppaTarget$subarea))
  }
  
  targetPopPerStratumUrban = poppaTarget$popUrb
  targetPopPerStratumRural = poppaTarget$popRur
  
  # generate 2 nArea x nPixels matrices for urban and rural strata integrating pixels with respect to population density to get county estimates
  getAreaStratumIntegrationMatrix = function(getUrban=TRUE) {
    areas = as.character(areas)
    
    if(adjustBy == "area") {
      mat = t(sapply(areas, function(area) {
        popMat$area == area
      }))
    } else {
      mat = t(sapply(areas, function(area) {
        popMat$subarea == area
      }))
    }
    
    mat = sweep(mat, 2, popMat$pop, "*")
    
    sweep(mat, 2, popMat$urban == getUrban, "*")
  }
  urbanIntegrationMat = getAreaStratumIntegrationMatrix()
  ruralIntegrationMat = getAreaStratumIntegrationMatrix(FALSE)
  
  # calculate number of people per stratum by integrating the population density surface
  urbanPopulations = rowSums(urbanIntegrationMat)
  ruralPopulations = rowSums(ruralIntegrationMat)
  
  # adjust each row of the integration matrices to get the correct expected number of children per stratum
  urbanIntegrationMat = sweep(urbanIntegrationMat, 1, targetPopPerStratumUrban / urbanPopulations, "*")
  urbanIntegrationMat[urbanPopulations == 0,] = 0
  ruralIntegrationMat = sweep(ruralIntegrationMat, 1, targetPopPerStratumRural / ruralPopulations, "*")
  ruralIntegrationMat[ruralPopulations == 0,] = 0
  
  # the column sums of the matrices give the correct modified population densities
  popMat$pop = colSums(urbanIntegrationMat) + colSums(ruralIntegrationMat)
  
  popMat
}

#' @describeIn generatePixelGrid Generate a population frame of a similar format to poppa argument of \code{\link{simPopPixel}} with a custom set of regions
#' @export
poppRegionFromPopMat = function(popMat, regions) {
  out = aggregate(popMat$pop, by=list(region=as.character(regions), urban=popMat$urban), FUN=sum, drop=FALSE)
  regions = sort(unique(out$region))
  poppr = data.frame(region=regions, popUrb=out[(length(regions) + 1):(2*length(regions)), 3], 
                     popRur=out[1:length(regions), 3])
  poppr$popUrb[is.na(poppr$popUrb)] = 0
  poppr$popRur[is.na(poppr$popRur)] = 0
  poppr$popTotal = poppr$popUrb + poppr$popRur
  poppr
}