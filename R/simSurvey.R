# functions for simulating surveys


#' Simulate DHS-like surveys from EA- or HH-level population information
#' 
#' Given a spatial risk model, simulate populations and population prevalences at the 
#' enumeration area level (represented as points), and aggregate to the pixel and 
#' administrative areal level.
#' 
#' @param popSim Simulated population. Must contain EA or HH level population information, i.e. 
#'               either an element named 'eaPop' or 'eaSimDat' containing a list of 
#'               populations simulated at the EA level. See \code{\link{simPopCustom}}
#' @param survDat simulated survey from \code{sampleClusterSurveys}
#' @param stratOrder if provided, will sort resulting strata to be in the same order as in 
#             stratOrder after simulating the survey. Otherwise, sorted alphabetically
#' @param stratName The name of the survey stratum variable to aggregate the number 
#'                  of clusters over. Defaults to 'area'.
#' @param HHperClust Fixed number of households per selected EA to include in the survey
#' @param n Number of simulated surveys. NOTE: storing household level information can 
#' get quite memory intensive. Recommended to simulate 1 household level population and 1 
#' associated survey at a time.
#' @param eaSampleStrat EA level sampling strategy. Either 'pps' or 'srs'. If 'pps', samples 
#' EAs with probability proportional to the number of households in the EA. Within each sampled 
#' EA, households are always sampled via srs.
#' @param clustpaList a list of data.frame objects, 1 for each survey to sample, providing 
#' information on the number of clusters to sample per stratum. Should contain elements:
#' \describe{
#'   \item{area}{name of area}
#'   \item{EAUrb}{number of urban enumeration areas in the area}
#'   \item{EARur}{number of rural enumeration areas in the area}
#'   \item{EATotal}{total number of enumeration areas in the area}
#'   \item{HHUrb}{number of urban households in the area}
#'   \item{HHRur}{number of rural households in the area}
#'   \item{HHTotal}{total number of households in the area}
#' }
#' @param fixPopPerHH Currently only 1 or NULL is supported. If 1, fixes the target population 
#'                    to be 1 in each simulated household (requires EA populations and 
#'                    num households are equal). If NULL, randomly distributes population 
#'                    among the households (default)
#' @param verbose If TRUE, prints progress and all warnings
#' @param seed If not NULL, the random number seed to set at the beginning of the function
#' @return The simulated survey data, household level population, or aggregated survey information from 
#' \code{sampleClusterSurveys}, \code{getHHpop}, and \code{getClustpaFromSurvey} respectively.
#' 
#' @details Functions for simulating multilevel household cluster surveys based on simulated populations 
#' from other SUMMER functions. For surveys with structures similar to the DHS or MICS. Information on 
#' spatial coordinates included. NOTE: storing household level information can 
#' get quite memory intensive. Recommended to simulate 1 household level population and 1 
#' associated survey at a time. Multiple surveys can be sampled from the same household level population
#' 
#' @author John Paige
#' @seealso \code{\link{simPopCustom}}, \code{\link{simPopSPDE}}, \code{\link{makePopIntegrationTab}}, \code{\link{adjustPopMat}}, \code{\link{simSPDE}}.
#' @examples 
#' \dontrun{
#' ## In this script we will create 5km resolution pixellated grid over Kenya, 
#' ## and generate tables of estimated (both target and general) population 
#' ## totals at the area (e.g. Admin-1) and subarea (e.g. Admin-2) levels. Then 
#' ## we will use that to simulate populations and associated surveys of these
#' 
#' # download Kenya GADM shapefiles from SUMMERdata github repository
#' githubURL <- paste0("https://github.com/paigejo/SUMMERdata/blob/main/data/", 
#'                     "kenyaMaps.rda?raw=true")
#' tempDirectory = "~/"
#' mapsFilename = paste0(tempDirectory, "/kenyaMaps.rda")
#' if(!file.exists(mapsFilename)) {
#'   download.file(githubURL,mapsFilename)
#' }
#' 
#' # load it in
#' # out = load(mapsFilename)
#' out = load(url(githubURL))
#' out
#' adm1@data$NAME_1 = as.character(adm1@data$NAME_1)
#' adm1@data$NAME_1[adm1@data$NAME_1 == "Trans Nzoia"] = "Trans-Nzoia"
#' adm1@data$NAME_1[adm1@data$NAME_1 == "Elgeyo-Marakwet"] = "Elgeyo Marakwet"
#' adm2@data$NAME_1 = as.character(adm2@data$NAME_1)
#' adm2@data$NAME_1[adm2@data$NAME_1 == "Trans Nzoia"] = "Trans-Nzoia"
#' adm2@data$NAME_1[adm2@data$NAME_1 == "Elgeyo-Marakwet"] = "Elgeyo Marakwet"
#' 
#' # some Admin-2 areas have the same name
#' adm2@data$NAME_2 = as.character(adm2@data$NAME_2)
#' adm2@data$NAME_2[(adm2@data$NAME_1 == "Bungoma") & 
#'                    (adm2@data$NAME_2 == "Lugari")] = "Lugari, Bungoma"
#' adm2@data$NAME_2[(adm2@data$NAME_1 == "Kakamega") & 
#'                    (adm2@data$NAME_2 == "Lugari")] = "Lugari, Kakamega"
#' adm2@data$NAME_2[(adm2@data$NAME_1 == "Meru") & 
#'                    (adm2@data$NAME_2 == "Igembe South")] = "Igembe South, Meru"
#' adm2@data$NAME_2[(adm2@data$NAME_1 == "Tharaka-Nithi") & 
#'                    (adm2@data$NAME_2 == "Igembe South")] = "Igembe South, Tharaka-Nithi"
#' 
#' # The spatial area of unknown 8 is so small, it causes problems unless its removed or 
#' # unioned with another subarea. Union it with neighboring Kakeguria:
#' newadm2 = adm2
#' unknown8I = which(newadm2$NAME_2 == "unknown 8")
#' newadm2$NAME_2[newadm2$NAME_2 %in% c("unknown 8", "Kapenguria")] <- 
#'   "Kapenguria + unknown 8"
#' admin2.IDs <- newadm2$NAME_2
#' 
#' newadm2@data = cbind(newadm2@data, NAME_2OLD = newadm2@data$NAME_2)
#' newadm2@data$NAME_2OLD = newadm2@data$NAME_2
#' newadm2@data$NAME_2 = admin2.IDs
#' newadm2$NAME_2 = admin2.IDs
#' temp <- terra::aggregate(as(newadm2, "SpatVector"), by="NAME_2")
#' 
#' library(sf)
#' temp <- sf::st_as_sf(temp)
#' temp <- sf::as_Spatial(temp)
#' 
#' tempData = newadm2@data[-unknown8I,]
#' tempData = tempData[order(tempData$NAME_2),]
#' newadm2 <- SpatialPolygonsDataFrame(temp, tempData, match.ID = F)
#' adm2 = newadm2
#' 
#' # download 2014 Kenya population density TIF file
#' 
#' githubURL <- paste0("https://github.com/paigejo/SUMMERdata/blob/main/data/", 
#'                     "Kenya2014Pop/worldpop_total_1y_2014_00_00.tif?raw=true")
#' popTIFFilename = paste0(tempDirectory, "/worldpop_total_1y_2014_00_00.tif")
#' if(!file.exists(popTIFFilename)) {
#'   download.file(githubURL,popTIFFilename)
#' }
#' 
#' # load it in
#' pop = terra::rast(popTIFFilename)
#' 
#' ver = terra::gdal(lib="proj")
#' PROJ6 <- as.numeric(substr(ver, 1, 1)) >= 6
#' 
#' # from lon/lat coords to easting/northing
#' if(!PROJ6) {
#'   crs(pop) = "+proj=longlat"
#' } else {
#'   crs(pop) = "EPSG:4326"
#' }
#' 
#' eastLim = c(-110.6405, 832.4544)
#' northLim = c(-555.1739, 608.7130)
#' 
#' ## Construct poppsubKenya, a table of urban/rural general population totals 
#' ## in each subarea. Technically, this is not necessary since we can load in 
#' ## poppsubKenya via data(kenyaPopulationData). First, we will need to calculate 
#' ## the areas in km^2 of the areas and subareas
#' 
#' # use Lambert equal area projection of areas (Admin-1) and subareas (Admin-2)
#' midLon = mean(adm1@bbox[1,])
#' midLat = mean(adm1@bbox[2,])
#' p4s = paste0("+proj=laea +x_0=0 +y_0=0 +lon_0=", midLon, 
#'              " +lat_0=", midLat, " +units=km")
#' 
#' adm1_sf = st_as_sf(adm1)
#' adm1proj_sf = st_transform(adm1_sf, p4s)
#' adm1proj = as(adm1proj_sf, "Spatial")
#' 
#' adm2_sf = st_as_sf(adm2)
#' adm2proj_sf = st_transform(adm2_sf, p4s)
#' adm2proj = as(adm2proj_sf, "Spatial")
#' 
#' # now calculate spatial area in km^2
#' admin1Areas = as.numeric(st_area(adm1proj_sf))
#' admin2Areas = as.numeric(st_area(adm2proj_sf))
#' 
#' areapaKenya = data.frame(area=adm1proj@data$NAME_1, spatialArea=admin1Areas)
#' areapsubKenya = data.frame(area=adm2proj@data$NAME_1, subarea=adm2proj@data$NAME_2, 
#'                            spatialArea=admin2Areas)
#' 
#' # Calculate general population totals at the subarea (Admin-2) x urban/rural 
#' # level and using 1km resolution population grid. Assign urbanicity by 
#' # thresholding population density based on estimated proportion population 
#' # urban/rural, making sure total area (Admin-1) urban/rural populations in 
#' # each area matches poppaKenya.
#' 
#' # NOTE: the following function will typically take ~15-20 minutes. Can speed up 
#' #       by setting kmRes to be higher, but we recommend fine resolution for 
#' #       this step, since it only needs to be done once.
#' system.time(poppsubKenya <- getPoppsub(
#'   kmRes=1, pop=pop, domainMapDat=adm0,
#'   eastLim=eastLim, northLim=northLim, mapProjection=projKenya,
#'   poppa = poppaKenya, areapa=areapaKenya, areapsub=areapsubKenya, 
#'   areaMapDat=adm1, subareaMapDat=adm2, 
#'   areaNameVar = "NAME_1", subareaNameVar="NAME_2"))
#' 
#' # Now generate a general population integration table at 5km resolution, 
#' # based on subarea (Admin-2) x urban/rural population totals. This takes 
#' # ~1 minute
#' system.time(popMatKenya <- makePopIntegrationTab(
#'   kmRes=5, pop=pop, domainMapDat=adm0,
#'   eastLim=eastLim, northLim=northLim, mapProjection=projKenya,
#'   poppa = poppaKenya, poppsub=poppsubKenya, 
#'   areaMapDat = adm1, subareaMapDat = adm2,
#'   areaNameVar = "NAME_1", subareaNameVar="NAME_2"))
#' 
#' ## Adjust popMat to be target (neonatal) rather than general population 
#' ## density. First create the target population frame
#' ## (these numbers are based on IPUMS microcensus data)
#' mothersPerHouseholdUrb = 0.3497151
#' childrenPerMotherUrb = 1.295917
#' mothersPerHouseholdRur = 0.4787696
#' childrenPerMotherRur = 1.455222
#' targetPopPerStratumUrban = easpaKenya$HHUrb * mothersPerHouseholdUrb * 
#'   childrenPerMotherUrb
#' targetPopPerStratumRural = easpaKenya$HHRur * mothersPerHouseholdRur * 
#'   childrenPerMotherRur
#' easpaKenyaNeonatal = easpaKenya
#' easpaKenyaNeonatal$popUrb = targetPopPerStratumUrban
#' easpaKenyaNeonatal$popRur = targetPopPerStratumRural
#' easpaKenyaNeonatal$popTotal = easpaKenyaNeonatal$popUrb + 
#'   easpaKenyaNeonatal$popRur
#' easpaKenyaNeonatal$pctUrb = 100 * easpaKenyaNeonatal$popUrb / 
#'   easpaKenyaNeonatal$popTotal
#' easpaKenyaNeonatal$pctTotal = 
#'   100 * easpaKenyaNeonatal$popTotal / sum(easpaKenyaNeonatal$popTotal)
#' 
#' # Generate the target population density by scaling the current 
#' # population density grid at the Admin1 x urban/rural level
#' popMatKenyaNeonatal = adjustPopMat(popMatKenya, easpaKenyaNeonatal)
#' 
#' # Generate neonatal population table from the neonatal population integration 
#' # matrix. This is technically not necessary for population simulation purposes, 
#' # but is here for illustrative purposes
#' poppsubKenyaNeonatal = poppRegionFromPopMat(popMatKenyaNeonatal, 
#'                                             popMatKenyaNeonatal$subarea)
#' poppsubKenyaNeonatal = 
#'   cbind(subarea=poppsubKenyaNeonatal$region, 
#'         area=adm2@data$NAME_1[match(poppsubKenyaNeonatal$region, adm2@data$NAME_2)], 
#'         poppsubKenyaNeonatal[,-1])
#' print(head(poppsubKenyaNeonatal))
#' 
#' ## Now we're ready to simulate neonatal populations along with neonatal 
#' ## mortality risks and prevalences
#' 
#' # use the following model to simulate the neonatal population based roughly 
#' # on Paige et al. (2020) neonatal mortality modeling for Kenya.
#' beta0=-2.9 # intercept
#' gamma=-1 # urban effect
#' rho=(1/3)^2 # spatial variance
#' effRange = 400 # effective spatial range in km
#' sigmaEpsilon=sqrt(1/2.5) # cluster (nugget) effect standard deviation
#' 
#' # Run a simulation! This produces multiple dense nEA x nsim and nPixel x nsim 
#' # matrices. In the future sparse matrices and chunk by chunk computations 
#' # may be incorporated.
#' simPop = simPopSPDE(nsim=1, easpa=easpaKenyaNeonatal, 
#'                     popMat=popMatKenya, targetPopMat=popMatKenyaNeonatal, 
#'                     poppsub=poppsubKenya, spdeMesh=kenyaMesh, 
#'                     margVar=rho, sigmaEpsilon=sigmaEpsilon, 
#'                     gamma=gamma, effRange=effRange, beta0=beta0, 
#'                     seed=12, inla.seed=12, nHHSampled=25, 
#'                     stratifyByUrban=TRUE, subareaLevel=TRUE, 
#'                     doFineScaleRisk=TRUE, doSmoothRisk=TRUE, 
#'                     min1PerSubarea=TRUE)
#' 
#' # get average absolute percent error relative to fine scale prevalence at Admin-2 level
#' tempDat = simPop$subareaPop$aggregationResults[c("region", "pFineScalePrevalence", 
#'                                                   "pFineScaleRisk", "pSmoothRisk")]
#' 100*mean(abs(tempDat$pFineScalePrevalence - tempDat$pFineScaleRisk) / 
#'            tempDat$pFineScalePrevalence)
#' 100*mean(abs(tempDat$pFineScalePrevalence - tempDat$pSmoothRisk) / 
#'            tempDat$pFineScalePrevalence)
#' 100*mean(abs(tempDat$pFineScaleRisk - tempDat$pSmoothRisk) / 
#'            tempDat$pFineScalePrevalence)
#' 
#' # verify number of EAs per area and subarea
#' cbind(aggregate(simPop$eaPop$eaSamples[,1], by=list(area=popMatKenya$area), FUN=sum), 
#'       trueNumEAs=easpaKenya$EATotal[order(easpaKenya$area)])
#' aggregate(simPop$eaPop$eaSamples[,1], by=list(area=popMatKenya$subarea), FUN=sum)
#' 
#' ## generate a survey assuming 1 target population member per household from the 
#' ## simulated population
#' thisEApop = simPop$eaPop$eaDatList[1]
#' 
#' ## get associated HH level population information
#' thisHHpop = getHHpop(thisEApop, fixPopPerHH=fixPopPerHH)
#' 
#' ## sample DHS survey for this population. Must define clustpaDHS data.frame 
#' ## first based on your desired survey design
#' survDHS = sampleClusterSurveys(1, thisHHpop, HHperClust=25, clustpaList=list(clustpaDHS))
#' }
#' @name simSurvey
NULL

#' @importFrom stats rmultinom
#' @describeIn simSurvey
#' Simulates household level population data given EA level population data
#' 
#' @export
getHHpop = function(popSim, fixPopPerHH=NULL, verbose=TRUE) {
  # first get ea level data from popSim
  if("eaPop" %in% names(popSim)) {
    eaPop = popSim$eaPop
    eaPopDat = eaPop$eaDatList
  } else if("eaDatList" %in% names(popSim)) {
    eaPopDat = popSim$eaDatList
    
  } else if(is.list(popSim)) {
    if("nHH" %in% names(popSim[[1]])) {
      eaPopDat = popSim
    } else {
      stop("popSim has no EA level information. Could set 'returnEAinfo' to TRUE in simPopCustom")
    }
  } else {
    stop("popSim has no EA level information. Could set 'returnEAinfo' to TRUE in simPopCustom")
  }
  
  if((length(eaPopDat) > 1) && verbose) {
    warning("length(eaPopDat) > 1, so there may be a lot of household level data, and memory issues accordingly...")
  }
  
  HHdat = list()
  for(i in 1:length(eaPopDat)) {
    eaDat = eaPopDat[[i]]
    
    # first get the number of households
    numHouseholds = eaDat$nHH
    
    # now expand the eaDat table to be in the long format, where each row is a house
    rowsLong = rep(1:nrow(eaDat), numHouseholds)
    thisHHdat = eaDat[rowsLong, ]
    thisHHdat$eaIs = rowsLong
    urbanHH = thisHHdat$urban
    
    # function for randomly spreading people among households in long form data:
    extendThisDat = function(xs, nHH) {
      revenMultinom = function(sizeK) {
        size = sizeK[1]
        k = sizeK[2]
        prob = rep(1/k, k)
        stats::rmultinom(1, size, prob)
      }
      unlist(apply(cbind(xs, nHH), 1, revenMultinom))
    }
    
    # generate how many of the target population are in each cluster
    if(is.null(fixPopPerHH)) {
      lived = extendThisDat(eaDat$N - eaDat$Z, numHouseholds)
      died = extendThisDat(eaDat$Z, numHouseholds)
    } else if(fixPopPerHH == 1) {
      extendDatEven = function(Ns, Zs) {
        
        spreadAmongHHs = function(thisRow) {
          thisN = thisRow[1]
          thisZ = thisRow[2]
          
          # spread population evenly among households
          # nHH = thisRow$nHH
          # hhI = sample(1:nHH, nHH, replace=FALSE)
          c(rep(1, thisZ), rep(0, thisN - thisZ))
        }
        c(unlist(apply(cbind(Ns, Zs), 1, spreadAmongHHs)))
      }
      died = extendDatEven(eaDat$N, eaDat$Z)
      lived = 1 - died
    } else {
      stop("If fixPopPerHH is not NULL it must be 1. Other values not currently supported")
    }
    
    thisHHdat$N = died + lived
    thisHHdat$Z = died
    thisHHdat$nHH = 1
    thisHHdat$pFineScalePrevalence = thisHHdat$Z / thisHHdat$N
    thisHHdat$pFineScalePrevalence[thisHHdat$N == 0] = 0 # NaN otherwise
    
    HHdat = c(HHdat, list(thisHHdat))
  }
  
  HHdat
}

#' @importFrom stats aggregate
#' @importFrom sampling UPmidzuno
#' @describeIn simSurvey
#' Simulates household cluster surveys from input population with EA or HH level info
#' 
#' @export
sampleClusterSurveys = function(n=NULL, popSim=NULL, HHperClust=25, fixPopPerHH=NULL, 
                                eaSampleStrat=c("pps", "srs"), clustpaList, 
                                seed=NULL, verbose=FALSE) {
  
  eaSampleStrat = match.arg(eaSampleStrat)
  
  # set random seed if supplied
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # initialize EA and HH level population info to NULL
  hhPopDat = NULL
  eaPopDat = NULL
  
  # first get ea level data from popSim if need be
  if("eaPop" %in% names(popSim)) {
    eaPop = popSim$eaPop
    eaPopDat = eaPop$eaDatList
  } else if("eaDatList" %in% names(popSim)) {
    eaPopDat = popSim$eaDatList
    
  } else if(is.list(popSim)) {
    if("nHH" %in% names(popSim[[1]])) {
      if(all(popSim[[1]]$nHH == 1)) {
        hhPopDat = popSim
      } else {
        eaPopDat = popSim
      }
    }
  } else {
    stop("popSim has no EA or HH level information. Could set 'returnEAinfo' to TRUE in simPopCustom")
  }
  
  
  # set n if unset
  if(is.null(n)) {
    if(!is.null(hhPopDat)) {
      n = length(hhPopDat)
    } else {
      n = length(eaPopDat)
    }
  }
  
  if((length(clustpaList) > 1) && (length(clustpaList) != n)) {
    stop(paste0("length mismatch between n (", n, ") and length(clustpaList) (", length(clustpaList), ")"))
  }
  
  # if popDat has length 1 and n > 1, sample multiple surveys from the same population
  if(!is.null(eaPopDat)) {
    resamplePop = ifelse((length(eaPopDat) == 1) && (n > 1), TRUE, FALSE)
  } else if(!is.null(hhPopDat)) {
    resamplePop = ifelse((length(hhPopDat) == 1) && (n > 1), TRUE, FALSE)
  }
  
  if(verbose) {
    print("Simulating surveys")
  }
  surveys = list()
  for(i in 1:n) {
    if(verbose) {
      print(paste0("Sampling survey ", i, "/", n))
    }
    
    # get number of clusters to sample per area
    thisClustpa = clustpaList[[min(c(length(clustpaList), i))]]
    
    # get HH level population information
    if(!is.null(eaPopDat)) {
      # if no HH info available, first get the EA level info, then use to get HH level info
      
      # either resample surveys from 1 population, or sample 1 survey per population
      if(resamplePop) {
        eaDat = eaPopDat[[1]]
      } else {
        eaDat = eaPopDat[[i]]
      }
      
      hhDat = getHHpop(list(eaDatList = list(eaDat)), fixPopPerHH = fixPopPerHH)[[1]]
      
    } else if(!is.null(hhPopDat)) {
      # we already have HH level info
      
      # either resample surveys from 1 population, or sample 1 survey per population
      if(resamplePop) {
        hhDat = hhPopDat[[1]]
      } else {
        hhDat = hhPopDat[[i]]
      }
      eaDat = NULL
    }
    
    # obtain info about the EAs
    uniqueEAIs = sort(unique(hhDat$eaIs))
    eaUrbs = hhDat$urban[match(uniqueEAIs, hhDat$eaIs)]
    eaAreas = hhDat$area[match(uniqueEAIs, hhDat$eaIs)]
    
    if(eaSampleStrat == "pps") {
      # calculate number of HHs per EA
      aggHHs = stats::aggregate(hhDat$nHH, by=list(hhDat$eaIs), FUN=sum)
      eaHHs = aggHHs$x
    }
    
    # sample EAs
    sampledEAIs = numeric(0)
    inclusionProbs = numeric(0)
    for(j in 1:nrow(thisClustpa)) {
      
      thisArea = thisClustpa$area[j]
      nUrbEA = thisClustpa$EAUrb[j]
      nRurEA = thisClustpa$EARur[j]
      thisEAIs = uniqueEAIs[eaAreas == thisArea]
      thisEAUrbs = eaUrbs[eaAreas == thisArea]
      thisEAIsUrb = thisEAIs[thisEAUrbs]
      thisEAIsRur = thisEAIs[!thisEAUrbs]
      
      if(eaSampleStrat == "pps") {
        thisEAhhs = eaHHs[eaAreas == thisArea]
        thisEAhhsUrb = thisEAhhs[thisEAUrbs]
        thisEAhhsRur = thisEAhhs[!thisEAUrbs]
      }
      
      # sample EAs and HHs for this area
      sampUrbEAIs = numeric(0)
      sampRurEAIs = numeric(0)
      inclusionProbsUrb = numeric(0)
      inclusionProbsRur = numeric(0)
      if(eaSampleStrat == "srs") {
        if(nUrbEA != 0) {
          sampUrbEAIs = sample(thisEAIsUrb, nUrbEA, replace=F)
          inclusionProbsUrb = rep(nUrbEA/length(thisEAIsUrb), nUrbEA)
        }
        if(nRurEA != 0) {
          sampRurEAIs = sample(thisEAIsRur, nRurEA, replace=F)
          inclusionProbsRur = rep(nRurEA/length(thisEAIsRur), nRurEA)
        }
      } else if(eaSampleStrat == "pps") {
        require(sampling)
        if(nUrbEA != 0) {
          # sampUrbEAIs = sample(thisEAIsUrb, nUrbEA, replace=F, prob=thisEAhhsUrb/sum(thisEAhhsUrb))
          inclusionProbsUrb = nUrbEA * thisEAhhsUrb/sum(thisEAhhsUrb)
          sampUrbEAIs = thisEAIsUrb[as.logical(sampling::UPmidzuno(inclusionProbsUrb))]
        }
        if(nRurEA != 0) {
          # sampRurEAIs = sample(thisEAIsRur, nRurEA, replace=F, prob=thisEAhhsRur/sum(thisEAhhsRur))
          inclusionProbsRur = nRurEA * thisEAhhsRur/sum(thisEAhhsRur)
          sampRurEAIs = thisEAIsRur[as.logical(sampling::UPmidzuno(inclusionProbsRur))]
        }
      } else {
        stop(paste0("eaSampleStrat '", eaSampleStrat, "' not supported"))
      }
      
      # concatenate urban and rural EAs sampled from this area
      thisSampEAIs = c(sampUrbEAIs, sampRurEAIs)
      thisInclusionProbs = c(inclusionProbsUrb, inclusionProbsRur)
      
      # concatenate to vector of all EAs sampled from all areas
      sampledEAIs = c(sampledEAIs, thisSampEAIs)
      inclusionProbs = c(inclusionProbs, thisInclusionProbs)
    }
    
    # subset hhDat to only EAs sampled
    hhSubdat = hhDat[hhDat$eaIs %in% sampledEAIs,]
    
    # sample HHs within chosen EAs
    hhIsTab = stats::aggregate(1:nrow(hhSubdat), by=list(eaIs=hhSubdat$eaIs), function(x) {
      sample(x, HHperClust, replace=FALSE)
    })
    hhIs = sort(c(as.matrix(hhIsTab[,-1])))
    
    hhDatSample = hhSubdat[hhIs,]
    
    # aggregate HH level data for only the EAs sampled
    aggTab = lapply(1:ncol(hhSubdat), function(j) {
      varName = names(hhSubdat)[j]
      stats::aggregate(hhSubdat[,j], by=list(hhSubdat$eaIs), FUN = function(x) {
        if(varName %in% c("N", "nHH", "Z")) {
          sum(x, na.rm=TRUE)
        } else if(is.numeric(x)) {
          mean(x, na.rm=TRUE)
        } else {
          x[1]
        }
      })$x
    })
    names(aggTab) = names(hhSubdat)
    aggTab = as.data.frame(aggTab)
    aggTab$pFineScalePrevalence = aggTab$Z/aggTab$N
    aggTab$pFineScalePrevalence[aggTab$N == 0] = 0
    aggTab$includeProbEA = inclusionProbs[match(sampledEAIs, aggTab$eaIs)]
    
    eaDatSampled = aggTab
    
    # do the same for the actual sampled HHs
    aggTab = lapply(1:ncol(hhDatSample), function(j) {
      varName = names(hhDatSample)[j]
      aggregate(hhDatSample[,j], by=list(hhDatSample$eaIs), FUN = function(x) {
        if(varName %in% c("N", "nHH", "Z")) {
          sum(x, na.rm=TRUE)
        } else if(is.numeric(x)) {
          mean(x, na.rm=TRUE)
        } else {
          x[1]
        }
      })$x
    })
    names(aggTab) = names(hhDatSample)
    aggTab = as.data.frame(aggTab)
    aggTab$pFineScalePrevalence = aggTab$Z/aggTab$N
    aggTab$pFineScalePrevalence[aggTab$N == 0] = 0
    aggTab$includeProbEA = inclusionProbs[match(sampledEAIs, aggTab$eaIs)]
    
    surveyDat = aggTab
    
    # calculate final sampling weight and number HHs in the full EA
    surveyDat$nHHsEA = eaDatSampled$nHH[match(surveyDat$eaIs, eaDatSampled$eaIs)]
    surveyDat$includeProbHH = surveyDat$nHH / surveyDat$nHHsEA
    surveyDat$samplingWeight = surveyDat$N / (surveyDat$includeProbHH * surveyDat$includeProbEA)
    
    
    # concatenate results
    surveys = c(surveys, list(surveyDat))
  }
  
  surveys
}

#' @importFrom stats aggregate
#' @describeIn simSurvey
#' Collects information about number of clusters per stratum into a data.frame
#' 
#' @export
getClustpaFromSurvey = function(survDat, stratName="area", stratOrder=sort(survDat[[stratName]]), HHperClust=25) {
  
  if("ns" %in% names(survDat)) {
    survDat$n = survDat$ns
  }
  
  nAreas = length(unique(survDat[[stratName]]))
  
  # first the number of sampled clusters
  nEAtabTmp = stats::aggregate(survDat$n, by=list(strat=survDat[[stratName]], urban=survDat$urban), FUN=length, drop=FALSE)
  nEAtabTmp[is.na(nEAtabTmp[,3]), 3] = 0
  
  nEAtab = data.frame(nEAtabTmp[1:nAreas, 1], EAUrb=nEAtabTmp[(nAreas+1):(2*nAreas), 3], EARur=nEAtabTmp[1:nAreas, 3])
  names(nEAtab)[1] = stratName
  nEAtab$EATotal = nEAtab$EAUrb + nEAtab$EARur
  
  # initialize clustpa
  clustpa = nEAtab
  
  # second the number of households
  clustpa$HHUrb = clustpa$EAUrb * HHperClust
  clustpa$HHRur = clustpa$EARur * HHperClust
  clustpa$HHTotal = clustpa$EATotal * HHperClust
  
  # third the number of people
  popTabTmp = stats::aggregate(survDat$n, by=list(strat=survDat[[stratName]], urban=survDat$urban), FUN=sum, drop=FALSE)
  popTabTmp[is.na(popTabTmp[,3]), 3] = 0
  # urbanToRuralI = c(1:27, 29, 31:47) # skip mombasa and nairobi
  popTab = data.frame(popTabTmp[1:nAreas, 1], popUrb=popTabTmp[(nAreas+1):(2*nAreas), 3], popRur=popTabTmp[1:nAreas, 3])
  names(popTab)[1] = stratName
  popTab$popTotal = popTab$popUrb + popTab$popRur
  
  # concatenate cluster level denominator info to clustpa info
  clustpa$popUrb = popTab$popUrb
  clustpa$popRur = popTab$popRur
  clustpa$popTotal = popTab$popTotal
  
  # sort if need be
  if(!is.null(stratOrder)) {
    ordering = order(stratOrder)
    clustpa = clustpa[ordering,]
  }
  
  clustpa
}