# Functions for simulating populations

#' Simulate populations and areal prevalences
#' 
#' Given a spatial risk model, simulate populations and population prevalences at the 
#' enumeration area level (represented as points), and aggregate to the pixel and 
#' administrative areal level.
#' 
#' @param nsim number of simulations
#' @param easpa data.frame of enumeration area, households, and target population per area stratified by urban/rural with variables:
#' \describe{
#'   \item{area}{name of area}
#'   \item{EAUrb}{number of urban enumeration areas in the area}
#'   \item{EARur}{number of rural enumeration areas in the area}
#'   \item{EATotal}{total number of enumeration areas in the area}
#'   \item{HHUrb}{number of urban households in the area}
#'   \item{HHRur}{number of rural households in the area}
#'   \item{HHTotal}{total number of households in the area}
#'   \item{popUrb}{total urban (target) population of area}
#'   \item{popRur}{total rural (target) population of area}
#'   \item{popTotal}{total (general) population of area}
#' }
#' @param popMat pixellated grid data frame with variables `lon`, `lat`, `pop`, `area`, `subareas` (if subareaLevel is TRUE), `urban` (if stratifyByUrban is TRUE), `east`, and `north`
#' @param targetPopMat same as popMat, but `pop` variable gives target rather than general population
#' @param poppsub data.frame of population per subarea separated by 
#' urban/rural using for population density grid normalization or urbanicity 
#' classification. Often based on extra fine scale population density grid. Has variables:
#' @param spdeMesh triangular mesh for the SPDE
#' @param margVar marginal variance of the spatial process, excluding cluster effects. 
#'          If 0, no spatial component is included
#' @param effRange effective spatial range for the SPDE model
#' @param beta0 intercept of logit model for risk
#' @param gamma effect of urban on logit scale for logit model for risk
#' @param sigmaEpsilon standard deviation on the logit scale for iid Gaussian EA level random effects in the risk model
#' @param seed random number generator seed
#' @param inla.seed seed input to inla.qsample. 0L sets seed intelligently, 
#'            > 0 sets a specific seed, < 0 keeps existing RNG
#' @param nHHSampled number of households sampled per enumeration area. Default is 25 to match DHS surveys
#' @param stratifyByUrban whether or not to stratify simulations by urban/rural classification
#' @param subareaLevel whether or not to aggregate the population by subarea
#' @param doFineScaleRisk whether or not to calculate the fine scale risk at each aggregation level in addition to the prevalence
#' @param min1PerSubarea if TRUE, ensures there is at least 1 EA per subarea. If subareas are particularly unlikely to 
#' have enumeration areas since they have a very low proportion of the population in an area, then setting this to TRUE may be 
#' computationally intensive.
#' @param uDraws nPixels x nsim dimension matrix of draws from the spatial risk field on logit scale
#' @param sigmaEpsilonDraws nsim length vector of draws of cluster effect logit scale SD (joint draws with uDraws)
#' @param validationPixelI CURRENTLY FOR TESTING PURPOSES ONLY a set of indices of pixels for which we want to simulate populations (used for pixel level validation)
#' @param validationClusterI CURRENTLY FOR TESTING PURPOSES ONLY a set of indices of cluster for which we want to simulate populations (used for cluster level validation)
#' @param clustersPerPixel CURRENTLY FOR TESTING PURPOSES ONLY Used for pixel level validation. Fixes the number of EAs per pixel.
#' @param returnEAinfo if TRUE, returns information on every individual EA (BAU) for each simulated population
#' @param epsc nEAs x nsim matrix of simulated EA (BAU) level iid effects representing fine scale variation in 
#'       risk. If NULL, they are simulated as iid Gaussian on a logit scale with 
#'       SD given by sigmaEpsilonDraws
#' list(pixelPop=outPixelLevel, subareaPop=outSubareaLevel, areaPop=outAreaLevel, uDraws=uDraws)
#' @return The simulated population aggregated to the enumeration area, 
#' pixel, subarea (generally Admin2), and area (generally Admin1) levels. Output includes:
#' \item{pixelPop}{A list of pixel level population aggregates}
#' \item{subareaPop}{A list of `subarea` level population aggregates}
#' \item{areaPop}{A list of `area` level population aggregates}
#' Each of these contains population numerator and denominator as well as prevalence and risk 
#' information aggregated to the appropriate level.
#' 
#' @details For population simulation and aggregation, we consider two models: the 
#' fine scale risk and the fine scale prevalence model. Both will be described in detail 
#' in a paper in preparation. In both models, enumeration areas (EAs) are simulated as 
#' point locations using an inhomogeneous Poisson process, with rate proportional to 
#' population density. EAs and populations are dispersed conditional on the (possibly 
#' approximately) known number of EAs, households, and target population at a particular 
#' areal level (these we call `areas`) using multilevel multinomial sampling, first 
#' sampling the EAs, then distributing households among the EAs, then the target population 
#' among the households. Any areal level below the `areas` we call `subareas`. For instance, 
#' the `areas` might be Admin1 if that is the smallest level at which the number of EAs, 
#' households, and people is known, and the `subareas` might be Admin2. The multilevel 
#' multinomial sampling may be stratified by urban/rural within the areas if the number of 
#' EAs, households, and people is also approximately known at that level.
#' 
#' Within each EA we assume a fixed probability of an event occurring, which is the `risk`. 
#' The `prevalence` is the empirical proportion of events within that EA. We assume EA 
#' level logit scale iid N(0, sigmaEpsilon^2) random effects in the risk model. When averaged 
#' with equal weights over all EAs in an areal unit, this forms the fine scale risk. When 
#' instead the population numerators and denominators are aggregated, and are used to 
#' calculate the empirical proportion of events occurring in an areal unit, the resulting 
#' quantity is the fine scale prevalence in that areal unit.
#' 
#' Note that these functions can be used for either simulating populations for simulation 
#' studies, or for generating predictions accounting for uncertainty in EA locations 
#' and fine scale variation occuring at the EA level due to EA level iid random effects. 
#' Required, however, is a seperatelt fit EA level spatial risk model 
#' and information on the spatial population density and the population frame.
#' 
#' @author John Paige
#' @references In preparation
#' @seealso \code{\link{simPopPixel}}, \code{\link{makeInterpPopMat}}, \code{\link{adjustPopMat}}, \code{\link{simSPDE}}.
#' @examples 
#' \dontrun{
#' ##### Now we make a model for the risk. We will use an SPDE model with these 
#' ##### parameters for the linear predictor on the logist scale, which are chosen 
#' ##### to be of practical interest:
#' beta0=-2.9 # intercept
#' gamma=-1 # urban effect
#' rho=(1/3)^2 # spatial variance
#' effRange = 400 # effective spatial range in km
#' sigmaEpsilon=sqrt(1/2.5) # cluster (nugget) effect standard deviation
#' 
#' # simulate the population! Note that this produces multiple dense nEA x nsim and nPixel x nsim 
#' # matrices. In the future sparse matrices will and chunk by chunk computations may be incorporated.
#' data(kenyaMaps)
#' data(kenyaPopulationData)
#' simPop = simPopSPDE(nsim=1, easpa=easpaKenyaNeonatal, 
#'                     popMat=popMatKenya, targetPopMat=popMatKenyaNeonatal, 
#'                     poppsub=poppsubKenya, spdeMesh=kenyaMesh, 
#'                     margVar=rho, sigmaEpsilonSq=sigmaEpsilon^2, 
#'                     gamma=gamma, effRange=effRange, beta0=beta0, 
#'                     seed=123, inla.seed=12, nHHSampled=25, 
#'                     stratifyByUrban=TRUE, subareaLevel=TRUE, 
#'                     doFineScaleRisk=TRUE, 
#'                     min1PerSubarea=TRUE)
#'                     
#' # plot EA level (aggregated to pixels by quilt.plot) simulated risk and prevalence
#' eaDat = simPop$pixelPop$eaDat[[1]]
#' require(fields)
#' quilt.plot(eaDat$lon, eaDat$lat, eaDat$Z) # population numerator
#' quilt.plot(eaDat$lon, eaDat$lat, eaDat$N) # population denominator
#' quilt.plot(eaDat$lon, eaDat$lat, eaDat$pFineScalePrevalence) # numerator/denominator
#' quilt.plot(eaDat$lon, eaDat$lat, eaDat$pFineScaleRisk)
#' quilt.plot(eaDat$lon, eaDat$lat, eaDat$pFineScaleRisk, zlim=c(0, .09))
#' 
#' # We can simulate the above population in a more general way 
#' # by first simulating from the SPDE model for risk, and then 
#' # inputting that risk to a more general population simulator:
#' set.seed(123)
#' pixelCoords = cbind(popMatKenya$east, popMatKenya$north)
#' 
#' SPDEArgs = list(coords=pixelCoords, nsim=1, margVar=rho, effRange=effRange, 
#'                 mesh=kenyaMesh, inla.seed=12L)
#' simVals = do.call("simSPDE", SPDEArgs)
#' 
#' # add in intercept
#' simVals = simVals + beta0
#' 
#' # add in urban effect
#' simVals = sweep(simVals, 1, gamma*popMatKenya$urban, "+")
#' 
#' # simulate nugget/cluster effect
#' totalEAs = sum(easpaKenyaNeonatal$EATotal)
#' epsc = matrix(stats::rnorm(totalEAs, sd=sigmaEpsilon), ncol=1)
#' 
#' # transform back to original scale for the pixel level probabilities
#' probsNoNug = expit(simVals)
#' 
#' # simulate the enumeration areas
#' uDraws = simVals
#' sigmaEpsilonDraws = sigmaEpsilon
#' 
#' # simulate EA and pixel level populations given the risk 
#' # simulated by the SPDE model
#' pixelPop = simPopPixel(uDraws=uDraws, sigmaEpsilonDraws=sigmaEpsilonDraws, 
#'                        easpa=easpaKenyaNeonatal, 
#'                        popMat=popMatKenya, targetPopMat=popMatKenyaNeonatal, 
#'                        stratifyByUrban=TRUE, 
#'                        doFineScaleRisk=TRUE, poppsub=poppsubKenya, 
#'                        subareaLevel=TRUE, min1PerSubarea=TRUE, 
#'                        returnEAinfo=TRUE, epsc=epsc)
#' }
#' @name simPop
NULL

#' @describeIn simPop
#' Simulate populations and population prevalences given census frame and population density 
#' information. Uses SPDE model for generating spatial risk and can include iid cluster 
#' level effect.
#' 
#' @export
simPopSPDE = function(nsim=1, easpa, popMat, targetPopMat, poppsub, spdeMesh, 
                      margVar=0.243, sigmaEpsilon=sqrt(0.463), 
                      gamma=0.009, effRange=406.51, beta0=-3.922, 
                      seed=NULL, inla.seed=-1L, nHHSampled=25, 
                      stratifyByUrban=TRUE, subareaLevel=TRUE, 
                      doFineScaleRisk=FALSE, 
                      min1PerSubarea=TRUE) {
  if(!is.null(seed))  {
    set.seed(seed)
    
    if(inla.seed < 0) {
      stop("seed specified, but not inla.seed. Set inla.seed to a positive integer to ensure reproducibility")
    }
  }
  
  if(nsim > 1) {
    warning("nsim > 1. eaDat will only be generated for first simulation")
  }
  
  totalEAs = sum(easpa$EATotal)
  totalHouseholds = sum(easpa$HHTotal)
  
  ### generate Binomial probabilities from transformed logit scale GP
  # generate SPDE simulations
  pixelCoords = cbind(popMat$east, popMat$north)
  
  if(margVar != 0) {
    SPDEArgs = list(coords=pixelCoords, nsim=nsim, margVar=margVar, effRange=effRange, 
                    mesh=spdeMesh, inla.seed=inla.seed)
    simVals = do.call("simSPDE", SPDEArgs)
  } else {
    simVals = matrix(rep(0, nrow(pixelCoords)), ncol=1)
  }
  
  # add in intercept
  simVals = simVals + beta0
  
  # add in urban effect
  simVals = sweep(simVals, 1, gamma*popMat$urban, "+")
  
  # simulate nugget/cluster effect
  epsc = matrix(stats::rnorm(totalEAs*nsim, sd=sigmaEpsilon), ncol=nsim)
  
  # transform back to original scale for the pixel level probabilities
  probsNoNug = expit(simVals)
  
  # simulate the enumeration areas
  uDraws = simVals
  sigmaEpsilonDraws = rep(sigmaEpsilon, nsim)
  
  print("Using SPDE model to simulate EA and pixel level populations")
  outPixelLevel = simPopPixel(uDraws=uDraws, sigmaEpsilonDraws=sigmaEpsilonDraws, easpa=easpa, 
                              popMat=popMat, targetPopMat=targetPopMat, 
                              stratifyByUrban=stratifyByUrban, 
                              doFineScaleRisk=doFineScaleRisk, poppsub=poppsub, 
                              subareaLevel=subareaLevel, min1PerSubarea=min1PerSubarea, 
                              returnEAinfo=TRUE, epsc=epsc)
  
  
  if(subareaLevel) {
    print("aggregating from pixel level to subarea level")
    outSubareaLevel = pixelPopToArea(pixelLevelPop=outPixelLevel, eaSamples=outPixelLevel$eaSamples, 
                                     areas=popMat$subarea, stratifyByUrban=stratifyByUrban, 
                                     targetPopMat=targetPopMat, doFineScaleRisk=doFineScaleRisk)
    
    print("aggregating from subarea level to area level")
    
    # get areas associated with each subarea for aggregation
    tempAreasFrom = popMat$subarea
    tempAreasTo = popMat$area
    areasFrom = sort(unique(tempAreasFrom))
    areasToI = match(areasFrom, tempAreasFrom)
    areasTo = tempAreasTo[areasToI]
    
    # do the aggregation from subareas to areas
    outAreaLevel = areaPopToArea(areaLevelPop=outSubareaLevel, areasFrom=areasFrom, areasTo=areasTo, 
                                 stratifyByUrban=stratifyByUrban, doFineScaleRisk=doFineScaleRisk)
  } else {
    outSubareaLevel = NULL
    
    print("aggregating from pixel level to area level")
    outAreaLevel = pixelPopToArea(pixelLevelPop=outPixelLevel, eaSamples=outPixelLevel$eaSamples, 
                                  areas=popMat$area, stratifyByUrban=stratifyByUrban, 
                                  doFineScaleRisk=doFineScaleRisk)
  }
  
  list(pixelPop=outPixelLevel, subareaPop=outSubareaLevel, areaPop=outAreaLevel, uDraws=uDraws)
}

#' Aggregate populations to the specified areal level
#' 
#' Takes simulated populations and aggregates 
#' them to the specified areal level. Also calculates the aggregated risk and prevalence.
#' 
#' @param pixelLevelPop pixel level population information that we want aggregate. In the same format as output from \code{\link{simPopPixel}}
#' @param eaSamples nPixel x nsim matrix of the number of enumeration areas per pixel sampled in the input pixel level population
#' @param areas character vector of length nPixels of area names over which we 
#'        want to aggregate. Can also be subareas
#' @param stratifyByUrban whether or not to stratify simulations by urban/rural classification
#' @param targetPopMat pixellated grid data frame with variables `lon`, `lat`, `pop` (target population), `area`, `subareas` (if subareaLevel is TRUE), `urban` (if stratifyByUrban is TRUE), `east`, and `north`
#' @param doFineScaleRisk whether or not to calculate the fine scale risk in addition to the prevalence
#' @param easpa see \code{\link{simPopSPDE}}
#' @param areaLevelPop output of \code{\link{simPopPixel}} containing pixel level information 
#'               about the population of interest
#' @param areasFrom character vector of length equal to the number of areas from which 
#'            we would like to aggregate containing the unique names of the areas. 
#'            Can also be subareas, but these are smaller than the "to areas", and 
#'            each "from area" but be entirely contained in a single "to area"
#' @param areasTo character vector of length equal to the number of subareas to which 
#'          we would like to aggregate containing the names of the areas associated 
#'          with each respective subarea. Can also be a different set of subareas, 
#'          but these are larger than the "from areas".
#' 
#' @author John Paige
#' @references In Preparation
#' @return A list containing elements `fineScalePrevalence` and `fineScaleRisk`. Each 
#' of these are in turn lists with aggregated prevalence and risk for the area of 
#' interest, containg the following elements, were paranethesis indicate the elements 
#' for the fineScaleRisk model rather than fineScalePrevalence:
#' \item{p}{Aggregated prevalence (risk), calculated as aggregate of Z divided by 
#' aggregate of N}
#' \item{Z}{Aggregated (expected) population numerator}
#' \item{N}{Aggregated (expected) population denominator}
#' \item{pUrban}{Aggregated prevalence (risk) in urban part of the area, calculated 
#' as aggregate of Z divided by aggregate of N}
#' \item{ZUrban}{Aggregated (expected) population numerator in urban part of the area}
#' \item{NUrban}{Aggregated (expected) population denominator in urban part of the area}
#' \item{pRural}{Aggregated prevalence (risk) in rural part of the area, calculated 
#' as aggregate of Z divided by aggregate of N}
#' \item{ZRural}{Aggregated (expected) population numerator in rural part of the area}
#' \item{NRural}{Aggregated (expected) population denominator in rural part of the area}
#' \item{A}{Aggregation matrix used to aggregate from pixel level to areal level}
#' \item{AUrban}{Aggregation matrix used to aggregate from pixel level to urban part of the areal level}
#' \item{ARural}{Aggregation matrix used to aggregate from pixel level to rural part of the areal level}
#' @seealso \code{\link{areaPopToArea}}
#' @name aggPop
#' @examples
#' \dontrun{
#' ##### Now we make a model for the risk. We will use an SPDE model with these 
#' ##### parameters for the linear predictor on the logist scale, which are chosen 
#' ##### to be of practical interest:
#' beta0=-2.9 # intercept
#' gamma=-1 # urban effect
#' rho=(1/3)^2 # spatial variance
#' effRange = 400 # effective spatial range in km
#' sigmaEpsilon=sqrt(1/2.5) # cluster (nugget) effect standard deviation
#' 
#' # simulate the population! Note that this produces multiple dense nEA x nsim and nPixel x nsim 
#' # matrices. In the future sparse matrices will and chunk by chunk computations may be incorporated.
#' simPop = simPopSPDE(nsim=1, easpa=easpaKenyaNeonatal, 
#'                     popMat=popMatKenya, targetPopMat=popMatKenyaNeonatal, 
#'                     poppsub=poppsubKenya, spdeMesh=kenyaMesh, 
#'                     margVar=rho, sigmaEpsilonSq=sigmaEpsilon^2, 
#'                     gamma=gamma, effRange=effRange, beta0=beta0, 
#'                     seed=123, inla.seed=12, nHHSampled=25, 
#'                     stratifyByUrban=TRUE, subareaLevel=TRUE, 
#'                     doFineScaleRisk=TRUE, 
#'                     min1PerSubarea=TRUE)
#' 
#' pixelPop = simPop$pixelPop
#' subareaPop = pixelPopToArea(pixelLevelPop=pixelPop, eaSamples=pixelPop$eaSamples, 
#'   areas=popMatKenya$subarea, stratifyByUrban=TRUE, 
#'   targetPopMat=popMatKenyaNeonatal, doFineScaleRisk=TRUE)
#' 
#' # get areas associated with each subarea for aggregation
#' tempAreasFrom = popMatKenya$subarea
#' tempAreasTo = popMatKenya$area
#' areasFrom = sort(unique(tempAreasFrom))
#' areasToI = match(areasFrom, tempAreasFrom)
#' areasTo = tempAreasTo[areasToI]
#' 
#' # do the aggregation from subareas to areas
#' outAreaLevel = areaPopToArea(areaLevelPop=subareaPop, 
#'   areasFrom=areasFrom, areasTo=areasTo, 
#'   stratifyByUrban=TRUE, doFineScaleRisk=TRUE)
#' }
NULL

#' @describeIn aggPop Aggregate from pixel to areal level
#' @export
pixelPopToArea = function(pixelLevelPop, eaSamples, areas, stratifyByUrban=TRUE, targetPopMat=NULL, 
                          doFineScaleRisk=!is.null(pixelLevelPop$fineScaleRisk$p), easpa=NULL) {
  
  # fine scale prevalence aggregation model
  nSamples = pixelLevelPop$fineScalePrevalence$N
  zSamples = pixelLevelPop$fineScalePrevalence$Z
  zSamples[is.na(zSamples)] = 0 # must set to zero temporarily so matrix multiplication works
  aggregatedResults = aggPixelPreds(Zg=zSamples, Ng=nSamples, areas=areas, targetPopMat=targetPopMat, 
                                    useDensity=FALSE, stratifyByUrban=stratifyByUrban, normalize=FALSE)
  
  # fine scale risk aggregation model
  if(doFineScaleRisk) {
    
    # in order to get valid count estimates, we also need the expected denominator per EA in each stratum:
    nPerEA = getExpectedNperEA(easpa, targetPopMat)
    nSamplesFineScaleRisk = sweep(eaSamples, 1, nPerEA, "*")
    zSamplesFineScaleRisk = pixelLevelPop$fineScaleRisk$p * nSamplesFineScaleRisk
    zSamplesFineScaleRisk[is.na(zSamplesFineScaleRisk)] = 0 # must set to zero temporarily so matrix multiplication works out
    aggregatedResultsFineScaleRisk = aggPixelPreds(Zg=zSamplesFineScaleRisk, Ng=nSamplesFineScaleRisk, areas=areas, targetPopMat=targetPopMat, 
                                                   useDensity=FALSE, stratifyByUrban=stratifyByUrban, normalize=FALSE)
  } else {
    nSamplesFineScaleRisk = NULL
    zSamplesFineScaleRisk = NULL
    aggregatedResultsFineScaleRisk = NULL
  }
  
  areaLevelPop=list(fineScalePrevalence=aggregatedResults, 
                    fineScaleRisk=aggregatedResultsFineScaleRisk)
  areaLevelPop
}

#' Helper function of \code{\link{pixelPopToArea}}
#' 
#' Aggregates population from the 
#' pixel level to the level of the area of interest.
#' 
#' @param Zg nPixel x nsim matrix of simulated response (population numerators) for each pixel and sample
#' @param Ng nPixel x nsim matrix of simulated counts (population denominators) for each pixel and sample
#' @param areas nPixel length character vector of areas (or subareas) 
#' @param urban nPixel length vector of indicators specifying whether or not pixels are urban or rural
#' @param targetPopMat same as in \code{\link{simPopPixel}}
#' @param useDensity whether to use population density as aggregation weights. 
#' @param stratifyByUrban whether or not to stratify simulations by urban/rural classification
#' @param normalize if TRUE, pixel level aggregation weights within specified area are normalized to sum to 1. This produces an 
#' average of the values in Zg rather than a sum. In general, should only be set to TRUE for smooth integrals of risk.
aggPixelPreds = function(Zg, Ng, areas, urban=targetPopMat$urban, targetPopMat=NULL, useDensity=FALSE, 
                         stratifyByUrban=TRUE, normalize=useDensity) {
  
  if(useDensity && !normalize) {
    stop("if useDensity is set to TRUE, normalize must be set to TRUE as well")
  }
  predsUrban = urban
  predsArea = areas
  
  # set NAs and pixels without any sample size to 0
  Ng[is.na(Ng)] = 0
  if(!useDensity) {
    # is useDensity is true, then Zg is really a set of probabilities, so no need to set to 0
    Zg[Ng == 0] = 0
  }
  
  # function to aggregate predictions over the 
  # population density grid. Use the following function to get numerical 
  # integration matrix for a given level of areal aggregation. returned 
  # matrices have dimension length(unique(areaNames)) x length(areaNames)
  # areaNames: 
  # urbanProportions: DEPRACATED vector giving proportion of population urban for each unique area in areaNames. 
  #                   If specified, ensure that urban and rural parts of the full integration 
  #                   matrix have the appropriate relative weights for each area. Used for population 
  #                   density based integration
  # normalize: whether or not to normalize the rows of the matrices to sum to 1 or to instead 
  #            contain only binary values (or non-binary values based on the binary values if 
  #            urbanProportions is not NULL)
  getIntegrationMatrix = function(areaNames, urbanProportions=NULL, normalize=FALSE) {
    
    if(useDensity) {
      popDensities = targetPopMat$pop
      densities = popDensities
    } else {
      equalDensities = rep(1, nrow(Zg))
      densities = equalDensities
    }
    
    uniqueNames = sort(unique(areaNames))
    getMatrixHelper = function(i, thisUrban=NULL, thisUseDensity=useDensity, thisNormalize=normalize) {
      areaI = areaNames == uniqueNames[i]
      
      if(thisUseDensity) {
        theseDensities = popDensities
      } else {
        theseDensities = equalDensities
      }
      
      # make sure we only include pixels in the given area and, if necessary, with the given urbanicity
      theseDensities[!areaI] = 0
      if(!is.null(thisUrban))
        theseDensities[predsUrban != thisUrban] = 0
      thisSum = sum(theseDensities)
      if(thisSum != 0 && thisNormalize)
        theseDensities * (1/thisSum)
      else if(thisSum == 0)
        rep(0, length(theseDensities))
      else
        theseDensities
    }
    
    if(!stratifyByUrban) {
      integrationMatrix = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper), ncol=length(uniqueNames)))
      
      integrationMatrix
    } else {
      integrationMatrixUrban = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper, thisUrban=TRUE), ncol=length(uniqueNames)))
      integrationMatrixRural = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper, thisUrban=FALSE), ncol=length(uniqueNames)))
      if(!is.null(urbanProportions)) {
        integrationMatrix = sweep(integrationMatrixUrban, 1, urbanProportions, "*") + sweep(integrationMatrixRural, 1, 1-urbanProportions, "*")
      } else {
        integrationMatrix = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper), ncol=length(uniqueNames)))
      }
      
      rownames(integrationMatrix) = uniqueNames
      rownames(integrationMatrixUrban) = uniqueNames
      rownames(integrationMatrixRural) = uniqueNames
      
      list(integrationMatrix=integrationMatrix, 
           integrationMatrixUrban=integrationMatrixUrban, 
           integrationMatrixRural=integrationMatrixRural)
    }
  }
  
  # Use the following function to perform the
  # aggregations
  getIntegratedPredictions = function(areaNames) {
    # get numerical integration matrix
    A = getIntegrationMatrix(areaNames, normalize=normalize)
    
    # aggregate the prediction and denominator matrices (for whole areas and also urban/rural strata if necessary)
    if(!stratifyByUrban) {
      ZAggregated = A %*% Zg
      NAggregated = A %*% Ng
      pAggregated = ZAggregated / NAggregated
      pAggregated[NAggregated == 0] = NA
      
      list(p=pAggregated, Z=ZAggregated, N=NAggregated, A=A)
    } else {
      AUrban = A$integrationMatrixUrban
      ARural = A$integrationMatrixRural
      A = A$integrationMatrix
      
      # first aggregate the numerator. The denominator will depend on the aggregation method
      ZAggregated = A %*% Zg
      ZAggregatedUrban = AUrban %*% Zg
      ZAggregatedRural = ARural %*% Zg
      
      if(useDensity) {
        # for population density aggregation, we integrate probabilities rather than aggregate 
        # empirical proportions
        NAggregated = NULL
        pAggregated = ZAggregated
        ZAggregated = NULL
        
        NAggregatedUrban = NULL
        pAggregatedUrban = ZAggregatedUrban
        ZAggregatedUrban = NULL
        
        NAggregatedRural = NULL
        pAggregatedRural = ZAggregatedRural
        ZAggregatedRural = NULL
      } else {
        # if we do not use density, we must also aggregate the denominator to calculate 
        # the aggregated empirical proportions
        NAggregated = A %*% Ng
        pAggregated = ZAggregated / NAggregated
        pAggregated[NAggregated == 0] = NA
        
        NAggregatedUrban = AUrban %*% Ng
        pAggregatedUrban = ZAggregatedUrban / NAggregatedUrban
        pAggregatedUrban[NAggregatedUrban == 0] = NA
        
        NAggregatedRural = ARural %*% Ng
        pAggregatedRural = ZAggregatedRural / NAggregatedRural
        pAggregatedRural[NAggregatedRural == 0] = NA
      }
      
      list(p=pAggregated, Z=ZAggregated, N=NAggregated, 
           pUrban=pAggregatedUrban, ZUrban=ZAggregatedUrban, NUrban=NAggregatedUrban, 
           pRural=pAggregatedRural, ZRural=ZAggregatedRural, NRural=NAggregatedRural, 
           A=A, AUrban=AUrban, ARural=ARural)
    }
  }
  
  areaMatrices = getIntegratedPredictions(predsArea)
  
  # return results
  areaMatrices
}

#' @describeIn aggPop Aggregate areal populations to another areal level
#' @export
areaPopToArea = function(areaLevelPop, areasFrom, areasTo, stratifyByUrban=TRUE, 
                         doFineScaleRisk=!is.null(areaLevelPop$fineScaleRisk$p)) {
  
  if(length(areasFrom) != length(unique(areasFrom))) {
    stop("areasFrom must contain only unique names of areas to which we want to aggregate")
  }
  
  uniqueNames = sort(unique(areasTo))
  
  # construct row of the aggregation matrix given toArea index
  getMatrixHelper = function(i) {
    # get areasTo associated with this fromArea
    thisToArea = uniqueNames[i]
    thisFromAreas = unique(areasFrom[areasTo == thisToArea])
    areaI = areasFrom %in% thisFromAreas
    
    areaI
  }
  
  # construct the aggregation matrix from areasFrom to areasTo
  A = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper), ncol=length(uniqueNames)))
  rownames(A) = uniqueNames
  colnames(A) = areasFrom
  
  ##### aggregate populations
  # fine scale prevalence aggregation model
  nSamples = areaLevelPop$fineScalePrevalence$N
  zSamples = areaLevelPop$fineScalePrevalence$Z
  zSamples[is.na(zSamples)] = 0 # must set to zero temporarily so matrix multiplication works out
  
  ZAggregated =  A %*% zSamples
  NAggregated =  A %*% nSamples
  pAggregated = ZAggregated / NAggregated
  pAggregated[NAggregated == 0] = NA
  thisA=A %*% areaLevelPop$fineScalePrevalence$A
  rownames(thisA) = uniqueNames
  
  if(stratifyByUrban) {
    nSamplesUrban = areaLevelPop$fineScalePrevalence$NUrban
    zSamplesUrban = areaLevelPop$fineScalePrevalence$ZUrban
    zSamplesUrban[is.na(zSamplesUrban)] = 0 # must set to zero temporarily so matrix multiplication works out
    
    nSamplesRural = areaLevelPop$fineScalePrevalence$NRural
    zSamplesRural = areaLevelPop$fineScalePrevalence$ZRural
    zSamplesRural[is.na(zSamplesRural)] = 0 # must set to zero temporarily so matrix multiplication works out
    
    ZAggregatedUrban =  A %*% zSamplesUrban
    NAggregatedUrban =  A %*% nSamplesUrban
    pAggregatedUrban = ZAggregatedUrban / NAggregatedUrban
    pAggregatedUrban[NAggregatedUrban == 0] = NA
    thisAUrban=A %*% areaLevelPop$fineScalePrevalence$AUrban
    rownames(thisAUrban) = uniqueNames
    
    ZAggregatedRural =  A %*% zSamplesRural
    NAggregatedRural =  A %*% nSamplesRural
    pAggregatedRural = ZAggregatedRural / NAggregatedRural
    pAggregatedRural[NAggregatedRural == 0] = NA
    thisARural=A %*% areaLevelPop$fineScalePrevalence$ARural
    rownames(thisARural) = uniqueNames
  } else {
    ZAggregatedUrban = NULL
    NAggregatedUrban = NULL
    pAggregatedUrban = NULL
    thisAUrban=NULL
    
    ZAggregatedRural = NULL
    NAggregatedRural = NULL
    pAggregatedRural = NULL
    thisARural=NULL
  }
  
  aggregatedResults = list(p=pAggregated, Z=ZAggregated, N=NAggregated, 
                           pUrban=pAggregatedUrban, ZUrban=ZAggregatedUrban, NUrban=NAggregatedUrban, 
                           pRural=pAggregatedRural, ZRural=ZAggregatedRural, NRural=NAggregatedRural, 
                           A=thisA, 
                           AUrban=thisAUrban, 
                           ARural=thisARural)
  
  # fine scale risk aggregation model
  nSamples = areaLevelPop$fineScaleRisk$N
  zSamples = areaLevelPop$fineScaleRisk$Z
  zSamples[is.na(zSamples)] = 0 # must set to zero temporarily so matrix multiplication works out
  
  ZAggregated =  A %*% zSamples
  NAggregated =  A %*% nSamples
  pAggregated = ZAggregated / NAggregated
  pAggregated[NAggregated == 0] = NA
  thisA=A %*% areaLevelPop$fineScalePrevalence$A
  rownames(thisA) = uniqueNames
  
  if(stratifyByUrban) {
    nSamplesUrban = areaLevelPop$fineScalePrevalence$NUrban
    zSamplesUrban = areaLevelPop$fineScalePrevalence$ZUrban
    zSamplesUrban[is.na(zSamplesUrban)] = 0 # must set to zero temporarily so matrix multiplication works out
    
    nSamplesRural = areaLevelPop$fineScalePrevalence$NRural
    zSamplesRural = areaLevelPop$fineScalePrevalence$ZRural
    zSamplesRural[is.na(zSamplesRural)] = 0 # must set to zero temporarily so matrix multiplication works out
    
    ZAggregatedUrban =  A %*% zSamplesUrban
    NAggregatedUrban =  A %*% nSamplesUrban
    pAggregatedUrban = ZAggregatedUrban / NAggregatedUrban
    pAggregatedUrban[NAggregatedUrban == 0] = NA
    thisAUrban=A %*% areaLevelPop$fineScalePrevalence$AUrban
    rownames(thisAUrban) = uniqueNames
    
    ZAggregatedRural =  A %*% zSamplesRural
    NAggregatedRural =  A %*% nSamplesRural
    pAggregatedRural = ZAggregatedRural / NAggregatedRural
    pAggregatedRural[NAggregatedRural == 0] = NA
    thisARural=A %*% areaLevelPop$fineScalePrevalence$ARural
    rownames(thisARural) = uniqueNames
  } else {
    ZAggregatedUrban = NULL
    NAggregatedUrban = NULL
    pAggregatedUrban = NULL
    thisAUrban=NULL
    
    ZAggregatedRural = NULL
    NAggregatedRural = NULL
    pAggregatedRural = NULL
    thisARural=NULL
  }
  
  aggregatedResultsRisk = list(p=pAggregated, Z=ZAggregated, N=NAggregated, 
                               pUrban=pAggregatedUrban, ZUrban=ZAggregatedUrban, NUrban=NAggregatedUrban, 
                               pRural=pAggregatedRural, ZRural=ZAggregatedRural, NRural=NAggregatedRural, 
                               A=thisA, 
                               AUrban=thisAUrban, 
                               ARural=thisARural)
  
  list(aggregatedResultsPrevalence=aggregatedResults, 
       aggregatedResultsRisk=aggregatedResultsRisk)
}

#' @describeIn simPop
#' Simulate populations and population prevalences given census frame and population density 
#' information. Uses custom spatial logit risk function and can include iid cluster 
#' level effect.
#' 
#' @export
simPopPixel = function(uDraws, sigmaEpsilonDraws, easpa, popMat, targetPopMat, 
                       stratifyByUrban=TRUE, validationPixelI=NULL, validationClusterI=NULL, 
                       clustersPerPixel=NULL, 
                       doFineScaleRisk=FALSE, poppsub=NULL, subareaLevel=FALSE, 
                       min1PerSubarea=TRUE, 
                       returnEAinfo=FALSE, epsc=NULL) {
  
  if(is.null(poppsub) && subareaLevel) {
    stop("if subareaLevel is TRUE, user must specify poppsub")
  }
  
  if(!is.null(validationPixelI) || !is.null(validationClusterI) || !is.null(clustersPerPixel)) {
    stop("validationPixelI, validationClusterI, and clustersPerPixel not yet fully implemented")
  }
  
  nDraws = ncol(uDraws)
  
  # set default inputs
  totalEAs = sum(easpa$EATotal)
  if(!is.null(clustersPerPixel)) {
    emptyPixels = clustersPerPixel == 0
    if(totalEAs != sum(clustersPerPixel))
      stop("sum(easpa$EATotal) != sum(clustersPerPixel)")
  }
  
  # get area names
  areas = sort(unique(popMat$area))
  if(any(areas != easpa$area))
    stop("area names and easpa do not match popMat or are not in the correct order")
  
  # determine if we care about subareas (smallest areas we care about. No info of EAs per subarea)
  subareaLevel = !is.null(popMat$subarea)
  
  ##### Line 1 (of the algorithm): take draws from the binomial process for each stratum (each row of easpa)
  # get probabilities for each pixel (or at least something proportional within each stratum)
  print("drawing EAs")
  pixelProbs = popMat$pop
  
  # take draws from the stratified binomial process for each posterior sample
  if(is.null(clustersPerPixel)) {
    if(subareaLevel) {
      eaSamples = rStratifiedMultnomialBySubarea(nDraws, popMat, easpa, stratifyByUrban, poppsub=poppsub, 
                                                 min1PerSubarea=min1PerSubarea)
    } else {
      eaSamples = rStratifiedMultnomial(nDraws, popMat, easpa, stratifyByUrban)
    }
  }
  
  if(!is.null(clustersPerPixel) && !exists("eaSamples")) {
    eaSamples = matrix(rep(clustersPerPixel, nDraws), ncol=nDraws)
  }
  
  # make matrix (or list) of pixel indices mapping matrices of EA values to matrices of pixel values
  if(!is.null(clustersPerPixel)) {
    pixelIndices = rep(1:nrow(popMat), times=clustersPerPixel) # this contains repetitions and has length == nEAs
    uniquePixelIndices = sort(unique(pixelIndices))
  } else {
    pixelIndexMat = matrix(rep(rep(1:nrow(popMat), nDraws), times=eaSamples), ncol=nDraws)
  }
  
  # determine which EAs are urban if necessary
  if(stratifyByUrban) {
    # urbanMat = matrix(rep(rep(popMat$urban, nDraws), times=c(eaSamples)), ncol=nDraws)
    if(!is.null(clustersPerPixel)) {
      urbanVals = popMat$urban[pixelIndices]
      uniqueUrbanVals = popMat$urban[uniquePixelIndices]
    }else {
      urbanMat = matrix(popMat$urban[pixelIndexMat], ncol=nDraws)
    }
  } else {
    urbanMat = NULL
  }
  
  # determine which EAs are from which area
  if(!is.null(clustersPerPixel)) {
    areaVals = popMat$area[pixelIndices]
    uniqueAreaVals = popMat$area[uniquePixelIndices]
  } else {
    areaMat = matrix(popMat$area[pixelIndexMat], ncol=nDraws)
  }
  
  ##### Line 2: draw cluster effects, epsilon
  # NOTE1: we assume there are many more EAs then sampled clusters, so that 
  #       the cluster effects for each EA, including those sampled, are iid
  print("simulating EA level risks, numerators, and denominators")
  if(is.null(epsc)) {
    epsc = matrix(stats::rnorm(totalEAs*nDraws, sd=rep(sigmaEpsilonDraws, each=totalEAs)), ncol=nDraws)
  }
  
  ##### Line 3: draw EA population denominators, N
  
  if(!is.null(clustersPerPixel)) {
    if(is.null(validationPixelI))
      stop("clustersPerPixel must only be set for validation, but validationPixelI is NULL")
    
    # in this case, every left out cluster has exactly 25 households. Simply sample target population 
    # with equal probability from each cluster/faux EA
    Ncs = sampleNMultilevelMultinomialFixed(clustersPerPixel, nDraws=nDraws, pixelIndices=pixelIndices, 
                                            urbanVals=urbanVals, areaVals=areaVals, easpa=easpa, popMat=popMat, stratifyByUrban=stratifyByUrban, 
                                            verbose=TRUE)
  } else {
    if(returnEAinfo) {
      out = sampleNMultilevelMultinomial(pixelIndexMat=pixelIndexMat, urbanMat=urbanMat, areaMat=areaMat, easpaList=list(easpa), 
                                         popMat=popMat, stratifyByUrban=stratifyByUrban, verbose=TRUE, returnEAinfo=returnEAinfo)
      householdDraws = out$householdDraws
      Ncs = out$targetPopDraws
    } else {
      Ncs <- sampleNMultilevelMultinomial(pixelIndexMat=pixelIndexMat, urbanMat=urbanMat, areaMat=areaMat, easpaList=list(easpa), 
                                       popMat=popMat, stratifyByUrban=stratifyByUrban, verbose=TRUE, returnEAinfo=returnEAinfo)
      householdDraws = NULL
    }
  }
  
  ##### do part of Line 7 in advance
  # calculate mu_{ic} for each EA in each pixel
  if(!is.null(clustersPerPixel)) {
    uc = uDraws[pixelIndices,]
    muc = expit(uc + epsc)
  } else {
    uc = matrix(uDraws[cbind(rep(rep(1:nrow(uDraws), nDraws), times=c(eaSamples)), rep(1:nDraws, each=totalEAs))], ncol=nDraws)
    muc = expit(uc + epsc)
  }
  
  # calculate Z_{ic} for each EA in each pixel
  Zc = matrix(stats::rbinom(n=totalEAs * nDraws, size=Ncs, prob=as.matrix(muc)), ncol=nDraws)
  
  ##### Line 4: Aggregate appropriate values from EAs to the grid cell level
  
  # function for aggregating values for each grid cell
  getPixelColumnFromEAs = function(i, vals, applyFun=sum, popWeightMatrix=NULL) {
    # calculate levels over which to aggregate
    if(!is.null(clustersPerPixel)) {
      indices = pixelIndices
    } else {
      indices = factor(as.character(pixelIndexMat[,i]))
    }
    
    # in this case (the LCPb model), we calculate weighted means within factor levels using popWeightMatrix
    if(!is.null(popWeightMatrix)) {
      stop("using popWeightMatrix is no longer support, since this is much slower than calculating normalized 
           weights separately and multiplying values by them outside this function")
      # applyFun = function(x) {stats::weighted.mean(x, popWeightMatrix[,i], na.rm=TRUE)}
      
      Data = data.frame(v=vals[,i], w=popWeightMatrix[,i])
      out = sapply(split(Data, indices), function(x) stats::weighted.mean(x$v,x$w))
    } else {
      if(!is.null(clustersPerPixel)) {
        # out = tapply(vals[,i], factor(as.character(pixelIndices)), FUN=applyFun)
        out = tapply(vals[,i], indices, FUN=applyFun)
      } else {
        out = tapply(vals[,i], indices, FUN=applyFun)
      }
    }
    
    if(!is.null(clustersPerPixel)) {
      returnValues = out
    } else {
      indices = as.numeric(names(out))
      
      returnValues = rep(NA, nrow(uDraws))
      returnValues[indices] = out
    }
    
    returnValues
  }
  
  ##### Line 5: We already did this, resulting in uDraws input
  
  ##### Line 6: aggregate population denominators for each grid cell to get N_{ig}
  print("Aggregating from EA level to the pixel level")
  Ng <- sapply(1:ncol(Ncs), getPixelColumnFromEAs, vals=Ncs)
  Ng[is.na(Ng)] = 0
  
  ##### Line 7: aggregate response for each grid cell to get Z_{ig}
  Zg <- sapply(1:ncol(Zc), getPixelColumnFromEAs, vals=Zc)
  
  ##### Line 8: Calculate empirical mortality proportions for each grid cell, p_{ig}. 
  #####         Whenever N_{ig} is 0, set p_{ig} to NA as well
  pg = Zg / Ng
  pg[Ng == 0] = NA
  
  ##### calculate results for the other models if necessary
  if(doFineScaleRisk) {
    fineScaleRisk = sapply(1:ncol(muc), getPixelColumnFromEAs, vals=muc, applyFun=function(x) {mean(x, na.rm=TRUE)})
    fineScaleRisk[!is.finite(fineScaleRisk)] = NA
    
    # in order to get valid count estimates, we also need the expected denominator per EA in each stratum:
    nPerEA = getExpectedNperEA(easpa, targetPopMat)
    nSamplesFineScaleRisk = sweep(eaSamples, 1, nPerEA, "*")
    zSamplesFineScaleRisk = fineScaleRisk * nSamplesFineScaleRisk
    zSamplesFineScaleRisk[is.na(zSamplesFineScaleRisk)] = 0
  } else {
    fineScaleRisk = NULL
    zSamplesFineScaleRisk = NULL
    nSamplesFineScaleRisk = NULL
  }
  
  ##### Extra steps: collect draws at each level and generate:
  ##### areas, preds, 
  pixelLevelPop = list(fineScalePrevalence=list(p=pg, Z=Zg, N=Ng), 
                       fineScaleRisk=list(p=fineScaleRisk, Z=zSamplesFineScaleRisk, N=nSamplesFineScaleRisk))
  
  if(!returnEAinfo) {
    pixelLevelPop
  } else {
    
    
    # return list of eaDat objects
    getEADat = function(i) {
      theseI = pixelIndexMat[,i]
      
      eaDat = data.frame(lon=popMat$lon[theseI], lat=popMat$lat[theseI], 
                         area=popMat$area[theseI], subarea=rep("temp", length(theseI)), 
                         urban=popMat$urban[theseI], east=popMat$east[theseI], north=popMat$north[theseI], 
                         popDensity=popMat$pop[theseI], popDensityTarget=targetPopMat$pop[theseI], pixelIs=theseI, 
                         nHH=householdDraws[,i], N=Ncs[,i], Z=Zc[,i], 
                         pFineScaleRisk=muc[,i], pFineScalePrevalence=Zc[,i]/Ncs[,i])
      if(subareaLevel) {
        eaDat$subarea = popMat$subarea[theseI]
      } else {
        eaDat$subarea = NULL
      }
      
      eaDat$pFineScalePrevalence[eaDat$n == 0] = NA
      
      eaDat
    }
    
    print("Constructing list of simulated EA data.frames")
    eaDatList = lapply(1:nDraws, getEADat)
    
    c(pixelLevelPop, list(eaDatList=eaDatList, eaSamples=eaSamples))
  }
}


#' Internal functions for population simulation
#' 
#' Functions for calculating valuable quantities and for drawing from important 
#' distributions for population simulation.
#' 
#' @param easpa census frame. See \code{\link{simPopPixel}} for details
#' @param popMat data.frame of pixellated grid of population densities. See \code{\link{simPopPixel}} for details
#' @param i index
#' @param urban if TRUE, calculate only for urban part of the area. If FALSE, for only rural part
#' @param stratifyByUrban whether or not to stratify calculations by urban/rural classification
#' @param validationPixelI CURRENTLY FOR TESTING PURPOSES ONLY a set of indices of pixels for which we want to simulate populations (used for pixel level validation)
#' @param n number of samples
#' @param poppsub population per subarea. See \code{\link{simPopPixel}} for details
#' @param min1PerSubarea whether or not to ensure there is at least 1 EA per subarea. See \code{\link{simPopPixel}} for details
#' @param method if min1PerSubarea is TRUE, the sampling method for the truncated multinomial to use with rmulitnom1. rmultinom1 automatically 
#'         switches between them depending on the number of expected samples. The methods are:
#' \describe{
#'   \item{mult1}{rejection sampling from multinomial plus 1 in each category}
#'   \item{mult}{rejection sampling from multinomial if any category has zero count}
#'   \item{indepMH}{independent Metropolis-Hastings using multinomial plus 1 distribution as proposal}
#' }
#' @param minSample the minimum number of samples per `chunk` of samples for truncated multinomial sampling. Defaults to 1
#' @param easpsub this could either be total EAs per subarea, or subarea crossed with urban or 
#'          rural if stratifyByUrban is TRUE
#' @param size multinomial size parameter. See \code{\link[stats]{rmultinom}}
#' @param prob multinomial probability vector parameter. See \code{\link[stats]{rmultinom}}
#' @param maxSize the maximum number of elements in a matrix drawn from the proposal distribution per sample chunk. 
#' @param maxExpectedSizeBeforeSwitch max expected number of samples / k, the number of categories, before switching method
#' @param init initial sample if method is `indepMH`
#' @param burnIn number of initial samples before samples are collected if method is `indepMH`
#' @param filterEvery store only every filterEvery samples if method is i`indepMH`
#' @param zeroProbZeroSamples if TRUE, set samples for parts of prob vector that are zero to zero. Otherwise they are set to one.
#' @param allowSizeLessThanK if TRUE, then if size < the number of categories (k), returns matrix where each 
#'                     column is vector of size ones and k - size zeros. If FALSE, throws an error if size < k
#' @param clustersPerPixel CURRENTLY FOR TESTING PURPOSES ONLY a vector of length nPixels specifying the number of clusters per pixel if they are fixed
#' @param pixelIndices a nEA x n matrix of pixel indices associated with each EA per simulation/draw
#' @param urbanVals a nEA x n matrix of urbanicities associated with each EA per simulation/draw
#' @param areaVals a nEA x n matrix of area names associated with each EA per simulation/draw
#' @param easpaList a list of length n with each element being of the format of easpa 
#'            giving the number of households and EAs 
#'            per stratum. It is assumed that the number of EAs per stratum is 
#'            the same in each list element. If easpaList is a data frame, 
#'            number of households per stratum is assumed constant
#' @param nDraws number of draws
#' @param pixelIndexMat matrix of pixel indices associated with each EA and draw
#' @param urbanMat matrix of urbanicities associated with each EA and draw
#' @param areaMat matrix of areas associated with each EA and draw
#' @param verbose whether to print progress as the function proceeds
#' @param returnEAinfo whether a data frame at the EA level is desired
#' 
#' @name simPopInternal
NULL

#' @describeIn simPopInternal Calculates expected denominator per enumeration area
getExpectedNperEA = function(easpa, popMat) {
  
  # calculate the expected denominator per enumeration area in each stratum. 
  nPerEAUrban = easpa$popUrb / easpa$EAUrb
  nPerEARural = easpa$popRur / easpa$EARur
  
  # expanded the expected denominator values victor to be of length equal 
  # to the number of grid cells
  uniqueAreas = sort(unique(popMat$area))
  outUrban = numeric(nrow(popMat))
  outRural = numeric(nrow(popMat))
  for(i in 1:length(uniqueAreas)) {
    urbanI = getSortIndices(i, urban=TRUE, popMat=popMat, stratifyByUrban=TRUE)
    ruralI = getSortIndices(i, urban=FALSE, popMat=popMat, stratifyByUrban=TRUE)
    outUrban[urbanI] = nPerEAUrban[i]
    outRural[ruralI] = nPerEARural[i]
  }
  
  outUrban + outRural
}

#' @describeIn simPopInternal For recombining separate multinomials into the draws over all pixels
getSortIndices = function(i, urban=TRUE, popMat, stratifyByUrban=TRUE, validationPixelI=NULL) {
  
  # get area names
  areas = sort(unique(popMat$area))
  
  # determine which pixels and how many EAs are in this stratum
  if(stratifyByUrban) {
    includeI = popMat$area == areas[i] & popMat$urban == urban
  }
  else {
    includeI = popMat$area == areas[i]
  }
  
  # include only indices included within validation if necessary
  if(!is.null(validationPixelI)) {
    
    # convert validationPixelI into a logical
    temp = rep(FALSE, length(includeI))
    temp[validationPixelI] = TRUE
    
    # include only indices we are interested in for the validation
    includeI = includeI & temp
  }
  
  which(includeI)
}

#' @describeIn simPopInternal Gives nPixels x n matrix of draws from the stratified multinomial with values 
#' corresponding to the value of |C^g| for each pixel, g (the number of EAs/pixel)
rStratifiedMultnomial = function(n, popMat, easpa, stratifyByUrban=TRUE) {
  
  # get area names
  areas = sort(unique(popMat$area))
  if(any(areas != easpa$area))
    stop("area names and easpa do not match popMat or are not in the correct order")
  
  # we will need to draw separate multinomial for each stratum. Start by 
  # creating matrix of all draws of |C^g|
  eaSamples = matrix(NA, nrow=nrow(popMat), ncol=n)
  
  # now draw multinomials
  if(stratifyByUrban) {
    # draw for each area crossed with urban/rural
    urbanSamples = do.call("rbind", lapply(1:length(areas), rMyMultinomial, n=n, urban=TRUE, 
                                           stratifyByUrban=stratifyByUrban, popMat=popMat, easpa=easpa))
    ruralSamples = do.call("rbind", lapply(1:length(areas), rMyMultinomial, n=n, urban=FALSE, 
                                           stratifyByUrban=stratifyByUrban, popMat=popMat, easpa=easpa))
    
    # get the indices used to recombine into the full set of draws
    urbanIndices = unlist(sapply(1:length(areas), getSortIndices, urban=TRUE, popMat=popMat, stratifyByUrban=stratifyByUrban))
    ruralIndices = unlist(sapply(1:length(areas), getSortIndices, urban=FALSE, popMat=popMat, stratifyByUrban=stratifyByUrban))
    
    # recombine into eaSamples
    eaSamples[urbanIndices,] = urbanSamples
    eaSamples[ruralIndices,] = ruralSamples
  } else {
    # draw for each area
    stratumSamples = rbind(sapply(1:length(areas), n=n, rMyMultinomial, 
                                  stratifyByUrban=stratifyByUrban, popMat=popMat, easpa=easpa))
    
    # get the indices used to recombine into the full set of draws
    stratumIndices = c(sapply(1:length(areas), getSortIndices, popMat=popMat, stratifyByUrban=stratifyByUrban))
    
    # recombine into eaSamples
    eaSamples[stratumIndices,] = stratumSamples
  }
  
  # return results
  eaSamples
}

#' @describeIn simPopInternal Gives nPixels x n matrix of draws from the stratified multinomial with values 
# corresponding to the number of EAs in each pixel
rStratifiedMultnomialBySubarea = function(n, popMat, easpa, stratifyByUrban=TRUE, poppsub=NULL, 
                                          min1PerSubarea=TRUE, minSample=1) {
  if(is.null(poppsub)) {
    # use popMat to calculate poppsub
    stop("Calculating poppsub with popMat not yet implemented")
    
    # poppsub: a table with the following variables:
    # subarea
    # area
    # popUrb
    # popRur
    # popTotal
  }
  
  # get area names
  areas = sort(unique(popMat$area))
  subareas = sort(unique(popMat$subarea))
  if(any(areas != easpa$area))
    stop("area names and easpa do not match popMat or are not in the correct order")
  
  # we will need to draw separate multinomial for each stratum. Start by 
  # creating matrix of all draws of |C^g|
  eaSamples = matrix(NA, nrow=nrow(popMat), ncol=n)
  
  # create temporary popMat, except with one row for each constituency
  popSubareaMat = popMat[1:length(subareas),]
  popSubareaMat$area = poppsub$area
  popSubareaMat$subarea = poppsub$subarea
  if(stratifyByUrban) {
    popSubareaMat$urban = FALSE
    popSubareaMat = rbind(popSubareaMat, popSubareaMat)
    popSubareaMat$urban[1:length(subareas)] = TRUE
    popSubareaMat$pop[1:length(subareas)] = poppsub$popUrb
    popSubareaMat$pop[(length(subareas) + 1):(2 * length(subareas))] = poppsub$popRur
  } else {
    popSubareaMat$pop = poppsub$popTotal
  }
  # browser()
  # now draw multinomials
  if(stratifyByUrban) {
    # draw for each constituency in each area crossed with urban/rural
    urbanSamplesCon = do.call("rbind", lapply(1:length(areas), rMyMultinomial, n=n, urban=TRUE, 
                                              stratifyByUrban=stratifyByUrban, popMat=popSubareaMat, easpa=easpa, 
                                              min1PerSubarea=min1PerSubarea, method="mult", 
                                              minSample=minSample))
    ruralSamplesCon = do.call("rbind", lapply(1:length(areas), rMyMultinomial, n=n, urban=FALSE, 
                                              stratifyByUrban=stratifyByUrban, popMat=popSubareaMat, easpa=easpa, 
                                              min1PerSubarea=min1PerSubarea, method="mult", 
                                              minSample=minSample))
    
    # get the indices used to recombine into the full set of draws for the subareas
    urbanIndicesCon = unlist(sapply(1:length(areas), getSortIndices, urban=TRUE, popMat=popSubareaMat, stratifyByUrban=stratifyByUrban))
    ruralIndicesCon = unlist(sapply(1:length(areas), getSortIndices, urban=FALSE, popMat=popSubareaMat, stratifyByUrban=stratifyByUrban)) - length(urbanIndicesCon)
    
    # recombine into eaSamples for the subareas
    urbanSamplesCon[urbanIndicesCon,] = urbanSamplesCon
    ruralSamplesCon[ruralIndicesCon,] = ruralSamplesCon
    
    # draw for each pixel crossed with urban/rural
    urbanSamples = do.call("rbind", lapply(1:length(subareas), rMyMultinomialSubarea, n=n, urban=TRUE, 
                                           stratifyByUrban=stratifyByUrban, popMat=popMat, easpsub=urbanSamplesCon))
    ruralSamples = do.call("rbind", lapply(1:length(subareas), rMyMultinomialSubarea, n=n, urban=FALSE, 
                                           stratifyByUrban=stratifyByUrban, popMat=popMat, easpsub=ruralSamplesCon))
    
    # get the indices used to recombine into the full set of draws
    tempPopMat = popMat
    tempPopMat$area = tempPopMat$subarea
    urbanIndices = unlist(sapply(1:length(subareas), getSortIndices, urban=TRUE, popMat=tempPopMat, stratifyByUrban=stratifyByUrban))
    ruralIndices = unlist(sapply(1:length(subareas), getSortIndices, urban=FALSE, popMat=tempPopMat, stratifyByUrban=stratifyByUrban))
    
    # recombine into eaSamples
    eaSamples[urbanIndices,] = urbanSamples
    eaSamples[ruralIndices,] = ruralSamples
  } else {
    # draw for each constituency in each area crossed with urban/rural
    samplesCon = do.call("rbind", lapply(1:length(areas), rMyMultinomial, n=n, urban=TRUE, 
                                         stratifyByUrban=stratifyByUrban, popMat=popSubareaMat, easpa=easpa, 
                                         min1PerSubarea=min1PerSubarea, method="mult", 
                                         minSample=minSample))
    
    # get the indices used to recombine into the full set of draws for the subareas
    indicesCon = unlist(sapply(1:length(areas), getSortIndices, popMat=popSubareaMat, stratifyByUrban=stratifyByUrban))
    
    # recombine into eaSamples for the subareas
    samplesCon[indicesCon,] = samplesCon
    
    # draw for each pixel in each constituency
    stratumSamples = rbind(sapply(1:length(subareas), n=n, rMyMultinomialSubarea, 
                                  stratifyByUrban=stratifyByUrban, popMat=popMat, easpsub=samplesCon))
    
    # get the indices used to recombine into the full set of draws
    tempPopMat = popMat
    tempPopMat$area = tempPopMat$subarea
    stratumIndices = c(sapply(1:length(subareas), getSortIndices, popMat=tempPopMat, stratifyByUrban=stratifyByUrban))
    
    # recombine into eaSamples
    eaSamples[stratumIndices,] = stratumSamples
  }
  
  # return results
  eaSamples
}

#' @describeIn simPopInternal 
rMyMultinomial = function(n, i, stratifyByUrban=TRUE, urban=TRUE, popMat=NULL, easpa=NULL, min1PerSubarea=FALSE, 
                          method=c("mult1", "mult", "indepMH"), minSample=1) {
  method = match.arg(method)
  
  # get area names
  areas = sort(unique(popMat$area))
  if(any(areas != easpa$area))
    stop("area names and easpa do not match popMat or are not in the correct order")
  
  # determine which pixels and how many EAs are in this stratum
  if(stratifyByUrban) {
    includeI = popMat$area == areas[i] & popMat$urban == urban
    nEA = ifelse(urban, easpa$EAUrb[i], easpa$EARur[i])
  }
  else {
    includeI = popMat$area == areas[i]
    nEA = ifelse(urban, easpa$EAUrb[i], easpa$EATotal[i])
  }
  
  # sample from the pixels if this stratum exists
  if(sum(includeI) == 0)
    return(matrix(nrow=0, ncol=n))
  thesePixelProbs = popMat$pop[includeI]
  if(any(thesePixelProbs > 0)) {
    if(!min1PerSubarea) {
      stats::rmultinom(n, nEA, prob=thesePixelProbs)
    } else {
      rmultinom1(n, nEA, prob=thesePixelProbs, method=method, allowSizeLessThanK=TRUE, minSample=minSample)
    }
  } else {
    matrix(0, nrow=length(thesePixelProbs), ncol=n)
  }
}

#' @describeIn simPopInternal 
rMyMultinomialSubarea = function(n, i, easpsub, stratifyByUrban=TRUE, urban=TRUE, popMat=NULL) {
  
  # get constituency names
  subareas = sort(unique(popMat$subarea))
  
  # determine which pixels and how many EAs are in this stratum
  if(stratifyByUrban) {
    includeI = popMat$subarea == subareas[i] & popMat$urban == urban
  }
  else {
    includeI = popMat$subarea == subareas[i]
  }
  nEA = easpsub[i,]
  
  # sample from the pixels if this stratum exists
  if(sum(includeI) == 0){
    if(any(nEA != 0))
      stop(paste0("no valid pixels to put EAs in for constituency ", as.character(subareas[i]), " and urban level ", urban))
    return(matrix(nrow=0, ncol=n))
  }
  thesePixelProbs = popMat$pop[includeI]
  sapply(nEA, stats::rmultinom, n=1, prob=thesePixelProbs)
}

#' @describeIn simPopInternal Random (truncated) multinomial draws conditional on the number of each type being at least one
rmultinom1 = function(n=1, size, prob, maxSize=5000*5000, method=c("mult1", "mult", "indepMH"), verbose=FALSE, minSample=100, 
                      maxExpectedSizeBeforeSwitch=1000*1e7, init=NULL, burnIn=floor(n/4), filterEvery=10, zeroProbZeroSamples=TRUE, 
                      allowSizeLessThanK=FALSE) {
  method = match.arg(method)
  prob = prob*(1/sum(prob))
  
  if(zeroProbZeroSamples && any(prob == 0)) {
    zero = prob == 0
    out = matrix(0, nrow=length(prob), ncol=n)
    
    if(sum(!zero) > 0) {
      out[!zero,] = rmultinom1(n, size, prob[!zero], maxSize, method, verbose, minSample, maxExpectedSizeBeforeSwitch, init, burnIn, 
                               filterEvery, zeroProbZeroSamples, allowSizeLessThanK)
    }
    
    return(out)
  }
  
  k = length(prob)
  if(allowSizeLessThanK && (size <= k)) {
    return(replicate(n, as.numeric(1:k %in% sample(1:k, size, replace=FALSE))))
  } else if(size < k) {
    stop("size < k but rmultinom1 requires at least 1 sample per multinomial type")
  }
  
  maxSamples = floor(maxSize / k)
  averageProbMult = prod((size/k)*prob)
  
  if(method != "indepMH")
    samples = matrix(NA, nrow=k, ncol=n)
  else
    samples = matrix(NA, nrow=k, ncol=round(n*filterEvery))
  if(method == "mult1") {
    averagex = 1 + (size-k)*prob
    averageProb = (size-k) / (prod(averagex))
    
    while(any(is.na(samples))) {
      # calculate the number of remaining samples
      samplesLeft = sum(apply(samples, 2, function(x) {any(is.na(x))}))
      
      # approximate expected number of samples so that, after some are rejected, we will 
      # have the right number of samples
      expectedSamples = ceiling(samplesLeft/averageProb)
      
      if(expectedSamples*k > maxExpectedSizeBeforeSwitch) {
        warning("too many samples expected with method=='mult1'. Switching to method=='indepMH'")
        return(rmultinom1(n, size, prob, maxSize, method="indepMH", verbose, minSample, 
                          maxExpectedSizeBeforeSwitch, init, burnIn, filterEvery, 
                          zeroProbZeroSamples, allowSizeLessThanK))
      }
      
      # sample expectedSamples times a fudge factor, but make sure we don't get past memory limit
      thisNumberOfSamples = max(minSample, min(maxSamples, expectedSamples * 1.1))
      if(verbose)
        print(paste0("Sampling ", thisNumberOfSamples, ". Sampled ", n-samplesLeft, "/", n, ". Expected remaining samples: ", expectedSamples))
      thisSamples = 1 + stats::rmultinom(thisNumberOfSamples, size-k, prob=prob)
      
      # calculate accept probabilities
      thisProbs = (size-k) / apply(thisSamples, 2, prod)
      if(verbose) {
        print(paste0("Max sampled accept prob: ", max(thisProbs), ". Mean sampled accept prob: ", mean(thisProbs)))
        print(paste0("Max theoretical accept prob: ", 1, ". Mean 'theoretical' accept prob: ", averageProb))
      }
      
      # reject relevant samples
      u = stats::runif(thisNumberOfSamples)
      thisSamples = thisSamples[,u<thisProbs]
      
      # remove excess samples if necessary
      totalSamples = ncol(thisSamples) + n - samplesLeft
      if(totalSamples > n) {
        thisSamples = thisSamples[,1:samplesLeft]
      }
      
      # add in accepted samples, if any
      if(ncol(thisSamples) > 0) {
        samples[,(n-samplesLeft+1):(n-samplesLeft+ncol(thisSamples))] = thisSamples
      } else {
        warning(paste0("no samples accepted this round out of ", thisNumberOfSamples, " total..."))
      }
    }
  } else if(method == "mult") {
    
    while(any(is.na(samples))) {
      # calculate the number of remaining samples
      samplesLeft = sum(apply(samples, 2, function(x) {any(is.na(x))}))
      
      # approximate expected number of samples so that, after some are rejected, we will 
      # have the right number of samples
      expectedSamples = ceiling(samplesLeft/averageProbMult)
      
      if(expectedSamples*k > maxExpectedSizeBeforeSwitch) {
        warning("too many samples expected with method=='mult'. Switching to method=='mult1'")
        return(rmultinom1(n, size, prob, maxSize, method="mult1", verbose, minSample, maxExpectedSizeBeforeSwitch, 
                          init, burnIn, filterEvery, 
                          zeroProbZeroSamples, allowSizeLessThanK))
      }
      
      # sample expectedSamples times a fudge factor, but make sure we don't get past memory limit
      thisNumberOfSamples = max(minSample, min(maxSamples, expectedSamples * 1.1))
      if(verbose)
        print(paste0("Sampling ", thisNumberOfSamples, ". Sampled ", n-samplesLeft, "/", n, ". Expected remaining samples: ", expectedSamples))
      thisSamples = matrix(stats::rmultinom(thisNumberOfSamples, size, prob=prob), ncol=thisNumberOfSamples)
      
      # reject relevant samples
      accept = apply(thisSamples, 2, function(x) {all(x>0)})
      thisSamples = matrix(thisSamples[,accept], nrow=length(prob))
      
      # remove excess samples if necessary
      totalSamples = ncol(thisSamples) + n - samplesLeft
      if(totalSamples > n) {
        thisSamples = matrix(thisSamples[,1:samplesLeft], nrow=length(prob))
      }
      
      # add in accepted samples, if any
      if(ncol(thisSamples) > 0) {
        samples[,(n-samplesLeft+1):(n-samplesLeft+ncol(thisSamples))] = thisSamples
      } else {
        warning(paste0("no samples accepted this round out of ", thisNumberOfSamples, " total..."))
      }
    }
  } else if(method == "indepMH") {
    # we use the mult1 method for independent proposals with independent Metropolis-Hastings
    # (https://www.statlect.com/fundamentals-of-statistics/Metropolis-Hastings-algorithm#hid9)
    
    # initialize at something reasonable, if not set by user
    if(is.null(init)) {
      init = 1 + floor(size*prob)
      while(sum(init) > size) {
        tooMuchBy = sum(init) - size
        numberReducible = sum(init > 1)
        reduceNumber = min(tooMuchBy, numberReducible)
        init[init > 1][1:reduceNumber] = init[init > 1][1:reduceNumber] - 1
      }
    }
    if(sum(init) != size)
      stop("sum(init) != size")
    
    # approximate target log-density
    lp <- log(prob)
    lf <- function(x) {
      if(any(x < 1) || sum(x) != size)
        return(-Inf)
      sum(lp*x - lfactorial(x))
    }
    
    # true proposal log-density
    lq = function(x) {
      if(sum(x) != size)
        return(-Inf)
      sum(lp*(x-1) - lfactorial(x-1)) + lfactorial(size-k)
    }
    
    # proposal function
    q <- function(x) {
      1 + stats::rmultinom(1, size-k, prob)
    }
    
    # do the sampling
    tmp <- init
    ar <- 0
    for (i in 1:burnIn) {
      proposal <- q(tmp)
      p <- exp((lf(proposal) - lq(proposal)) - (lf(tmp) - lq(tmp)))
      if (stats::runif(1) < p) {
        tmp <- proposal
      }
    }
    for (i in 1:ncol(samples)) {
      proposal <- q(tmp)
      p <- exp((lf(proposal) - lq(proposal)) - (lf(tmp) - lq(tmp)))
      if (stats::runif(1) < p) {
        tmp <- proposal
        ar <- ar + 1
      }
      samples[,i] <- tmp
    }
    
    # calculated acceptance percentage
    if(verbose) {
      print(paste0("acceptance percentage: ", ar/ncol(samples)))
    }
    
    # filter out samples to reduce autocorrelation
    samples = samples[,seq(from=1, to=ncol(samples), by=filterEvery)]
  }
  
  samples
}

#' @describeIn simPopInternal Take multilevel multinomial draws first from joint distribution of 
#' number of households per EA given the total per stratum, and then from the joint 
#' distribution of the total target population per household given 
#' the total per stratum
sampleNMultilevelMultinomial = function(nDraws = ncol(pixelIndexMat), pixelIndexMat=NULL, urbanMat=NULL, areaMat=NULL, easpaList, 
                                        popMat, stratifyByUrban=TRUE, verbose=TRUE, returnEAinfo=FALSE) {
  
  if((is.null(areaMat) || is.null(urbanMat)) && is.null(pixelIndexMat)) {
    stop("user must either supply pixelIndexMat or both areaMat and urbanMat")
  }
  if(is.null(areaMat)) {
    areaMat = matrix(popMat$area[pixelIndexMat], ncol=nDraws)
  }
  if(is.null(urbanMat)) {
    urbanMat = matrix(popMat$urban[pixelIndexMat], ncol=nDraws)
  }
  
  # start by drawing the totals, then divide households amongst EAs, then divide target population amongst households. 
  # Make sure there are at least 25 households per EA (only divide the rest randomly)
  
  ##### Draw the totals
  
  # get the total number of enumeration areas per stratum (this does not change between draws)
  areas = easpaList[[1]]$area
  totalEAsUrban = easpaList[[1]]$EAUrb
  totalEAsRural = easpaList[[1]]$EARur
  totalEAs = easpaList[[1]]$EATotal
  nEAs = sum(totalEAs)
  
  ##### draw the total target population per enumeration area
  
  targetPopDraws = matrix(nrow=nEAs, ncol=nDraws)
  if(returnEAinfo) {
    householdDraws = matrix(nrow=nEAs, ncol=nDraws)
  }
  
  # Draw the number of households per stratum area that will be randomly distributed (total minus the minimum 25)
  if(stratifyByUrban) {
    totalHouseholdsUrban = sweep(sapply(easpaList, function(x) {x$HHUrb}), 1, -25*totalEAsUrban, "+")
    totalHouseholdsRural = sweep(sapply(easpaList, function(x) {x$HHRur}), 1, -25*totalEAsRural, "+")
    totalChildrenUrban = sapply(easpaList, function(x) {x$popUrb})
    totalChildrenRural = sapply(easpaList, function(x) {x$popRur})
  } else {
    totalHouseholds = sweep(sapply(easpaList, function(x) {x$HHTotal}), 1, -25*totalEAs, "+")
    totalChildren = sapply(easpaList, function(x) {x$popHHTotal})
  }
  
  # distribute the households throughout the enumeration areas with multinomial distribution, then 
  # distribute the target population amongst the households, also with a multinomial distribution
  for(i in 1:length(areas)) {
    thisArea = areas[i]
    
    # print progress if in verbose mode
    if(verbose) {
      print(paste0("drawing Ns for each EA for area ", thisArea, " (", i, "/", length(areas), ")"))
    }
    
    # draw households per EA (make sure there are any rural EAs)
    if(stratifyByUrban) {
      householdDrawsUrban <- matrix(sapply(totalHouseholdsUrban[i,], stats::rmultinom, n=1, prob=rep(1/totalEAsUrban[i], totalEAsUrban[i])), nrow=totalEAsUrban[i], ncol=nDraws) + 25
      if(totalEAsRural[i] != 0) {
        householdDrawsRural <- matrix(sapply(totalHouseholdsRural[i,], stats::rmultinom, n=1, prob=rep(1/totalEAsRural[i], totalEAsRural[i])), nrow=totalEAsRural[i], ncol=nDraws) + 25
      }
      
      # if we must return EA info, we must return the household draws for each EA:
      if(returnEAinfo) {
        if(totalEAsUrban[i] != 0) {
          householdDraws[areaMat==thisArea & urbanMat] = householdDrawsUrban
        }
        if(totalEAsRural[i] != 0) {
          householdDraws[areaMat==thisArea & !urbanMat] = householdDrawsRural
        }
      }
    } else {
      householdDraws = matrix(sapply(totalHouseholds[i,], stats::rmultinom, n=1, prob=rep(1/totalEAs[i], totalEAs[i])), nrow=totalEAs[i], ncol=nDraws) + 25
    }
    
    # drawing target population per EA, with probability proportional to the number of households
    if(stratifyByUrban) {
      if(totalEAsUrban[i] != 0) {
        probsUrban = sweep(householdDrawsUrban, 2, 1 / colSums(householdDrawsUrban), "*")
        targetPopDraws[areaMat==thisArea & urbanMat] = sapply(1:nDraws, function(j) {stats::rmultinom(1, totalChildrenUrban[i,j], probsUrban[,j])})
      }
      
      if(totalEAsRural[i] != 0) {
        probsRural = sweep(householdDrawsRural, 2, 1 / colSums(householdDrawsRural), "*")
        targetPopDraws[areaMat==thisArea & !urbanMat] = sapply(1:nDraws, function(j) {stats::rmultinom(1, totalChildrenRural[i,j], probsRural[,j])})
      }
    } else {
      probs = sweep(householdDraws, 2, 1 / colSums(householdDraws), "*")
      targetPopDraws[areaMat==thisArea] = sapply(1:nDraws, function(j) {stats::rmultinom(1, totalChildren[i,j], probs[,j])})
    }
  }
  
  ##### Return results
  if(!returnEAinfo) {
    targetPopDraws
  } else {
    list(householdDraws=householdDraws, targetPopDraws=targetPopDraws)
  }
}

#' @describeIn simPopInternal Same as sampleNMultilevelMultinomial, except the number of EAs per pixel is fixed
sampleNMultilevelMultinomialFixed = function(clustersPerPixel, nDraws=ncol(pixelIndices), pixelIndices=NULL, 
                                             urbanVals=NULL, areaVals=NULL, easpa, popMat, stratifyByUrban=TRUE, 
                                             verbose=TRUE) {
  
  # set default inputs
  if((is.null(areaVals) || is.null(urbanVals)) && is.null(pixelIndices)) {
    stop("user must either supply pixelIndices or both areaVals and urbanVals")
  }
  if(is.null(areaVals)) {
    areaVals = matrix(popMat$area[pixelIndices], ncol=nDraws)
  }
  if(is.null(urbanVals)) {
    urbanVals = matrix(popMat$urban[pixelIndices], ncol=nDraws)
  }
  
  # start by drawing the totals, then divide households amongst EAs, then divide target population amongst households. 
  # Make sure there are at least 25 households per EA (only divide the rest randomly)
  
  ##### Draw the totals
  
  # get the total number of enumeration areas per stratum (this does not change between draws)
  areas = easpa$area
  totalEAsUrban = easpa$EAUrb
  totalEAsRural = easpa$EARur
  totalEAs = easpa$EATotal
  nEAs = sum(totalEAs)
  
  if(nEAs != sum(clustersPerPixel)) {
    stop("sum(easpa$EATotal) != sum(clustersPerPixel)")
  }
  
  ##### draw the total target population per EA
  targetPopDraws = matrix(nrow=nEAs, ncol=nDraws)
  
  # Draw the number of households per stratum area that will be randomly distributed (total minus the minimum 25)
  if(stratifyByUrban) {
    totalHouseholdsUrban = easpa$HHUrb -25*totalEAsUrban
    totalHouseholdsRural = easpa$HHRur -25*totalEAsRural
    totalChildrenUrban = easpa$popUrb
    totalChildrenRural = easpa$popRur
  } else {
    totalHouseholds = easpa$HHTotal - 25*totalEAs
    totalChildren = easpa$popHHTotal
  }
  
  # distribute the households throughout the enumeration areas with multinomial distribution, then 
  # distribute the target population amongst the households, also with a multinomial distribution
  for(i in 1:length(areas)) {
    thisArea = areas[i]
    
    # print progress if in verbose mode
    if(verbose) {
      print(paste0("drawing Ns for each EA for area ", thisArea, " (", i, "/", length(areas), ")"))
    }
    
    # draw households per EA (make sure there are any rural EAs)
    if(stratifyByUrban) {
      if(totalEAsUrban[i] != 0) {
        if(any(totalHouseholdsUrban != 0)) {
          householdDrawsUrban = stats::rmultinom(n=nDraws, size=totalHouseholdsUrban[i], prob=rep(1/totalEAsUrban[i], totalEAsUrban[i])) + 25
        } else {
          householdDrawsUrban = matrix(rep(25, totalEAsUrban[i]*nDraws), ncol=nDraws)
        }
      }
      
      if(totalEAsRural[i] != 0) {
        if(any(totalHouseholdsRural != 0)) {
          householdDrawsRural = stats::rmultinom(n=nDraws, size=totalHouseholdsRural[i], prob=rep(1/totalEAsRural[i], totalEAsRural[i])) + 25
        } else {
          householdDrawsRural = matrix(rep(25, totalEAsRural[i]*nDraws), ncol=nDraws)
        }
      }
    } else {
      if(any(totalHouseholdsRural != 0)) {
        householdDraws = stats::rmultinom(n=nDraws, size=totalHouseholds[i], prob=rep(1/totalEAs[i], totalEAs[i])) + 25
      } else {
        householdDraws = matrix(rep(25, totalEAs[i]*nDraws), ncol=nDraws)
      }
    }
    
    # drawing target population per EA, with probability proportional to the number of households
    if(stratifyByUrban) {
      
      if(totalEAsUrban[i] != 0) {
        probsUrban = sweep(householdDrawsUrban, 2, 1 / colSums(householdDrawsUrban), "*")
        targetPopDraws[areaVals==thisArea & urbanVals] = sapply(1:nDraws, function(j) {stats::rmultinom(1, totalChildrenUrban[i], probsUrban[,j])})
      }
      
      if(totalEAsRural[i] != 0) {
        probsRural = sweep(householdDrawsRural, 2, 1 / colSums(householdDrawsRural), "*")
        targetPopDraws[areaVals==thisArea & !urbanVals] = sapply(1:nDraws, function(j) {stats::rmultinom(1, totalChildrenRural[i], probsRural[,j])})
      }
      
    } else {
      probs = sweep(householdDraws, 2, 1 / colSums(householdDraws), "*")
      targetPopDraws[areaVals==thisArea] = sapply(1:nDraws, function(j) {stats::rmultinom(1, totalChildren[i], probs[,j])})
    }
  }
  
  ##### Return results
  targetPopDraws
}