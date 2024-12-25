# Functions for simulating populations

#' Simulate populations and areal prevalences
#' 
#' Given a spatial risk model, simulate populations and population prevalences at the 
#' enumeration area level (represented as points), and aggregate to the pixel and 
#' administrative areal level.
#' 
#' 
#' `r lifecycle::badge("experimental")`
#' 
#' @param nsim Number of simulations
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
#' @param pop.mat Pixellated grid data frame with variables `lon`, `lat`, `pop`, `area`, `subareas` (if subarea.level is TRUE), `urban` (if stratify.by.urban is TRUE), `east`, and `north`
#' @param target.pop.mat Same as pop.mat, but `pop` variable gives target rather than general population
#' @param poppsub data.frame of population per subarea separated by 
#' urban/rural using for population density grid normalization or urbanicity 
#' classification. Often based on extra fine scale population density grid. Has variables:
#' @param spde.mesh Triangular mesh for the SPDE
#' @param marg.var Marginal variance of the spatial process, excluding cluster effects. 
#'          If 0, no spatial component is included
#' @param eff.range Effective spatial range for the SPDE model
#' @param beta0 Intercept of logit model for risk
#' @param gamma Effect of urban on logit scale for logit model for risk
#' @param sigma.epsilon Standard deviation on the logit scale for iid Gaussian EA level random effects in the risk model
#' @param seed Random number generator seed
#' @param inla.seed Seed input to inla.qsample. 0L sets seed intelligently, 
#'            > 0 sets a specific seed, < 0 keeps existing RNG
#' @param n.HH.sampled Number of households sampled per enumeration area. Default is 25 to match DHS surveys
#' @param stratify.by.urban Whether or not to stratify simulations by urban/rural classification
#' @param subarea.level Whether or not to aggregate the population by subarea
#' @param do.fine.scale.risk Whether or not to calculate the fine scale risk at each aggregation level in addition to the prevalence
#' @param do.smooth.risk Whether or not to calculate the smooth risk at each aggregation level in addition to the prevalence
#' @param do.smooth.risk.logistic.approx Whether to use logistic approximation when calculating smooth risk. See 
#' \code{\link{logitNormMean}} for details.
#' @param min1.per.subarea If TRUE, ensures there is at least 1 EA per subarea. If subareas are particularly unlikely to 
#' have enumeration areas since they have a very low proportion of the population in an area, then setting this to TRUE may be 
#' computationally intensive.
#' @param logit.risk.draws nIntegrationPoints x nsim dimension matrix of draws from the pixel leve risk field on logit scale, leaving out 
#' potential nugget/cluster/EA level effects.
#' @param sigma.epsilon.draws nsim length vector of draws of cluster effect logit scale SD (joint draws with logit.risk.draws)
#' @param validation.pixel.I CURRENTLY FOR TESTING PURPOSES ONLY a set of indices of pixels for which we want to simulate populations (used for pixel level validation)
#' @param validation.cluster.I CURRENTLY FOR TESTING PURPOSES ONLY a set of indices of cluster for which we want to simulate populations (used for cluster level validation)
#' @param clusters.per.pixel CURRENTLY FOR TESTING PURPOSES ONLY Used for pixel level validation. Fixes the number of EAs per pixel.
#' @param return.EA.info If TRUE, returns information on every individual EA (BAU) for each simulated population
#' @param epsc nEAs x nsim matrix of simulated EA (BAU) level iid effects representing fine scale variation in 
#'       risk. If NULL, they are simulated as iid Gaussian on a logit scale with 
#'       SD given by sigma.epsilon.draws
#' list(pixelPop=outPixelLevel, subareaPop=outSubareaLevel, areaPop=outAreaLevel, logit.risk.draws=logit.risk.draws)
#' @return The simulated population aggregated to the enumeration area, 
#' pixel, subarea (generally Admin2), and area (generally Admin1) levels. Output includes:
#' \item{pixelPop}{A list of pixel level population aggregates}
#' \item{subareaPop}{A list of `subarea` level population aggregates}
#' \item{areaPop}{A list of `area` level population aggregates}
#' Each of these contains population numerator and denominator as well as prevalence and risk 
#' information aggregated to the appropriate level.
#' 
#' @details For population simulation and aggregation, we consider three models: smooth  
#' risk, fine scale risk, and the fine scale prevalence. All will be described in detail 
#' in a paper in preparation. In the smooth risk model, pixel level risks are integrated 
#' with respect to target population density when producing areal estimates on a prespecified 
#' set of integration points. The target population may be, for example, neonatals rather 
#' than the general population. In the fine scale models, enumeration areas (EAs) are simulated as 
#' point locations and iid random effects in the EA level risk are allowed. EAs and populations are dispersed conditional on the (possibly 
#' approximately) known number of EAs, households, and target population at a particular 
#' areal level (these we call `areas`) using multilevel multinomial sampling, first 
#' sampling the EAs, then distributing households among the EAs, then the target population 
#' among the households. Any areal level below the `areas` we call `subareas`. For instance, 
#' the `areas` might be Admin-1 if that is the smallest level at which the number of EAs, 
#' households, and people is known, and the `subareas` might be Admin-2. The multilevel 
#' multinomial sampling may be stratified by urban/rural within the areas if the number of 
#' EAs, households, and people is also approximately known at that level.
#' 
#' Within each EA we assume a fixed probability of an event occurring, which is the fine scale `risk`. 
#' The fine scale `prevalence` is the empirical proportion of events within that EA. We assume EA 
#' level logit scale iid N(0, sigma.epsilon^2) random effects in the risk model. When averaged 
#' with equal weights over all EAs in an areal unit, this forms the fine scale risk. When 
#' instead the population numerators and denominators are aggregated, and are used to 
#' calculate the empirical proportion of events occurring in an areal unit, the resulting 
#' quantity is the fine scale prevalence in that areal unit.
#' 
#' Note that these functions can be used for either simulating populations for simulation 
#' studies, or for generating predictions accounting for uncertainty in EA locations 
#' and fine scale variation occurring at the EA level due to EA level iid random effects. 
#' Required, however, is a separately fit EA level spatial risk model 
#' and information on the spatial population density and the population frame.
#' 
#' @author John Paige
#' @references In preparation
#' @seealso \code{\link{simPopCustom}}, \code{\link{makePopIntegrationTab}}, \code{\link{adjustPopMat}}, \code{\link{simSPDE}}.
#' @examples 
#' \dontrun{
#' ## In this script we will create 5km resolution pixellated grid over Kenya, 
#' ## and generate tables of estimated (both target and general) population 
#' ## totals at the area (e.g. Admin-1) and subarea (e.g. Admin-2) levels. Then 
#' ## we will use that to simulate populations of 
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
#' out = load(mapsFilename)
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
#' system.time(pop.matKenya <- makePopIntegrationTab(
#'   kmRes=5, pop=pop, domainMapDat=adm0,
#'   eastLim=eastLim, northLim=northLim, mapProjection=projKenya,
#'   poppa = poppaKenya, poppsub=poppsubKenya, 
#'   areaMapDat = adm1, subareaMapDat = adm2,
#'   areaNameVar = "NAME_1", subareaNameVar="NAME_2"))
#' 
#' ## Adjust pop.mat to be target (neonatal) rather than general population 
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
#' pop.matKenyaNeonatal = adjustPopMat(pop.matKenya, easpaKenyaNeonatal)
#' 
#' # Generate neonatal population table from the neonatal population integration 
#' # matrix. This is technically not necessary for population simulation purposes, 
#' # but is here for illustrative purposes
#' poppsubKenyaNeonatal = poppRegionFromPopMat(pop.matKenyaNeonatal, 
#'                                             pop.matKenyaNeonatal$subarea)
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
#' eff.range = 400 # effective spatial range in km
#' sigma.epsilon=sqrt(1/2.5) # cluster (nugget) effect standard deviation
#' 
#' # Run a simulation! This produces multiple dense nEA x nsim and nPixel x nsim 
#' # matrices. In the future sparse matrices and chunk by chunk computations 
#' # may be incorporated.
#' simPop = simPopSPDE(nsim=1, easpa=easpaKenyaNeonatal, 
#'                     pop.mat=pop.matKenya, target.pop.mat=pop.matKenyaNeonatal, 
#'                     poppsub=poppsubKenya, spde.mesh=kenyaMesh, 
#'                     marg.var=rho, sigma.epsilon=sigma.epsilon, 
#'                     gamma=gamma, eff.range=eff.range, beta0=beta0, 
#'                     seed=12, inla.seed=12, n.HH.sampled=25, 
#'                     stratify.by.urban=TRUE, subarea.level=TRUE, 
#'                     do.fine.scale.risk=TRUE, do.smooth.risk=TRUE, 
#'                     min1.per.subarea=TRUE)
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
#' cbind(aggregate(simPop$eaPop$ea.samples[,1], by=list(area=pop.matKenya$area), FUN=sum), 
#'       trueNumEAs=easpaKenya$EATotal[order(easpaKenya$area)])
#' aggregate(simPop$eaPop$ea.samples[,1], by=list(area=pop.matKenya$subarea), FUN=sum)
#' 
#' ## plot simulated population
#' # directory for plotting 
#' # (mapPlot takes longer when it doesn't save to a file)
#' tempDirectory = "~/"
#' 
#' # pixel level plots. Some pixels have no simulated EAs, in which case they will be 
#' # plotted as white. Expected noisy looking plots of fine scale risk and prevalence 
#' # due to EAs being discrete, as compared to a very smooth plot of smooth risk.
#' zlim = c(0, quantile(probs=.995, c(simPop$pixelPop$pFineScalePrevalence, 
#'                                    simPop$pixelPop$pFineScaleRisk, 
#'                                    simPop$pixelPop$pSmoothRisk), na.rm=TRUE))
#' pdf(file=paste0(tempDirectory, "simPopSPDEPixel.pdf"), width=8, height=8)
#' par(mfrow=c(2,2))
#' plot(adm2, asp=1)
#' points(simPop$eaPop$eaDatList[[1]]$lon, simPop$eaPop$eaDatList[[1]]$lat, pch=".", col="blue")
#' plot(adm2, asp=1)
#' quilt.plot(pop.matKenya$lon, pop.matKenya$lat, simPop$pixelPop$pFineScalePrevalence, 
#'            zlim=zlim, add=TRUE, FUN=function(x) {mean(x, na.rm=TRUE)})
#' plot(adm2, asp=1)
#' quilt.plot(pop.matKenya$lon, pop.matKenya$lat, simPop$pixelPop$pFineScaleRisk, 
#'            zlim=zlim, add=TRUE, FUN=function(x) {mean(x, na.rm=TRUE)})
#' quilt.plot(pop.matKenya$lon, pop.matKenya$lat, simPop$pixelPop$pSmoothRisk, 
#'            zlim=zlim, FUN=function(x) {mean(x, na.rm=TRUE)}, asp=1)
#' plot(adm2, add=TRUE)
#' dev.off()
#' 
#' range(simPop$eaPop$eaDatList[[1]]$N)
#' 
#' # areal (Admin-1) level: these results should look essentially identical
#' 
#' tempDat = simPop$areaPop$aggregationResults[c("region", "pFineScalePrevalence", 
#'                                                "pFineScaleRisk", "pSmoothRisk")]
#' pdf(file=paste0(tempDirectory, "simPopSPDEAdmin-1.pdf"), width=7, height=7)
#' mapPlot(tempDat, 
#'         variables=c("pFineScalePrevalence", "pFineScaleRisk", "pSmoothRisk"), 
#'         geo=adm1, by.geo="NAME_1", by.data="region", is.long=FALSE)
#' dev.off()
#' 
#' # subareal (Admin-2) level: these results should look subtley different 
#' # depending on the type of prevalence/risk considered
#' tempDat = simPop$subareaPop$aggregationResults[c("region", "pFineScalePrevalence", 
#'                                                   "pFineScaleRisk", "pSmoothRisk")]
#' pdf(file=paste0(tempDirectory, "simPopSPDEAdmin-2.pdf"), width=7, height=7)
#' mapPlot(tempDat, 
#'         variables=c("pFineScalePrevalence", "pFineScaleRisk", "pSmoothRisk"), 
#'         geo=adm2, by.geo="NAME_2", by.data="region", is.long=FALSE)
#' dev.off()
#' }
#' @name simPop
NULL
#' @importFrom stats cor
#' @describeIn simPop
#' Simulate populations and population prevalences given census frame and population density 
#' information. Uses SPDE model for generating spatial risk and can include iid cluster 
#' level effect.
#' 
#' @export
simPopSPDE = function(nsim=1, easpa, pop.mat, target.pop.mat, poppsub, spde.mesh, 
                      marg.var=0.243, sigma.epsilon=sqrt(0.463), 
                      gamma=0.009, eff.range=406.51, beta0=-3.922, 
                      seed=NULL, inla.seed=-1L, n.HH.sampled=25, 
                      stratify.by.urban=TRUE, subarea.level=TRUE, 
                      do.fine.scale.risk=FALSE, do.smooth.risk=FALSE, 
                      do.smooth.risk.logistic.approx=FALSE, 
                      min1.per.subarea=TRUE) {
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
  pixelCoords = cbind(pop.mat$east, pop.mat$north)
  
  if(marg.var != 0) {
    SPDEArgs = list(coords=pixelCoords, nsim=nsim, marg.var=marg.var, eff.range=eff.range, 
                    mesh=spde.mesh, inla.seed=inla.seed)
    simVals = do.call("simSPDE", SPDEArgs)
  } else {
    simVals = matrix(rep(0, nrow(pixelCoords)), ncol=1)
  }
  
  # add in intercept
  simVals = simVals + beta0
  
  # add in urban effect
  simVals = sweep(simVals, 1, gamma*pop.mat$urban, "+")
  
  # simulate nugget/cluster effect
  epsc = matrix(stats::rnorm(totalEAs*nsim, sd=sigma.epsilon), ncol=nsim)
  
  # transform back to original scale for the pixel level probabilities
  probsNoNug = expit(simVals)
  
  # simulate the enumeration areas
  logit.risk.draws = simVals
  sigma.epsilon.draws = rep(sigma.epsilon, nsim)
  
  print("Using SPDE model to simulate EA and pixel level populations")
  outPixelLevel = simPopCustom(logit.risk.draws=logit.risk.draws, sigma.epsilon.draws=sigma.epsilon.draws, easpa=easpa, 
                              pop.mat=pop.mat, target.pop.mat=target.pop.mat, 
                              stratify.by.urban=stratify.by.urban, do.smooth.risk=do.smooth.risk, 
                              do.smooth.risk.logistic.approx=do.smooth.risk.logistic.approx, 
                              do.fine.scale.risk=do.fine.scale.risk, poppsub=poppsub, 
                              subarea.level=subarea.level, min1.per.subarea=min1.per.subarea, 
                              return.EA.info=TRUE, epsc=epsc)
  eaPop = list(eaDatList=outPixelLevel$eaDatList, ea.samples=outPixelLevel$ea.samples)
  outPixelLevel$eaDatList = NULL
  outPixelLevel$ea.samples = NULL
  
  if(subarea.level) {
    print("aggregating from pixel level to subarea level")
    outSubareaLevel = pixelPopToArea(pixel.level.pop=outPixelLevel, ea.samples=eaPop$ea.samples, 
                                     areas=pop.mat$subarea, stratify.by.urban=stratify.by.urban, 
                                     target.pop.mat=target.pop.mat, do.fine.scale.risk=do.fine.scale.risk, 
                                     do.smooth.risk=do.smooth.risk)
    
    print("aggregating from subarea level to area level")
    
    # get areas associated with each subarea for aggregation
    tempAreasFrom = pop.mat$subarea
    tempAreasTo = pop.mat$area
    areas.from = sort(unique(tempAreasFrom))
    areas.toI = match(areas.from, tempAreasFrom)
    areas.to = tempAreasTo[areas.toI]
    
    # do the aggregation from subareas to areas
    outAreaLevel = areaPopToArea(area.level.pop=outSubareaLevel, areas.from=areas.from, areas.to=areas.to, 
                                 stratify.by.urban=stratify.by.urban, do.fine.scale.risk=do.fine.scale.risk, 
                                 do.smooth.risk=do.smooth.risk)
  } else {
    outSubareaLevel = NULL
    
    print("aggregating from pixel level to area level")
    outAreaLevel = pixelPopToArea(pixel.level.pop=outPixelLevel, ea.samples=eaPop$ea.samples, 
                                  areas=pop.mat$area, stratify.by.urban=stratify.by.urban, 
                                  do.fine.scale.risk=do.fine.scale.risk, do.smooth.risk=do.smooth.risk)
  }
  
  list(eaPop=eaPop, pixelPop=outPixelLevel, subareaPop=outSubareaLevel, areaPop=outAreaLevel, logit.risk.draws=logit.risk.draws)
}

#' Aggregate populations to the specified areal level
#' 
#' Takes simulated populations and aggregates 
#' them to the specified areal level. Also calculates the aggregated risk and prevalence.
#' 
#' 
#' `r lifecycle::badge("experimental")`
#' 
#' @param pixel.level.pop pixel level population information that we want aggregate. In the same format as output from \code{\link{simPopCustom}}
#' @param ea.samples nIntegrationPoint x nsim matrix of the number of enumeration areas per pixel sampled in the input pixel level population
#' @param areas character vector of length nIntegrationPoints of area names over which we 
#'        want to aggregate. Can also be subareas
#' @param stratify.by.urban whether or not to stratify simulations by urban/rural classification
#' @param target.pop.mat pixellated grid data frame with variables `lon`, `lat`, `pop` (target population), `area`, `subareas` (if subarea.level is TRUE), `urban` (if stratify.by.urban is TRUE), `east`, and `north`
#' @param do.fine.scale.risk whether or not to calculate the fine scale risk in addition to the prevalence. See details
#' @param do.smooth.risk Whether or not to calculate the smooth risk in addition to the prevalence. See details
#' @param area.level.pop output of \code{\link{simPopCustom}} containing pixel level information 
#'               about the population of interest
#' @param areas.from character vector of length equal to the number of areas from which 
#'            we would like to aggregate containing the unique names of the areas. 
#'            Can also be subareas, but these are smaller than the "to areas", and 
#'            each "from area" must be entirely contained in a single "to area"
#' @param areas.to character vector of length equal to the number of areas from which 
#'          we would like to aggregate containing the names of the areas containing 
#'          with each respective `from' area. Can also be a set of subareas, 
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
#' eff.range = 400 # effective spatial range in km
#' sigma.epsilon=sqrt(1/2.5) # cluster (nugget) effect standard deviation
#' 
#' # simulate the population! Note that this produces multiple dense 
#' # nEA x nsim and nIntegrationPoint x nsim matrices. In the future 
#' # sparse matrices will and chunk by chunk computations may be incorporated.
#' simPop = simPopSPDE(nsim=1, easpa=easpaKenyaNeonatal, 
#'                     pop.mat=pop.matKenya, target.pop.mat=pop.matKenyaNeonatal, 
#'                     poppsub=poppsubKenya, spde.mesh=kenyaMesh, 
#'                     marg.var=rho, sigma.epsilonSq=sigma.epsilon^2, 
#'                     gamma=gamma, eff.range=eff.range, beta0=beta0, 
#'                     seed=123, inla.seed=12, n.HH.sampled=25, 
#'                     stratify.by.urban=TRUE, subarea.level=TRUE, 
#'                     do.fine.scale.risk=TRUE, 
#'                     min1.per.subarea=TRUE)
#' 
#' pixelPop = simPop$pixelPop
#' subareaPop = pixelPopToArea(pixel.level.pop=pixelPop, ea.samples=pixelPop$ea.samples, 
#'   areas=pop.matKenya$subarea, stratify.by.urban=TRUE, 
#'   target.pop.mat=pop.matKenyaNeonatal, do.fine.scale.risk=TRUE)
#' 
#' # get areas associated with each subarea for aggregation
#' tempAreasFrom = pop.matKenya$subarea
#' tempAreasTo = pop.matKenya$area
#' areas.from = sort(unique(tempAreasFrom))
#' areas.toI = match(areas.from, tempAreasFrom)
#' areas.to = tempAreasTo[areas.toI]
#' 
#' # do the aggregation from subareas to areas
#' outAreaLevel = areaPopToArea(area.level.pop=subareaPop, 
#'   areas.from=areas.from, areas.to=areas.to, 
#'   stratify.by.urban=TRUE, do.fine.scale.risk=TRUE)
#' }
NULL

#' @describeIn aggPop Aggregate from pixel to areal level
#' @export
pixelPopToArea = function(pixel.level.pop, ea.samples, areas, stratify.by.urban=TRUE, target.pop.mat=NULL, 
                          do.fine.scale.risk=!is.null(pixel.level.pop$fineScaleRisk$p), 
                          do.smooth.risk=!is.null(pixel.level.pop$smoothRisk$p)) {
  
  # fine scale prevalence aggregation model
  nSamples = pixel.level.pop$NFineScalePrevalence
  zSamples = pixel.level.pop$ZFineScalePrevalence
  zSamples[is.na(zSamples)] = 0 # must set to zero temporarily so matrix multiplication works
  out = aggPixelPreds(Zg=zSamples, Ng=nSamples, areas=areas, target.pop.mat=target.pop.mat, 
                      use.density=FALSE, stratify.by.urban=stratify.by.urban, normalize=FALSE)
  aggregationResults = out$aggregationResults
  aggregationMatrices = out$aggregationMatrices
  names(aggregationResults)[-1] = paste(names(aggregationResults)[-1], "FineScalePrevalence", sep="")
  names(aggregationMatrices) = paste(names(aggregationMatrices), "FineScalePrevalence", sep="")
  
  # fine scale risk aggregation model
  if(do.fine.scale.risk) {
    nSamplesFineScaleRisk = pixel.level.pop$NFineScaleRisk
    zSamplesFineScaleRisk = pixel.level.pop$ZFineScaleRisk
    zSamplesFineScaleRisk[is.na(zSamplesFineScaleRisk)] = 0 # must set to zero temporarily so matrix multiplication works out
    out = aggPixelPreds(Zg=zSamplesFineScaleRisk, Ng=nSamplesFineScaleRisk, areas=areas, target.pop.mat=target.pop.mat, 
                        use.density=FALSE, stratify.by.urban=stratify.by.urban, normalize=FALSE)
    aggregationResultsFineScaleRisk = out$aggregationResults
    aggregationMatricesFineScaleRisk = out$aggregationMatrices
    names(aggregationResultsFineScaleRisk)[-1] = paste(names(aggregationResultsFineScaleRisk)[-1], "FineScaleRisk", sep="")
    names(aggregationMatricesFineScaleRisk) = paste(names(aggregationMatricesFineScaleRisk), "FineScaleRisk", sep="")
    
    # aggregationResults = merge(aggregationResults, aggregationResultsFineScaleRisk, by="region")
    aggregationResults = c(aggregationResults, aggregationResultsFineScaleRisk[-1])
    aggregationMatrices = c(aggregationMatrices, aggregationMatricesFineScaleRisk)
  }
  
  if(do.smooth.risk) {
    # NOTE: although use.density is set to FALSE, that's only because the density is already 
    # directly incorporated into nSamplesSmoothRisk
    nSamplesSmoothRisk = pixel.level.pop$NSmoothRisk
    zSamplesSmoothRisk = pixel.level.pop$ZSmoothRisk
    zSamplesSmoothRisk[is.na(zSamplesSmoothRisk)] = 0 # must set to zero temporarily so matrix multiplication works out
    out = aggPixelPreds(Zg=zSamplesSmoothRisk, Ng=nSamplesSmoothRisk, areas=areas, target.pop.mat=target.pop.mat, 
                        use.density=FALSE, stratify.by.urban=stratify.by.urban, normalize=FALSE)
    aggregationResultsSmoothRisk = out$aggregationResults
    aggregationMatricesSmoothRisk = out$aggregationMatrices
    names(aggregationResultsSmoothRisk)[-1] = paste(names(aggregationResultsSmoothRisk)[-1], "SmoothRisk", sep="")
    names(aggregationMatricesSmoothRisk) = paste(names(aggregationMatricesSmoothRisk), "SmoothRisk", sep="")
    
    # aggregationResults = merge(aggregationResults, aggregationResultsSmoothRisk, by="region")
    aggregationResults = c(aggregationResults, aggregationResultsSmoothRisk[-1])
    aggregationMatrices = c(aggregationMatrices, aggregationMatricesSmoothRisk)
  }
  
  list(aggregationResults=aggregationResults, aggregationMatrices=aggregationMatrices)
}

#' Helper function of \code{\link{pixelPopToArea}}
#' 
#' Aggregates population from the 
#' pixel level to the level of the area of interest.
#' 
#' @param Zg nIntegrationPoint x nsim matrix of simulated response (population numerators) for each pixel and sample
#' @param Ng nIntegrationPoint x nsim matrix of simulated counts (population denominators) for each pixel and sample
#' @param areas nIntegrationPoint length character vector of areas (or subareas) 
#' @param urban nIntegrationPoint length vector of indicators specifying whether or not pixels are urban or rural
#' @param target.pop.mat same as in \code{\link{simPopCustom}}
#' @param use.density whether to use population density as aggregation weights. 
#' @param stratify.by.urban whether or not to stratify simulations by urban/rural classification
#' @param normalize if TRUE, pixel level aggregation weights within specified area are normalized to sum to 1. This produces an 
#' average of the values in Zg rather than a sum. In general, should only be set to TRUE for smooth integrals of risk.
aggPixelPreds = function(Zg, Ng, areas, urban=target.pop.mat$urban, target.pop.mat=NULL, use.density=FALSE, 
                         stratify.by.urban=TRUE, normalize=use.density) {
  
  if(use.density && !normalize) {
    stop("if use.density is set to TRUE, normalize must be set to TRUE as well")
  }
  predsUrban = urban
  predsArea = areas
  
  # set NAs and pixels without any sample size to 0
  Ng[is.na(Ng)] = 0
  if(!use.density) {
    # is use.density is true, then Zg is really a set of probabilities, so no need to set to 0
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
    
    if(use.density) {
      popDensities = target.pop.mat$pop
      densities = popDensities
    } else {
      equalDensities = rep(1, nrow(Zg))
      densities = equalDensities
    }
    
    uniqueNames = sort(unique(areaNames))
    getMatrixHelper = function(i, thisUrban=NULL, thisUseDensity=use.density, thisNormalize=normalize) {
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
    
    if(!stratify.by.urban) {
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
    if(!stratify.by.urban) {
      ZAggregated = A %*% Zg
      NAggregated = A %*% Ng
      pAggregated = ZAggregated / NAggregated
      pAggregated[NAggregated == 0] = NA
      
      aggregationResults = list(p=pAggregated, Z=ZAggregated, N=NAggregated)
      aggregationMatrices = list(A=A, AUrban=NULL, ARural=NULL)
    } else {
      AUrban = A$integrationMatrixUrban
      ARural = A$integrationMatrixRural
      A = A$integrationMatrix
      
      # first aggregate the numerator. The denominator will depend on the aggregation method
      ZAggregated = A %*% Zg
      ZAggregatedUrban = AUrban %*% Zg
      ZAggregatedRural = ARural %*% Zg
      
      if(use.density) {
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
      
      aggregationResults = list(region=sort(unique(areas)), p=pAggregated, Z=ZAggregated, N=NAggregated, 
                                pUrban=pAggregatedUrban, ZUrban=ZAggregatedUrban, NUrban=NAggregatedUrban, 
                                pRural=pAggregatedRural, ZRural=ZAggregatedRural, NRural=NAggregatedRural)
      aggregationMatrices = list(A=A, AUrban=AUrban, ARural=ARural)
    }
    
    list(aggregationResults=aggregationResults, aggregationMatrices=aggregationMatrices)
  }
  
  areaResults = getIntegratedPredictions(predsArea)
  
  # return results
  areaResults
}

#' @describeIn aggPop Aggregate areal populations to another areal level
#' @export
areaPopToArea = function(area.level.pop, areas.from, areas.to, 
                         stratify.by.urban=TRUE, 
                         do.fine.scale.risk=!is.null(area.level.pop$aggregationResults$pFineScaleRisk), 
                         do.smooth.risk=!is.null(area.level.pop$aggregationResults$pSmoothRisk)) {
  
  if(length(areas.from) != length(unique(areas.from))) {
    stop("areas.from must contain only unique names of areas to which we want to aggregate")
  }
  
  uniqueNames = sort(unique(areas.to))
  
  # construct row of the aggregation matrix given toArea index
  getMatrixHelper = function(i) {
    # get areas.to associated with this fromArea
    thisToArea = uniqueNames[i]
    thisFromAreas = unique(areas.from[areas.to == thisToArea])
    areaI = areas.from %in% thisFromAreas
    
    areaI
  }
  
  # construct the aggregation matrix from areas.from to areas.to
  A = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper), ncol=length(uniqueNames)))
  rownames(A) = uniqueNames
  colnames(A) = areas.from
  
  ##### aggregate populations
  # fine scale prevalence aggregation model
  getaggregationResults = function(resultNameRoot="FineScalePrevalence") {
    if(! (paste("N", resultNameRoot, sep="") %in% names(area.level.pop$aggregationResults))) {
      stop(paste0(resultNameRoot, " was not computed in input area.level.pop"))
    }
    
    nSamples = area.level.pop$aggregationResults[[paste("N", resultNameRoot, sep="")]]
    zSamples = area.level.pop$aggregationResults[[paste("Z", resultNameRoot, sep="")]]
    zSamples[is.na(zSamples)] = 0 # must set to zero temporarily so matrix multiplication works out
    
    ZAggregated =  A %*% zSamples
    NAggregated =  A %*% nSamples
    pAggregated = ZAggregated / NAggregated
    pAggregated[NAggregated == 0] = NA
    thisA=A %*% area.level.pop$aggregationMatrices[[paste("A", resultNameRoot, sep="")]]
    rownames(thisA) = uniqueNames
    
    if(stratify.by.urban) {
      nSamplesUrban = area.level.pop$aggregationResults[[paste("NUrban", resultNameRoot, sep="")]]
      zSamplesUrban = area.level.pop$aggregationResults[[paste("ZUrban", resultNameRoot, sep="")]]
      zSamplesUrban[is.na(zSamplesUrban)] = 0 # must set to zero temporarily so matrix multiplication works out
      
      nSamplesRural = area.level.pop$aggregationResults[[paste("NRural", resultNameRoot, sep="")]]
      zSamplesRural = area.level.pop$aggregationResults[[paste("ZRural", resultNameRoot, sep="")]]
      zSamplesRural[is.na(zSamplesRural)] = 0 # must set to zero temporarily so matrix multiplication works out
      
      ZAggregatedUrban =  A %*% zSamplesUrban
      NAggregatedUrban =  A %*% nSamplesUrban
      pAggregatedUrban = ZAggregatedUrban / NAggregatedUrban
      pAggregatedUrban[NAggregatedUrban == 0] = NA
      thisAUrban=A %*% area.level.pop$aggregationMatrices[[paste("AUrban", resultNameRoot, sep="")]]
      rownames(thisAUrban) = uniqueNames
      
      ZAggregatedRural =  A %*% zSamplesRural
      NAggregatedRural =  A %*% nSamplesRural
      pAggregatedRural = ZAggregatedRural / NAggregatedRural
      pAggregatedRural[NAggregatedRural == 0] = NA
      thisARural=A %*% area.level.pop$aggregationMatrices[[paste("ARural", resultNameRoot, sep="")]]
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
    
    aggregationResults = list(region=uniqueNames, p=pAggregated, Z=ZAggregated, N=NAggregated, 
                             pUrban=pAggregatedUrban, ZUrban=ZAggregatedUrban, NUrban=NAggregatedUrban, 
                             pRural=pAggregatedRural, ZRural=ZAggregatedRural, NRural=NAggregatedRural)
    aggregationMatrices = list(A=thisA, AUrban=thisAUrban, ARural=thisARural)
    
    capitalResultNameRoot = resultNameRoot
    # capitalResultNameRoot = paste(toupper(substr(capitalResultNameRoot, 1, 1)), substr(capitalResultNameRoot, 2, nchar(capitalResultNameRoot)), sep="")
    names(aggregationResults)[-1] = paste(names(aggregationResults)[-1], capitalResultNameRoot, sep="")
    names(aggregationMatrices) = paste(names(aggregationMatrices)[-1], capitalResultNameRoot, sep="")
    
    list(aggregationResults=aggregationResults, aggregationMatrices=aggregationMatrices)
  }
  
  # fine scale prevalence model 
  out = getaggregationResults("FineScalePrevalence")
  aggregationResults = out$aggregationResults
  aggregationMatrices = out$aggregationMatrices
  
  # fine scale risk aggregation model
  if(do.fine.scale.risk) {
    out = getaggregationResults("FineScaleRisk")
    resFineScaleRisk = out$aggregationResults
    resAggregationMatrices = out$aggregationMatrices
    
    # aggregationResults = merge(aggregationResults, resFineScaleRisk, by="region")
    aggregationResults = c(aggregationResults, resFineScaleRisk)
    aggregationMatrices = c(aggregationMatrices, resAggregationMatrices)
  }
  
  if(do.smooth.risk) {
    out = getaggregationResults("SmoothRisk")
    resSmoothRisk = out$aggregationResults
    resAggregationMatrices = out$aggregationMatrices
    
    # aggregationResults = merge(aggregationResults, resSmoothRisk, by="region")
    aggregationResults = c(aggregationResults, resSmoothRisk)
    aggregationMatrices = c(aggregationMatrices, resAggregationMatrices)
  }
  
  list(aggregationResults=aggregationResults, aggregationMatrices=aggregationMatrices)
}

#' 
#' `r lifecycle::badge("experimental")`
#' 
#' @describeIn simPop
#' Simulate populations and population prevalences given census frame and population density 
#' information. Uses custom spatial logit risk function and can include iid cluster 
#' level effect.
#' 
#' @export
simPopCustom = function(logit.risk.draws, sigma.epsilon.draws, easpa, pop.mat, target.pop.mat, 
                        stratify.by.urban=TRUE, validation.pixel.I=NULL, validation.cluster.I=NULL, 
                        clusters.per.pixel=NULL, 
                        do.fine.scale.risk=FALSE, do.smooth.risk=FALSE, 
                        do.smooth.risk.logistic.approx=FALSE, 
                        poppsub=NULL, subarea.level=FALSE, 
                        min1.per.subarea=TRUE, 
                        return.EA.info=FALSE, epsc=NULL) {
  
  if(is.null(poppsub) && subarea.level) {
    stop("if subarea.level is TRUE, user must specify poppsub")
  }
  
  if(!is.null(validation.pixel.I) || !is.null(validation.cluster.I) || !is.null(clusters.per.pixel)) {
    stop("validation.pixel.I, validation.cluster.I, and clusters.per.pixel not yet fully implemented")
  }
  
  ndraws = ncol(logit.risk.draws)
  
  # set default inputs
  totalEAs = sum(easpa$EATotal)
  if(!is.null(clusters.per.pixel)) {
    emptyPixels = clusters.per.pixel == 0
    if(totalEAs != sum(clusters.per.pixel))
      stop("sum(easpa$EATotal) != sum(clusters.per.pixel)")
  }
  
  # get area names
  areas = sort(unique(pop.mat$area))
  if(any(areas != easpa$area))
    stop("area names and easpa do not match pop.mat or are not in the correct order")
  
  # determine if we care about subareas (smallest areas we care about. No info of EAs per subarea)
  subarea.level = !is.null(pop.mat$subarea)
  
  ##### Line 1 (of the algorithm): take draws from the binomial process for each stratum (each row of easpa)
  # get probabilities for each pixel (or at least something proportional within each stratum)
  print("drawing EAs")
  pixelProbs = pop.mat$pop
  
  # take draws from the stratified binomial process for each posterior sample
  if(is.null(clusters.per.pixel)) {
    if(subarea.level) {
      ea.samples = rStratifiedMultnomialBySubarea(ndraws, pop.mat, easpa, stratify.by.urban, poppsub=poppsub, 
                                                 min1.per.subarea=min1.per.subarea)
    } else {
      ea.samples = rStratifiedMultnomial(ndraws, pop.mat, easpa, stratify.by.urban)
    }
  }
  
  if(!is.null(clusters.per.pixel) && !exists("ea.samples")) {
    ea.samples = matrix(rep(clusters.per.pixel, ndraws), ncol=ndraws)
  }
  
  # make matrix (or list) of pixel indices mapping matrices of EA values to matrices of pixel values
  if(!is.null(clusters.per.pixel)) {
    pixel.indices = rep(1:nrow(pop.mat), times=clusters.per.pixel) # this contains repetitions and has length == nEAs
    uniquePixelIndices = sort(unique(pixel.indices))
  } else {
    pixel.index.mat = matrix(rep(rep(1:nrow(pop.mat), ndraws), times=ea.samples), ncol=ndraws)
  }
  
  # determine which EAs are urban if necessary
  if(stratify.by.urban) {
    # urban.mat = matrix(rep(rep(pop.mat$urban, ndraws), times=c(ea.samples)), ncol=ndraws)
    if(!is.null(clusters.per.pixel)) {
      urbanVals = pop.mat$urban[pixel.indices]
      uniqueUrbanVals = pop.mat$urban[uniquePixelIndices]
    } else {
      # urban.mat = matrix(pop.mat$urban[pixel.index.mat], ncol=ndraws)
      urban.mat = NULL # don't calculate here for memory's sake
    }
  } else {
    urban.mat = NULL
  }
  
  # determine which EAs are from which area
  if(!is.null(clusters.per.pixel)) {
    areaVals = pop.mat$area[pixel.indices]
    uniqueAreaVals = pop.mat$area[uniquePixelIndices]
  } else {
    # area.mat = matrix(pop.mat$area[pixel.index.mat], ncol=ndraws)
    area.mat = NULL # don't calculate here for memory's sake
  }
  
  ##### Line 2: draw cluster effects, epsilon
  # NOTE1: we assume there are many more EAs then sampled clusters, so that 
  #       the cluster effects for each EA, including those sampled, are iid
  print("simulating EA level risks, numerators, and denominators")
  if(is.null(epsc)) {
    epsc = matrix(stats::rnorm(totalEAs*ndraws, sd=rep(sigma.epsilon.draws, each=totalEAs)), ncol=ndraws)
  }
  
  ##### Line 3: draw EA population denominators, N
  
  if(!is.null(clusters.per.pixel)) {
    if(is.null(validation.pixel.I))
      stop("clusters.per.pixel must only be set for validation, but validation.pixel.I is NULL")
    
    # in this case, every left out cluster has exactly 25 households. Simply sample target population 
    # with equal probability from each cluster/faux EA
    Ncs = sampleNMultilevelMultinomialFixed(clusters.per.pixel, ndraws=ndraws, pixel.indices=pixel.indices, 
                                            urbanVals=urbanVals, areaVals=areaVals, easpa=easpa, pop.mat=pop.mat, stratify.by.urban=stratify.by.urban, 
                                            verbose=TRUE)
  } else {
    if(return.EA.info) {
      out = sampleNMultilevelMultinomial(pixel.index.mat=pixel.index.mat, easpa.list=list(easpa), 
                                         pop.mat=pop.mat, stratify.by.urban=stratify.by.urban, verbose=TRUE, return.EA.info=return.EA.info)
      householdDraws = out$householdDraws
      Ncs = out$targetPopDraws
    } else {
      Ncs <- sampleNMultilevelMultinomial(pixel.index.mat=pixel.index.mat, easpa.list=list(easpa), 
                                       pop.mat=pop.mat, stratify.by.urban=stratify.by.urban, verbose=TRUE, return.EA.info=return.EA.info)
      householdDraws = NULL
    }
  }
  
  ##### do part of Line 7 in advance
  # calculate mu_{ic} for each EA in each pixel
  if(!is.null(clusters.per.pixel)) {
    uc = logit.risk.draws[pixel.indices,]
    muc = expit(uc + epsc)
  } else {
    uc = matrix(logit.risk.draws[cbind(rep(rep(1:nrow(logit.risk.draws), ndraws), times=c(ea.samples)), rep(1:ndraws, each=totalEAs))], ncol=ndraws)
    muc = expit(uc + epsc)
  }
  
  # calculate Z_{ic} for each EA in each pixel
  Zc = matrix(stats::rbinom(n=totalEAs * ndraws, size=Ncs, prob=as.matrix(muc)), ncol=ndraws)
  
  ##### Line 4: Aggregate appropriate values from EAs to the grid cell level
  
  # function for aggregating values for each grid cell
  getPixelColumnFromEAs = function(i, vals, applyFun=sum, popWeightMatrix=NULL) {
    # calculate levels over which to aggregate
    if(!is.null(clusters.per.pixel)) {
      indices = pixel.indices
    } else {
      indices = factor(as.character(pixel.index.mat[,i]))
    }
    
    # in this case (the LCPb model), we calculate weighted means within factor levels using popWeightMatrix
    if(!is.null(popWeightMatrix)) {
      stop("using popWeightMatrix is no longer support, since this is much slower than calculating normalized 
           weights separately and multiplying values by them outside this function")
      # applyFun = function(x) {stats::weighted.mean(x, popWeightMatrix[,i], na.rm=TRUE)}
      
      Data = data.frame(v=vals[,i], w=popWeightMatrix[,i])
      out = sapply(split(Data, indices), function(x) stats::weighted.mean(x$v,x$w))
    } else {
      if(!is.null(clusters.per.pixel)) {
        # out = tapply(vals[,i], factor(as.character(pixel.indices)), FUN=applyFun)
        out = tapply(vals[,i], indices, FUN=applyFun)
      } else {
        out = tapply(vals[,i], indices, FUN=applyFun)
      }
    }
    
    if(!is.null(clusters.per.pixel)) {
      returnValues = out
    } else {
      indices = as.numeric(names(out))
      
      returnValues = rep(NA, nrow(logit.risk.draws))
      returnValues[indices] = out
    }
    
    returnValues
  }
  
  ##### Line 5: We already did this, resulting in logit.risk.draws input
  
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
  if(do.fine.scale.risk || do.smooth.risk) {
    nPerEA = getExpectedNperEA(easpa, target.pop.mat)
  }
  
  if(do.fine.scale.risk) {
    fineScaleRisk = sapply(1:ncol(muc), getPixelColumnFromEAs, vals=muc, applyFun=function(x) {mean(x, na.rm=TRUE)})
    fineScaleRisk[!is.finite(fineScaleRisk)] = NA
    
    # in order to get valid count estimates, we also need the expected denominator per EA in each stratum:
    nSamplesFineScaleRisk = sweep(ea.samples, 1, nPerEA, "*")
    zSamplesFineScaleRisk = fineScaleRisk * nSamplesFineScaleRisk
    zSamplesFineScaleRisk[is.na(zSamplesFineScaleRisk)] = 0
  } else {
    fineScaleRisk = NULL
    zSamplesFineScaleRisk = NULL
    nSamplesFineScaleRisk = NULL
  }
  
  if(do.smooth.risk) {
    # integrate out the cluster effect in the risk
    smoothRisk = matrix(SUMMER::logitNormMean(cbind(c(as.matrix(logit.risk.draws)), rep(sigma.epsilon.draws, each=nrow(logit.risk.draws))), logisticApprox=do.smooth.risk.logistic.approx), nrow=nrow(logit.risk.draws))
    
    # get expected numerator and denominators
    nSamplesSmoothRisk = matrix(rep(target.pop.mat$pop, ndraws), ncol=ndraws)
    zSamplesSmoothRisk = sweep(smoothRisk, 1, target.pop.mat$pop, "*")
  } else {
    smoothRisk = NULL
    zSamplesSmoothRisk = NULL
    nSamplesSmoothRisk = NULL
  }
  
  ##### Extra steps: collect draws at each level and generate:
  ##### areas, preds, 
  pixel.level.pop = list(pFineScalePrevalence=pg, ZFineScalePrevalence=Zg, NFineScalePrevalence=Ng, 
                       pFineScaleRisk=fineScaleRisk, ZFineScaleRisk=zSamplesFineScaleRisk, NFineScaleRisk=nSamplesFineScaleRisk, 
                       pSmoothRisk=smoothRisk, ZSmoothRisk=zSamplesSmoothRisk, NSmoothRisk=nSamplesSmoothRisk)
  
  if(!return.EA.info) {
    pixel.level.pop
  } else {
    
    
    # return list of eaDat objects
    getEADat = function(i) {
      theseI = pixel.index.mat[,i]
      
      eaDat = data.frame(lon=pop.mat$lon[theseI], lat=pop.mat$lat[theseI], 
                         area=pop.mat$area[theseI], subarea=rep("temp", length(theseI)), 
                         urban=pop.mat$urban[theseI], east=pop.mat$east[theseI], north=pop.mat$north[theseI], 
                         popDensity=pop.mat$pop[theseI], popDensityTarget=target.pop.mat$pop[theseI], pixelIs=theseI, 
                         nHH=householdDraws[,i], N=Ncs[,i], Z=Zc[,i], 
                         pFineScalePrevalence=Zc[,i]/Ncs[,i])
      if(do.fine.scale.risk) {
        eaDat$pFineScaleRisk=muc[,i]
      }
      if(do.smooth.risk) {
        eaDat$pSmoothRisk=smoothRisk[theseI,i]
      }
      
      if(subarea.level) {
        eaDat$subarea = pop.mat$subarea[theseI]
      } else {
        eaDat$subarea = NULL
      }
      
      eaDat$pFineScalePrevalence[eaDat$N == 0] = NA
      
      eaDat
    }
    
    print("Constructing list of simulated EA data.frames")
    eaDatList = lapply(1:ndraws, getEADat)
    
    c(pixel.level.pop, list(eaDatList=eaDatList, ea.samples=ea.samples))
  }
}


#' Internal functions for population simulation
#' 
#' Functions for calculating valuable quantities and for drawing from important 
#' distributions for population simulation.
#' 
#' `r lifecycle::badge("experimental")`
#' 
#' 
#' @param easpa Census frame. See \code{\link{simPopCustom}} for details
#' @param pop.mat data.frame of pixellated grid of population densities. See \code{\link{simPopCustom}} for details
#' @param i Index
#' @param urban If TRUE, calculate only for urban part of the area. If FALSE, for only rural part
#' @param stratify.by.urban whether or not to stratify calculations by urban/rural classification
#' @param validation.pixel.I CURRENTLY FOR TESTING PURPOSES ONLY a set of indices of pixels for which we want to simulate populations (used for pixel level validation)
#' @param n Number of samples
#' @param poppsub Population per subarea. See \code{\link{simPopCustom}} for details
#' @param min1.per.subarea Whether or not to ensure there is at least 1 EA per subarea. See \code{\link{simPopCustom}} for details
#' @param method If min1.per.subarea is TRUE, the sampling method for the truncated multinomial to use with rmulitnom1. rmultinom1 automatically 
#'         switches between them depending on the number of expected samples. The methods are:
#' \describe{
#'   \item{mult1}{rejection sampling from multinomial plus 1 in each category}
#'   \item{mult}{rejection sampling from multinomial if any category has zero count}
#'   \item{indepMH}{independent Metropolis-Hastings using multinomial plus 1 distribution as proposal}
#' }
#' @param min.sample The minimum number of samples per `chunk` of samples for truncated multinomial sampling. Defaults to 1
#' @param easpsub This could either be total EAs per subarea, or subarea crossed with urban or 
#'          rural if stratify.by.urban is TRUE
#' @param size Multinomial size parameter. See \code{\link[stats]{rmultinom}}
#' @param prob Multinomial probability vector parameter. See \code{\link[stats]{rmultinom}}
#' @param max.size The maximum number of elements in a matrix drawn from the proposal distribution per sample chunk. 
#' @param max.expected.size.before.switch Max expected number of samples / k, the number of categories, before switching method
#' @param init Initial sample if method is `indepMH`
#' @param burnin Number of initial samples before samples are collected if method is `indepMH`
#' @param filter.every Store only every filter.every samples if method is i`indepMH`
#' @param zero.prob.zero.samples If TRUE, set samples for parts of prob vector that are zero to zero. Otherwise they are set to one.
#' @param allow.size.less.than.K If TRUE, then if size < the number of categories (k), returns matrix where each 
#'                     column is vector of size ones and k - size zeros. If FALSE, throws an error if size < k
#' @param clusters.per.pixel CURRENTLY FOR TESTING PURPOSES ONLY a vector of length nIntegrationPoints specifying the number of clusters per pixel if they are fixed
#' @param pixel.indices A nEA x n matrix of pixel indices associated with each EA per simulation/draw
#' @param urbanVals A nEA x n matrix of urbanicities associated with each EA per simulation/draw
#' @param areaVals A nEA x n matrix of area names associated with each EA per simulation/draw
#' @param easpa.list A list of length n with each element being of the format of easpa 
#'            giving the number of households and EAs 
#'            per stratum. It is assumed that the number of EAs per stratum is 
#'            the same in each list element. If easpa.list is a data frame, 
#'            number of households per stratum is assumed constant
#' @param ndraws Number of draws
#' @param pixel.index.mat Matrix of pixel indices associated with each EA and draw. Not 
#' required by getExpectedNperEA unless level == "EA"
#' @param urban.mat Matrix of urbanicities associated with each EA and draw
#' @param area.mat Matrix of areas associated with each EA and draw
#' @param verbose Whether to print progress as the function proceeds
#' @param return.EA.info Whether a data frame at the EA level is desired
#' @param min.HH.per.EA The minimum number of households per EA (defaults to 25, since 
#' that is the number of households sampled per DHS cluster)
#' @param fix.HH.per.EA If not NULL, the fixed number of households per EA
#' @param fix.pop.per.HH If not NULL, the fixed target population per household
#' @param level Whether to calculate results at the integration grid or EA level
#' @name simPopInternal
NULL

#' @describeIn simPopInternal Calculates expected denominator per enumeration area.
getExpectedNperEA = function(easpa, pop.mat, level=c("grid", "EA"), pixel.index.mat=NULL) {
  level = match.arg(level)
  
  # calculate the expected denominator per enumeration area in each stratum. 
  nPerEAUrban = easpa$popUrb / easpa$EAUrb
  nPerEARural = easpa$popRur / easpa$EARur
  
  # expanded the expected denominator values vector to be of length equal 
  # to the number of grid cells
  uniqueAreas = sort(unique(pop.mat$area))
  outUrban = numeric(nrow(pop.mat))
  outRural = numeric(nrow(pop.mat))
  for(i in 1:length(uniqueAreas)) {
    urbanI = getSortIndices(i, urban=TRUE, pop.mat=pop.mat, stratify.by.urban=TRUE)
    ruralI = getSortIndices(i, urban=FALSE, pop.mat=pop.mat, stratify.by.urban=TRUE)
    outUrban[urbanI] = nPerEAUrban[i]
    outRural[ruralI] = nPerEARural[i]
  }
  
  out = outUrban + outRural
  
  if(level == "EA") {
    if(is.null(pixel.index.mat)) {
      stop("if level == EA, must specify pixel.index.mat")
    }
    out = matrix(out[pixel.index.mat], ncol=ncol(pixel.index.mat))
  }
  
  out
}

#' @describeIn simPopInternal For recombining separate multinomials into the draws over all grid points
getSortIndices = function(i, urban=TRUE, pop.mat, stratify.by.urban=TRUE, validation.pixel.I=NULL) {
  
  # get area names
  areas = sort(unique(pop.mat$area))
  
  # determine which pixels and how many EAs are in this stratum
  if(stratify.by.urban) {
    includeI = (pop.mat$area == areas[i]) & (pop.mat$urban == urban)
  }
  else {
    includeI = pop.mat$area == areas[i]
  }
  
  # include only indices included within validation if necessary
  if(!is.null(validation.pixel.I)) {
    
    # convert validation.pixel.I into a logical
    temp = rep(FALSE, length(includeI))
    temp[validation.pixel.I] = TRUE
    
    # include only indices we are interested in for the validation
    includeI = includeI & temp
  }
  
  which(includeI)
}

#' @describeIn simPopInternal Gives nIntegrationPoints x n matrix of draws from the stratified multinomial with values 
#' corresponding to the value of |C^g| for each pixel, g (the number of EAs/pixel)
rStratifiedMultnomial = function(n, pop.mat, easpa, stratify.by.urban=TRUE) {
  
  # get area names
  areas = sort(unique(pop.mat$area))
  if(any(areas != easpa$area))
    stop("area names and easpa do not match pop.mat or are not in the correct order")
  
  # we will need to draw separate multinomial for each stratum. Start by 
  # creating matrix of all draws of |C^g|
  ea.samples = matrix(NA, nrow=nrow(pop.mat), ncol=n)
  
  # now draw multinomials
  if(stratify.by.urban) {
    # draw for each area crossed with urban/rural
    urbanSamples = do.call("rbind", lapply(1:length(areas), rMyMultinomial, n=n, urban=TRUE, 
                                           stratify.by.urban=stratify.by.urban, pop.mat=pop.mat, easpa=easpa))
    ruralSamples = do.call("rbind", lapply(1:length(areas), rMyMultinomial, n=n, urban=FALSE, 
                                           stratify.by.urban=stratify.by.urban, pop.mat=pop.mat, easpa=easpa))
    
    # get the indices used to recombine into the full set of draws
    urbanIndices = unlist(sapply(1:length(areas), getSortIndices, urban=TRUE, pop.mat=pop.mat, stratify.by.urban=stratify.by.urban))
    ruralIndices = unlist(sapply(1:length(areas), getSortIndices, urban=FALSE, pop.mat=pop.mat, stratify.by.urban=stratify.by.urban))
    
    # recombine into ea.samples
    ea.samples[urbanIndices,] = urbanSamples
    ea.samples[ruralIndices,] = ruralSamples
  } else {
    # draw for each area
    stratumSamples = rbind(sapply(1:length(areas), n=n, rMyMultinomial, 
                                  stratify.by.urban=stratify.by.urban, pop.mat=pop.mat, easpa=easpa))
    
    # get the indices used to recombine into the full set of draws
    stratumIndices = c(sapply(1:length(areas), getSortIndices, pop.mat=pop.mat, stratify.by.urban=stratify.by.urban))
    
    # recombine into ea.samples
    ea.samples[stratumIndices,] = stratumSamples
  }
  
  # return results
  ea.samples
}

#' 
#' `r lifecycle::badge("experimental")`
#' 
#' @describeIn simPopInternal Gives nIntegrationPoints x n matrix of draws from the stratified multinomial with values 
# corresponding to the number of EAs in each pixel
rStratifiedMultnomialBySubarea = function(n, pop.mat, easpa, stratify.by.urban=TRUE, poppsub=NULL, 
                                          min1.per.subarea=TRUE, min.sample=1) {
  if(is.null(poppsub)) {
    # use pop.mat to calculate poppsub
    stop("Calculating poppsub with pop.mat not yet implemented")
    
    # poppsub: a table with the following variables:
    # subarea
    # area
    # popUrb
    # popRur
    # popTotal
  }
  
  # get area names
  areas = sort(unique(pop.mat$area))
  subareas = sort(unique(pop.mat$subarea))
  if(any(areas != easpa$area))
    stop("area names and easpa do not match pop.mat or are not in the correct order")
  
  # we will need to draw separate multinomial for each stratum
  
  # create temporary pop.mat, except with one row for each constituency
  popSubarea.mat = pop.mat[1:length(subareas),]
  popSubarea.mat$area = poppsub$area
  popSubarea.mat$subarea = poppsub$subarea
  if(stratify.by.urban) {
    popSubarea.mat$urban = FALSE
    popSubarea.mat = rbind(popSubarea.mat, popSubarea.mat)
    popSubarea.mat$urban[1:length(subareas)] = TRUE
    popSubarea.mat$pop[1:length(subareas)] = poppsub$popUrb
    popSubarea.mat$pop[(length(subareas) + 1):(2 * length(subareas))] = poppsub$popRur
  } else {
    popSubarea.mat$pop = poppsub$popTotal
  }
  
  # now draw multinomials
  if(stratify.by.urban) {
    # draw for each constituency in each area crossed with urban/rural
    urbanSamplesCon = do.call("rbind", lapply(1:length(areas), rMyMultinomial, n=n, urban=TRUE, 
                                              stratify.by.urban=stratify.by.urban, pop.mat=popSubarea.mat, easpa=easpa, 
                                              min1.per.subarea=min1.per.subarea, method="mult", 
                                              min.sample=min.sample))
    ruralSamplesCon = do.call("rbind", lapply(1:length(areas), rMyMultinomial, n=n, urban=FALSE, 
                                              stratify.by.urban=stratify.by.urban, pop.mat=popSubarea.mat, easpa=easpa, 
                                              min1.per.subarea=min1.per.subarea, method="mult", 
                                              min.sample=min.sample))
    
    # get the indices used to recombine into the full set of draws for the subareas
    urbanIndicesCon = unlist(sapply(1:length(areas), getSortIndices, urban=TRUE, pop.mat=popSubarea.mat, stratify.by.urban=stratify.by.urban))
    ruralIndicesCon = unlist(sapply(1:length(areas), getSortIndices, urban=FALSE, pop.mat=popSubarea.mat, stratify.by.urban=stratify.by.urban)) - length(urbanIndicesCon)
    
    # recombine into ea.samples for the subareas
    urbanSamplesCon[urbanIndicesCon,] = urbanSamplesCon
    ruralSamplesCon[ruralIndicesCon,] = ruralSamplesCon
    
    # draw for each pixel crossed with urban/rural
    urbanSamples = do.call("rbind", lapply(1:length(subareas), rMyMultinomialSubarea, n=n, urban=TRUE, 
                                           stratify.by.urban=stratify.by.urban, pop.mat=pop.mat, easpsub=urbanSamplesCon))
    ruralSamples = do.call("rbind", lapply(1:length(subareas), rMyMultinomialSubarea, n=n, urban=FALSE, 
                                           stratify.by.urban=stratify.by.urban, pop.mat=pop.mat, easpsub=ruralSamplesCon))
    
    # get the indices used to recombine into the full set of draws
    tempPopMat = pop.mat
    tempPopMat$area = tempPopMat$subarea
    urbanIndices = unlist(sapply(1:length(subareas), getSortIndices, urban=TRUE, pop.mat=tempPopMat, stratify.by.urban=stratify.by.urban))
    ruralIndices = unlist(sapply(1:length(subareas), getSortIndices, urban=FALSE, pop.mat=tempPopMat, stratify.by.urban=stratify.by.urban))
    
    # create matrix of all draws of |C^g| and recombine urban/rural draws into ea.samples
    ea.samples = matrix(NA, nrow=nrow(pop.mat), ncol=n)
    ea.samples[urbanIndices,] = urbanSamples
    ea.samples[ruralIndices,] = ruralSamples
  } else {
    # draw for each constituency in each area crossed with urban/rural
    samplesCon = do.call("rbind", lapply(1:length(areas), rMyMultinomial, n=n, urban=TRUE, 
                                         stratify.by.urban=stratify.by.urban, pop.mat=popSubarea.mat, easpa=easpa, 
                                         min1.per.subarea=min1.per.subarea, method="mult", 
                                         min.sample=min.sample))
    
    # get the indices used to recombine into the full set of draws for the subareas
    indicesCon = unlist(sapply(1:length(areas), getSortIndices, pop.mat=popSubarea.mat, stratify.by.urban=stratify.by.urban))
    
    # recombine into ea.samples for the subareas
    samplesCon[indicesCon,] = samplesCon
    
    # draw for each pixel in each constituency
    stratumSamples = rbind(sapply(1:length(subareas), n=n, rMyMultinomialSubarea, 
                                  stratify.by.urban=stratify.by.urban, pop.mat=pop.mat, easpsub=samplesCon))
    
    # get the indices used to recombine into the full set of draws
    tempPopMat = pop.mat
    tempPopMat$area = tempPopMat$subarea
    stratumIndices = c(sapply(1:length(subareas), getSortIndices, pop.mat=tempPopMat, stratify.by.urban=stratify.by.urban))
    
    # create matrix of all draws of |C^g|
    ea.samples = matrix(NA, nrow=nrow(pop.mat), ncol=n)
    ea.samples[stratumIndices,] = stratumSamples
  }
  
  # return results
  ea.samples
}

#' 
#' `r lifecycle::badge("experimental")`
#' 
#' @describeIn simPopInternal 
rMyMultinomial = function(n, i, stratify.by.urban=TRUE, urban=TRUE, pop.mat=NULL, easpa=NULL, min1.per.subarea=FALSE, 
                          method=c("mult1", "mult", "indepMH"), min.sample=1) {
  method = match.arg(method)
  
  # get area names
  areas = sort(unique(pop.mat$area))
  if(any(areas != easpa$area))
    stop("area names and easpa do not match pop.mat or are not in the correct order")
  
  # determine which pixels and how many EAs are in this stratum
  if(stratify.by.urban) {
    includeI = (pop.mat$area == areas[i]) & (pop.mat$urban == urban)
    nEA = ifelse(urban, easpa$EAUrb[i], easpa$EARur[i])
  }
  else {
    includeI = pop.mat$area == areas[i]
    nEA = ifelse(urban, easpa$EAUrb[i], easpa$EATotal[i])
  }
  
  # sample from the pixels if this stratum exists
  if(sum(includeI) == 0)
    return(matrix(nrow=0, ncol=n))
  thesePixelProbs = pop.mat$pop[includeI]
  if(any(thesePixelProbs > 0)) {
    if(!min1.per.subarea) {
      stats::rmultinom(n, nEA, prob=thesePixelProbs)
    } else {
      rmultinom1(n, nEA, prob=thesePixelProbs, method=method, allow.size.less.than.K=TRUE, min.sample=min.sample)
    }
  } else {
    matrix(0, nrow=length(thesePixelProbs), ncol=n)
  }
}

#' 
#' `r lifecycle::badge("experimental")`
#' 
#' @describeIn simPopInternal 
rMyMultinomialSubarea = function(n, i, easpsub, stratify.by.urban=TRUE, urban=TRUE, pop.mat=NULL) {
  
  # get constituency names
  subareas = sort(unique(pop.mat$subarea))
  
  # determine which pixels and how many EAs are in this stratum
  if(stratify.by.urban) {
    includeI = (pop.mat$subarea == subareas[i]) & (pop.mat$urban == urban)
  }
  else {
    includeI = pop.mat$subarea == subareas[i]
  }
  nEA = easpsub[i,]
  
  # sample from the pixels if this stratum exists
  if(sum(includeI) == 0){
    if(any(nEA != 0))
      stop(paste0("no valid pixels to put EAs in for subarea ", as.character(subareas[i]), " and urban level ", urban))
    return(matrix(nrow=0, ncol=n))
  }
  
  thesePixelProbs = pop.mat$pop[includeI]
  nonZeroEAs = nEA != 0
  out = matrix(0, nrow=sum(includeI), ncol=n)
  if(any(nonZeroEAs)) {
    out[,nonZeroEAs] = sapply(nEA[nonZeroEAs], stats::rmultinom, n=1, prob=thesePixelProbs)
  }
  out
}

#' 
#' `r lifecycle::badge("experimental")`
#' 
#' @describeIn simPopInternal Random (truncated) multinomial draws conditional on the number of each type being at least one
rmultinom1 = function(n=1, size, prob, max.size=8000*8000, method=c("mult1", "mult", "indepMH"), verbose=FALSE, min.sample=100, 
                      max.expected.size.before.switch=1000*1e7, init=NULL, burnin=floor(n/4), filter.every=10, zero.prob.zero.samples=TRUE, 
                      allow.size.less.than.K=FALSE) {
  method = match.arg(method)
  prob = prob*(1/sum(prob))
  
  if(zero.prob.zero.samples && any(prob == 0)) {
    zero = prob == 0
    out = matrix(0, nrow=length(prob), ncol=n)
    
    if(sum(!zero) > 0) {
      out[!zero,] = rmultinom1(n, size, prob[!zero], max.size, method, verbose, min.sample, max.expected.size.before.switch, init, burnin, 
                               filter.every, zero.prob.zero.samples, allow.size.less.than.K)
    }
    
    return(out)
  }
  
  k = length(prob)
  if(allow.size.less.than.K && (size <= k)) {
    return(replicate(n, as.numeric(1:k %in% sample(1:k, size, replace=FALSE))))
  } else if(size < k) {
    stop("size < k but rmultinom1 requires at least 1 sample per multinomial type")
  }
  
  maxSamples = floor(max.size / k)
  averageProbMult = prod((size/k)*prob)
  
  if(method != "indepMH")
    samples = matrix(NA, nrow=k, ncol=n)
  else
    samples = matrix(NA, nrow=k, ncol=round(n*filter.every))
  if(method == "mult1") {
    averagex = 1 + (size-k)*prob
    averageProb = (size-k) / (prod(averagex))
    
    while(any(is.na(samples))) {
      # calculate the number of remaining samples
      samplesLeft = sum(apply(samples, 2, function(x) {any(is.na(x))}))
      
      # approximate expected number of samples so that, after some are rejected, we will 
      # have the right number of samples. Kappa is approximated for when p_1=p_2=...=p_k
      expectedSamples = ceiling(samplesLeft/averageProb)
      # logMStar = lfactorial(size) - lfactorial(size-k) + sum(log(prob)) - log(size-k+1)
      # avgAcceptProb = exp(-logMStar)
      # logKappaApprox = lchoose(size + k - 1, k - 1) - lchoose(size - 1, k - 1)
      # logM = logMStar -logKappaApprox
      # avgAcceptProb = exp(-logM)
      # expectedSamples = ceiling(samplesLeft/avgAcceptProb)
      
      if(expectedSamples*k > max.expected.size.before.switch) {
        warning("too many samples expected with method=='mult1'. Switching to method=='indepMH'")
        return(rmultinom1(n, size, prob, max.size, method="indepMH", verbose, min.sample, 
                          max.expected.size.before.switch, init, burnin, filter.every, 
                          zero.prob.zero.samples, allow.size.less.than.K))
      }
      
      # sample expectedSamples times a fudge factor, but make sure we don't get past memory limit
      thisNumberOfSamples = max(min.sample, min(maxSamples, expectedSamples * 1.1))
      if(verbose)
        print(paste0("Sampling ", thisNumberOfSamples, ". Sampled ", n-samplesLeft, "/", n, ". Expected remaining samples: ", expectedSamples))
      thisSamples = matrix(1 + stats::rmultinom(thisNumberOfSamples, size-k, prob=prob), nrow=length(prob))
      
      # calculate accept probabilities
      thisProbs = exp(log(size-k) - apply(log(thisSamples), 2, sum))
      # thisProbs = (size-k) / apply(thisSamples, 2, prod)
      if(verbose) {
        print(paste0("Max sampled accept prob: ", max(thisProbs), ". Mean sampled accept prob: ", mean(thisProbs)))
        print(paste0("Expected number of samples based on avg sampled acceptance prob: ", ceiling(samplesLeft/mean(thisProbs))))
      }
      
      # reject relevant samples
      u = stats::runif(thisNumberOfSamples)
      thisSamples = matrix(thisSamples[,u<thisProbs], nrow=length(prob))
      
      # remove excess samples if necessary
      totalSamples = ncol(thisSamples) + n - samplesLeft
      if(totalSamples > n) {
        thisSamples = matrix(thisSamples[,1:samplesLeft], nrow=length(prob))
      }
      
      # add in accepted samples, if any
      if(ncol(thisSamples) > 0) {
        samples[,(n-samplesLeft+1):(n-samplesLeft+ncol(thisSamples))] = thisSamples
      } else {
        warning(paste0("no samples accepted this round out of ", thisNumberOfSamples, " total. Redrawing..."))
      }
    }
  } else if(method == "mult") {
    
    while(any(is.na(samples))) {
      # calculate the number of remaining samples
      samplesLeft = sum(apply(samples, 2, function(x) {any(is.na(x))}))
      
      # approximate expected number of samples so that, after some are rejected, we will 
      # have the right number of samples
      expectedSamples = ceiling(samplesLeft/averageProbMult)
      
      if(expectedSamples*k > max.expected.size.before.switch) {
        warning("too many samples expected with method=='mult'. Switching to method=='mult1'")
        return(rmultinom1(n, size, prob, max.size, method="mult1", verbose, min.sample, max.expected.size.before.switch, 
                          init, burnin, filter.every, 
                          zero.prob.zero.samples, allow.size.less.than.K))
      }
      
      # sample expectedSamples times a fudge factor, but make sure we don't get past memory limit
      thisNumberOfSamples = max(min.sample, min(maxSamples, expectedSamples * 1.1))
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
        warning(paste0("no samples accepted this round out of ", thisNumberOfSamples, " total. Redrawing..."))
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
    for (i in 1:burnin) {
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
    
    # estimate autocorrelation in samples (if it is na then all samples have the same value there)
    calcCor = apply(samples, 1, function(x) {stats::cor(x[-1], x[-length(x)])})
    calcCor[is.na(calcCor)] = 1
    print(paste0("estimated mean (max) lag 1 correlation after filtering: ", mean(calcCor^filter.every), " (", max(calcCor^filter.every), ")"))
    
    # filter out samples to reduce autocorrelation
    samples = samples[,seq(from=1, to=ncol(samples), by=filter.every)]
  }
  
  samples
}

#' 
#' `r lifecycle::badge("experimental")`
#' 
#' @describeIn simPopInternal Take multilevel multinomial draws first from joint distribution of 
#' number of households per EA given the total per stratum, and then from the joint 
#' distribution of the total target population per household given 
#' the total per stratum
sampleNMultilevelMultinomial = function(ndraws = ncol(pixel.index.mat), pixel.index.mat=NULL, urban.mat=NULL, area.mat=NULL, easpa.list, 
                                        pop.mat, stratify.by.urban=TRUE, verbose=TRUE, return.EA.info=FALSE, 
                                        min.HH.per.EA=25, fix.HH.per.EA=NULL, fix.pop.per.HH=NULL) {
  
  if(length(easpa.list) == 1) {
    easpa.list = replicate(ndraws, easpa.list[[1]], simplify=FALSE)
  }
  areas = easpa.list[[1]]$area
  
  if((is.null(area.mat) || is.null(urban.mat)) && is.null(pixel.index.mat)) {
    stop("user must either supply pixel.index.mat or both area.mat and urban.mat")
  }
  if(is.null(area.mat)) {
    areaIs = match(pop.mat$area, sort(unique(areas))) # convert from character to indices to save memory
    area.mat = matrix(areaIs[pixel.index.mat], ncol=ndraws)
  }
  if(is.null(urban.mat)) {
    urban.mat = matrix(pop.mat$urban[pixel.index.mat], ncol=ndraws)
  }
  
  # start by drawing the totals, then divide households amongst EAs, then divide target population amongst households. 
  # Make sure there are at least 25 households per EA (only divide the rest randomly)
  
  ##### Draw the totals
  
  # get the total number of enumeration areas per stratum (this does not change between draws)
  areasIs = match(areas, sort(unique(areas))) # convert from character to indices
  totalEAsUrban = easpa.list[[1]]$EAUrb
  totalEAsRural = easpa.list[[1]]$EARur
  totalEAs = easpa.list[[1]]$EATotal
  nEAs = sum(totalEAs)
  
  ##### draw the total target population per enumeration area
  
  targetPopDraws = matrix(nrow=nEAs, ncol=ndraws)
  if(return.EA.info) {
    householdDraws = matrix(nrow=nEAs, ncol=ndraws)
  }
  
  # Draw the number of households per stratum area that will be randomly distributed (total minus the minimum min.HH.per.EA)
  if(stratify.by.urban) {
    totalHouseholdsUrban = sweep(matrix(sapply(easpa.list, function(x) {x$HHUrb}), ncol=length(easpa.list)), 1, -min.HH.per.EA*totalEAsUrban, "+")
    totalHouseholdsRural = sweep(matrix(sapply(easpa.list, function(x) {x$HHRur}), ncol=length(easpa.list)), 1, -min.HH.per.EA*totalEAsRural, "+")
    totalChildrenUrban = matrix(sapply(easpa.list, function(x) {x$popUrb}), ncol=length(easpa.list))
    totalChildrenRural = matrix(sapply(easpa.list, function(x) {x$popRur}), ncol=length(easpa.list))
  } else {
    totalHouseholds = sweep(sapply(matrix(easpa.list, function(x) {x$HHTotal}), ncol=length(easpa.list)), 1, -min.HH.per.EA*totalEAs, "+")
    totalChildren = matrix(sapply(easpa.list, function(x) {x$popHHTotal}), ncol=length(easpa.list))
  }
  
  # distribute the households throughout the enumeration areas with multinomial distribution, then 
  # distribute the target population amongst the households, also with a multinomial distribution
  for(i in 1:length(areasIs)) {
    thisArea = areasIs[i]
    thisAreaL = area.mat==thisArea
    
    # print progress if in verbose mode
    if(verbose) {
      print(paste0("drawing Ns for each EA for area ", thisArea, " (", i, "/", length(areas), ")"))
    }
    
    # draw households per EA (make sure there are any rural EAs)
    if(stratify.by.urban) {
      if(totalEAsUrban[i] != 0) {
        householdDrawsUrban <- matrix(sapply(totalHouseholdsUrban[i,], stats::rmultinom, n=1, prob=rep(1/totalEAsUrban[i], totalEAsUrban[i])), nrow=totalEAsUrban[i], ncol=ndraws) + min.HH.per.EA
      }
      if(totalEAsRural[i] != 0) {
        householdDrawsRural <- matrix(sapply(totalHouseholdsRural[i,], stats::rmultinom, n=1, prob=rep(1/totalEAsRural[i], totalEAsRural[i])), nrow=totalEAsRural[i], ncol=ndraws) + min.HH.per.EA
      }
      
      # if we must return EA info, we must return the household draws for each EA:
      if(return.EA.info) {
        if(totalEAsUrban[i] != 0) {
          householdDraws[thisAreaL & urban.mat] = householdDrawsUrban
        }
        if(totalEAsRural[i] != 0) {
          householdDraws[thisAreaL & !urban.mat] = householdDrawsRural
        }
      }
    } else {
      householdDraws = matrix(sapply(totalHouseholds[i,], stats::rmultinom, n=1, prob=rep(1/totalEAs[i], totalEAs[i])), nrow=totalEAs[i], ncol=ndraws) + min.HH.per.EA
    }
    
    # drawing target population per EA
    if(is.null(fix.pop.per.HH)) {
      # draw target pop per EA with probability proportional to the number of households
      if(stratify.by.urban) {
        if(totalEAsUrban[i] != 0) {
          probsUrban = sweep(householdDrawsUrban, 2, 1 / colSums(householdDrawsUrban), "*")
          targetPopDraws[thisAreaL & urban.mat] = sapply(1:ndraws, function(j) {stats::rmultinom(1, totalChildrenUrban[i,j], probsUrban[,j])})
        }
        
        if(totalEAsRural[i] != 0) {
          probsRural = sweep(householdDrawsRural, 2, 1 / colSums(householdDrawsRural), "*")
          targetPopDraws[thisAreaL & !urban.mat] = sapply(1:ndraws, function(j) {stats::rmultinom(1, totalChildrenRural[i,j], probsRural[,j])})
        }
      } else {
        probs = sweep(householdDraws, 2, 1 / colSums(householdDraws), "*")
        targetPopDraws[thisAreaL] = sapply(1:ndraws, function(j) {stats::rmultinom(1, totalChildren[i,j], probs[,j])})
      }
    } else {
      # set target pop per EA based on fixed number per household
      if(stratify.by.urban) {
        if(totalEAsUrban[i] != 0) {
          targetPopDraws[thisAreaL & urban.mat] = fix.pop.per.HH*householdDrawsUrban
        }
        if(totalEAsRural[i] != 0) {
          targetPopDraws[thisAreaL & !urban.mat] = fix.pop.per.HH*householdDrawsRural
        }
      } else {
        targetPopDraws[thisAreaL] = fix.pop.per.HH*householdDraws
      }
    }
  }
  
  ##### Return results
  if(!return.EA.info) {
    targetPopDraws
  } else {
    list(householdDraws=householdDraws, targetPopDraws=targetPopDraws)
  }
}

#' 
#' `r lifecycle::badge("experimental")`
#' 
#' @describeIn simPopInternal Same as sampleNMultilevelMultinomial, except the number of EAs per pixel is fixed
sampleNMultilevelMultinomialFixed = function(clusters.per.pixel, ndraws=ncol(pixel.indices), pixel.indices=NULL, 
                                             urbanVals=NULL, areaVals=NULL, easpa, pop.mat, stratify.by.urban=TRUE, 
                                             verbose=TRUE) {
  
  # set default inputs
  if((is.null(areaVals) || is.null(urbanVals)) && is.null(pixel.indices)) {
    stop("user must either supply pixel.indices or both areaVals and urbanVals")
  }
  if(is.null(areaVals)) {
    areaVals = matrix(pop.mat$area[pixel.indices], ncol=ndraws)
  }
  if(is.null(urbanVals)) {
    urbanVals = matrix(pop.mat$urban[pixel.indices], ncol=ndraws)
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
  
  if(nEAs != sum(clusters.per.pixel)) {
    stop("sum(easpa$EATotal) != sum(clusters.per.pixel)")
  }
  
  ##### draw the total target population per EA
  targetPopDraws = matrix(nrow=nEAs, ncol=ndraws)
  
  # Draw the number of households per stratum area that will be randomly distributed (total minus the minimum 25)
  if(stratify.by.urban) {
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
    if(stratify.by.urban) {
      if(totalEAsUrban[i] != 0) {
        if(any(totalHouseholdsUrban != 0)) {
          householdDrawsUrban = stats::rmultinom(n=ndraws, size=totalHouseholdsUrban[i], prob=rep(1/totalEAsUrban[i], totalEAsUrban[i])) + 25
        } else {
          householdDrawsUrban = matrix(rep(25, totalEAsUrban[i]*ndraws), ncol=ndraws)
        }
      }
      
      if(totalEAsRural[i] != 0) {
        if(any(totalHouseholdsRural != 0)) {
          householdDrawsRural = stats::rmultinom(n=ndraws, size=totalHouseholdsRural[i], prob=rep(1/totalEAsRural[i], totalEAsRural[i])) + 25
        } else {
          householdDrawsRural = matrix(rep(25, totalEAsRural[i]*ndraws), ncol=ndraws)
        }
      }
    } else {
      if(any(totalHouseholdsRural != 0)) {
        householdDraws = stats::rmultinom(n=ndraws, size=totalHouseholds[i], prob=rep(1/totalEAs[i], totalEAs[i])) + 25
      } else {
        householdDraws = matrix(rep(25, totalEAs[i]*ndraws), ncol=ndraws)
      }
    }
    
    # drawing target population per EA, with probability proportional to the number of households
    if(stratify.by.urban) {
      
      if(totalEAsUrban[i] != 0) {
        probsUrban = sweep(householdDrawsUrban, 2, 1 / colSums(householdDrawsUrban), "*")
        targetPopDraws[areaVals==thisArea & urbanVals] = sapply(1:ndraws, function(j) {stats::rmultinom(1, totalChildrenUrban[i], probsUrban[,j])})
      }
      
      if(totalEAsRural[i] != 0) {
        probsRural = sweep(householdDrawsRural, 2, 1 / colSums(householdDrawsRural), "*")
        targetPopDraws[areaVals==thisArea & !urbanVals] = sapply(1:ndraws, function(j) {stats::rmultinom(1, totalChildrenRural[i], probsRural[,j])})
      }
      
    } else {
      probs = sweep(householdDraws, 2, 1 / colSums(householdDraws), "*")
      targetPopDraws[areaVals==thisArea] = sapply(1:ndraws, function(j) {stats::rmultinom(1, totalChildren[i], probs[,j])})
    }
  }
  
  ##### Return results
  targetPopDraws
}