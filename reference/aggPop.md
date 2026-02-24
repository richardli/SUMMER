# Aggregate populations to the specified areal level

Takes simulated populations and aggregates them to the specified areal
level. Also calculates the aggregated risk and prevalence.

## Usage

``` r
pixelPopToArea(
  pixel.level.pop,
  ea.samples,
  areas,
  stratify.by.urban = TRUE,
  target.pop.mat = NULL,
  do.fine.scale.risk = !is.null(pixel.level.pop$fineScaleRisk$p),
  do.smooth.risk = !is.null(pixel.level.pop$smoothRisk$p)
)

areaPopToArea(
  area.level.pop,
  areas.from,
  areas.to,
  stratify.by.urban = TRUE,
  do.fine.scale.risk = !is.null(area.level.pop$aggregationResults$pFineScaleRisk),
  do.smooth.risk = !is.null(area.level.pop$aggregationResults$pSmoothRisk)
)
```

## Arguments

- pixel.level.pop:

  pixel level population information that we want aggregate. In the same
  format as output from
  [`simPopCustom`](https://richardli.github.io/SUMMER/reference/simPop.md)

- ea.samples:

  nIntegrationPoint x nsim matrix of the number of enumeration areas per
  pixel sampled in the input pixel level population

- areas:

  character vector of length nIntegrationPoints of area names over which
  we want to aggregate. Can also be subareas

- stratify.by.urban:

  whether or not to stratify simulations by urban/rural classification

- target.pop.mat:

  pixellated grid data frame with variables `lon`, `lat`, `pop` (target
  population), `area`, `subareas` (if subarea.level is TRUE), `urban`
  (if stratify.by.urban is TRUE), `east`, and `north`

- do.fine.scale.risk:

  whether or not to calculate the fine scale risk in addition to the
  prevalence. See details

- do.smooth.risk:

  Whether or not to calculate the smooth risk in addition to the
  prevalence. See details

- area.level.pop:

  output of
  [`simPopCustom`](https://richardli.github.io/SUMMER/reference/simPop.md)
  containing pixel level information about the population of interest

- areas.from:

  character vector of length equal to the number of areas from which we
  would like to aggregate containing the unique names of the areas. Can
  also be subareas, but these are smaller than the "to areas", and each
  "from area" must be entirely contained in a single "to area"

- areas.to:

  character vector of length equal to the number of areas from which we
  would like to aggregate containing the names of the areas containing
  with each respective \`from' area. Can also be a set of subareas, but
  these are larger than the "from areas".

## Value

A list containing elements `fineScalePrevalence` and `fineScaleRisk`.
Each of these are in turn lists with aggregated prevalence and risk for
the area of interest, containg the following elements, were paranethesis
indicate the elements for the fineScaleRisk model rather than
fineScalePrevalence:

- p:

  Aggregated prevalence (risk), calculated as aggregate of Z divided by
  aggregate of N

- Z:

  Aggregated (expected) population numerator

- N:

  Aggregated (expected) population denominator

- pUrban:

  Aggregated prevalence (risk) in urban part of the area, calculated as
  aggregate of Z divided by aggregate of N

- ZUrban:

  Aggregated (expected) population numerator in urban part of the area

- NUrban:

  Aggregated (expected) population denominator in urban part of the area

- pRural:

  Aggregated prevalence (risk) in rural part of the area, calculated as
  aggregate of Z divided by aggregate of N

- ZRural:

  Aggregated (expected) population numerator in rural part of the area

- NRural:

  Aggregated (expected) population denominator in rural part of the area

- A:

  Aggregation matrix used to aggregate from pixel level to areal level

- AUrban:

  Aggregation matrix used to aggregate from pixel level to urban part of
  the areal level

- ARural:

  Aggregation matrix used to aggregate from pixel level to rural part of
  the areal level

## Details

**\[experimental\]**

## Functions

- `pixelPopToArea()`: Aggregate from pixel to areal level

- `areaPopToArea()`: Aggregate areal populations to another areal level

## References

Paige, John, Geir-Arne Fuglstad, Andrea Riebler, and Jon Wakefield.
"Spatial aggregation with respect to a population distribution: Impact
on inference." Spatial Statistics 52 (2022): 100714.

## See also

`areaPopToArea`

## Author

John Paige

## Examples

``` r
if (FALSE) { # \dontrun{
# download Kenya GADM shapefiles from SUMMERdata github repository
githubURL <- paste0("https://github.com/paigejo/SUMMERdata/blob/main/data/", 
                    "kenyaMaps.rda?raw=true")
tempDirectory = "~/"
mapsFilename = paste0(tempDirectory, "/kenyaMaps.rda")
if(!file.exists(mapsFilename)) {
  download.file(githubURL,mapsFilename)
}

# load it in
out = load(mapsFilename)
out
kenyaMesh <- fmesher::fm_as_fm(kenyaMesh)

adm1@data$NAME_1 = as.character(adm1@data$NAME_1)
adm1@data$NAME_1[adm1@data$NAME_1 == "Trans Nzoia"] = "Trans-Nzoia"
adm1@data$NAME_1[adm1@data$NAME_1 == "Elgeyo-Marakwet"] = "Elgeyo Marakwet"
adm2@data$NAME_1 = as.character(adm2@data$NAME_1)
adm2@data$NAME_1[adm2@data$NAME_1 == "Trans Nzoia"] = "Trans-Nzoia"
adm2@data$NAME_1[adm2@data$NAME_1 == "Elgeyo-Marakwet"] = "Elgeyo Marakwet"

# some Admin-2 areas have the same name
adm2@data$NAME_2 = as.character(adm2@data$NAME_2)
adm2@data$NAME_2[(adm2@data$NAME_1 == "Bungoma") & 
                   (adm2@data$NAME_2 == "Lugari")] = "Lugari, Bungoma"
adm2@data$NAME_2[(adm2@data$NAME_1 == "Kakamega") & 
                   (adm2@data$NAME_2 == "Lugari")] = "Lugari, Kakamega"
adm2@data$NAME_2[(adm2@data$NAME_1 == "Meru") & 
                   (adm2@data$NAME_2 == "Igembe South")] = "Igembe South, Meru"
adm2@data$NAME_2[(adm2@data$NAME_1 == "Tharaka-Nithi") & 
                   (adm2@data$NAME_2 == "Igembe South")] = "Igembe South, Tharaka-Nithi"

# The spatial area of unknown 8 is so small, it causes problems unless its removed or 
# unioned with another subarea. Union it with neighboring Kakeguria:
newadm2 = adm2
unknown8I = which(newadm2$NAME_2 == "unknown 8")
newadm2$NAME_2[newadm2$NAME_2 %in% c("unknown 8", "Kapenguria")] <- 
  "Kapenguria + unknown 8"
admin2.IDs <- newadm2$NAME_2

newadm2@data = cbind(newadm2@data, NAME_2OLD = newadm2@data$NAME_2)
newadm2@data$NAME_2OLD = newadm2@data$NAME_2
newadm2@data$NAME_2 = admin2.IDs
newadm2$NAME_2 = admin2.IDs
temp <- terra::aggregate(as(newadm2, "SpatVector"), by="NAME_2")

library(sf)
temp <- sf::st_as_sf(temp)
temp <- sf::as_Spatial(temp)

tempData = newadm2@data[-unknown8I,]
tempData = tempData[order(tempData$NAME_2),]
newadm2 <- sp::SpatialPolygonsDataFrame(temp, tempData, match.ID = F)
adm2 = newadm2

# download 2014 Kenya population density TIF file

githubURL <- paste0("https://github.com/paigejo/SUMMERdata/blob/main/data/", 
                    "Kenya2014Pop/worldpop_total_1y_2014_00_00.tif?raw=true")
popTIFFilename = paste0(tempDirectory, "/worldpop_total_1y_2014_00_00.tif")
if(!file.exists(popTIFFilename)) {
  download.file(githubURL,popTIFFilename)
}

# load it in
pop = terra::rast(popTIFFilename)

east.lim = c(-110.6405, 832.4544)
north.lim = c(-555.1739, 608.7130)

## Construct poppsubKenya, a table of urban/rural general population totals 
## in each subarea. Technically, this is not necessary since we can load in 
## poppsubKenya via data(kenyaPopulationData). First, we will need to calculate 
## the areas in km^2 of the areas and subareas

# use Lambert equal area projection of areas (Admin-1) and subareas (Admin-2)
midLon = mean(adm1@bbox[1,])
midLat = mean(adm1@bbox[2,])
p4s = paste0("+proj=laea +x_0=0 +y_0=0 +lon_0=", midLon, 
             " +lat_0=", midLat, " +units=km")

adm1_sf = st_as_sf(adm1)
adm1proj_sf = st_transform(adm1_sf, p4s)
adm1proj = as(adm1proj_sf, "Spatial")

adm2_sf = st_as_sf(adm2)
adm2proj_sf = st_transform(adm2_sf, p4s)
adm2proj = as(adm2proj_sf, "Spatial")

# now calculate spatial area in km^2
admin1Areas = as.numeric(st_area(adm1proj_sf))
admin2Areas = as.numeric(st_area(adm2proj_sf))

areapaKenya = data.frame(area=adm1proj@data$NAME_1, spatialArea=admin1Areas)
areapsubKenya = data.frame(area=adm2proj@data$NAME_1, subarea=adm2proj@data$NAME_2, 
                           spatialArea=admin2Areas)

# Calculate general population totals at the subarea (Admin-2) x urban/rural 
# level and using 1km resolution population grid. Assign urbanicity by 
# thresholding population density based on estimated proportion population 
# urban/rural, making sure total area (Admin-1) urban/rural populations in 
# each area matches poppaKenya.

data(kenyaPopulationData)
pop.matKenya <- makePopIntegrationTab(
  km.res=5, pop=pop, domain.map.dat=adm0,
  east.lim=east.lim, north.lim=north.lim, map.projection=projKenya,
  poppa = poppaKenya, poppsub=poppsubKenya, 
  area.map.dat = adm1, subarea.map.dat = adm2,
  areaNameVar = "NAME_1", subareaNameVar="NAME_2")

##### Now we make a model for the risk. We will use an SPDE model with these 
##### parameters for the linear predictor on the logist scale, which are chosen 
##### to be of practical interest:
beta0=-2.9 # intercept
gamma=-1 # urban effect
rho=(1/3)^2 # spatial variance
eff.range = 400 # effective spatial range in km
sigma.epsilon=sqrt(1/2.5) # cluster (nugget) effect standard deviation

# simulate the population! Note that this produces multiple dense 
# nEA x nsim and nIntegrationPoint x nsim matrices. In the future 
# sparse matrices will and chunk by chunk computations may be incorporated.
simPop = simPopSPDE(nsim=1, easpa=easpaKenyaNeonatal, 
                    pop.mat=pop.matKenya, target.pop.mat=pop.matKenya, 
                    poppsub=poppsubKenya, spde.mesh=kenyaMesh, 
                    marg.var=rho, sigma.epsilon=sigma.epsilon, 
                    gamma=gamma, eff.range=eff.range, beta0=beta0, 
                    seed=123, inla.seed=12, n.HH.sampled=25, 
                    stratify.by.urban=TRUE, subarea.level=TRUE, 
                    do.fine.scale.risk=TRUE, 
                    min1.per.subarea=TRUE)

pixelPop = simPop$pixelPop
subareaPop = pixelPopToArea(pixel.level.pop=pixelPop, ea.samples=pixelPop$ea.samples, 
  areas=pop.matKenya$subarea, stratify.by.urban=TRUE, 
  target.pop.mat=pop.matKenya, do.fine.scale.risk=TRUE)

# get areas associated with each subarea for aggregation
tempAreasFrom = pop.matKenya$subarea
tempAreasTo = pop.matKenya$area
areas.from = sort(unique(tempAreasFrom))
areas.toI = match(areas.from, tempAreasFrom)
areas.to = tempAreasTo[areas.toI]

# do the aggregation from subareas to areas
outAreaLevel = areaPopToArea(area.level.pop=subareaPop, 
  areas.from=areas.from, areas.to=areas.to, 
  stratify.by.urban=TRUE, do.fine.scale.risk=TRUE)
} # }
```
