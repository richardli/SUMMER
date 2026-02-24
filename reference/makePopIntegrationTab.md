# Generating pixellated populations, and population frames

**\[experimental\]**

## Usage

``` r
makePopIntegrationTab(
  km.res = 5,
  pop,
  domain.map.dat,
  east.lim,
  north.lim,
  map.projection,
  area.map.dat,
  subarea.map.dat,
  areaNameVar = "NAME_1",
  subareaNameVar = "NAME_2",
  poppa = NULL,
  poppsub = NULL,
  stratify.by.urban = TRUE,
  areapa = NULL,
  areapsub = NULL,
  custom.subset.polygons = NULL,
  area.polygon.subset.I = NULL,
  subarea.polygon.subset.I = NULL,
  mean.neighbor = 50,
  delta = 0.1,
  return.popp.tables = FALSE,
  set.na.to.zero = TRUE,
  fix.zero.pop.density.subareas = FALSE,
  extract.method = "bilinear"
)

getPoppsub(
  km.res = 1,
  pop,
  domain.map.dat,
  east.lim,
  north.lim,
  map.projection,
  poppa,
  areapa = NULL,
  areapsub,
  subarea.map.dat,
  subareaNameVar = "NAME_2",
  stratify.by.urban = TRUE,
  area.map.dat = NULL,
  areaNameVar = "NAME_1",
  area.polygon.subset.I = NULL,
  subarea.polygon.subset.I = NULL,
  custom.subset.polygons = NULL,
  mean.neighbor = 50,
  delta = 0.1,
  set.na.to.zero = TRUE,
  fix.zero.pop.density.subareas = FALSE
)

adjustPopMat(
  pop.mat,
  poppa.target = NULL,
  adjust.by = c("area", "subarea"),
  stratify.by.urban = TRUE
)
```

## Arguments

- km.res:

  The resolution of the pixelated grid in km

- pop:

  Population density raster

- domain.map.dat:

  A shapefile representing the full spatial domain (e.g. country)

- east.lim:

  Range in km easting over the spatial domain under the input projection

- north.lim:

  Range in km northing over the spatial domain under the input
  projection

- map.projection:

  A projection function taking longitude and latitude and returning
  easting and northing in km. Or the inverse if inverse is set to TRUE.
  For example,
  [`projKenya`](https://richardli.github.io/SUMMER/reference/projKenya.md).
  Check https://epsg.io/ for example for best projection EPSG codes for
  specific countries

- area.map.dat:

  SpatialPolygonsDataFrame object with area level map information

- subarea.map.dat:

  SpatialPolygonsDataFrame object with subarea level map information

- areaNameVar:

  The name of the area variable associated with `area.map.dat@data` and
  `subarea.map.dat@data`

- subareaNameVar:

  The name of the subarea variable associated with
  `subarea.map.dat@data`

- poppa:

  data.frame of population per area separated by urban/rural. If
  `poppsub` is not included, this is used for normalization of
  populations associated with population integration points. Contains
  variables:

  area

  :   name of area

  popUrb

  :   total urban (general) population of area

  popRur

  :   total rural (general) population of area

  popTotal

  :   total (general) population of area

  pctUrb

  :   percentage of population in the area that is urban (between 0 and
      100)

- poppsub:

  data.frame of population per subarea separated by urban/rural using
  for population normalization or urbanicity classification. Often based
  on extra fine scale population density grid. Has variables:

  subarea

  :   name of subarea

  area

  :   name of area

  popUrb

  :   total urban (general) population of subarea

  popRur

  :   total rural (general) population of subarea

  popTotal

  :   total (general) population of subarea

  pctUrb

  :   percentage of population in the subarea that is urban (between 0
      and 100)

- stratify.by.urban:

  Whether to stratify the pixellated grid by urban/rural. If TRUE,
  renormalizes population densities within areas or subareas crossed
  with urban/rural

- areapa:

  A list with variables:

  area

  :   name of area

  spatialArea

  :   spatial area of the subarea (e.g. in km^2)

- areapsub:

  A list with variables:

  subarea

  :   name of subarea

  spatialArea

  :   spatial area of the subarea (e.g. in km^2)

- custom.subset.polygons:

  'SpatialPolygonsDataFrame' or 'SpatialPolygons' object to subset the
  grid over. This option can help reduce computation time relative to
  constructing the whole grid and subsetting afterwards.
  `area.polygon.subset.I` or `subarea.polygon.subset.I` can be used when
  subsetting by areas or subareas in `area.map.dat` or
  `subarea.map.dat`. Must be in latitude/longitude projection
  "EPSG:4326"

- area.polygon.subset.I:

  Index in area.map.dat for a specific area to subset the grid over.
  This option can help reduce computation time relative to constructing
  the whole grid and subsetting afterwards

- subarea.polygon.subset.I:

  FOR EXPERIMENTAL PURPOSES ONLY. Index in subarea.map.dat for a
  specific area to subset the grid over. This option can help reduce
  computation time relative to constructing the whole grid and
  subsetting afterwards

- mean.neighbor:

  For determining what area or subarea points are nearest to if they do
  not directly fall into an area. See
  [`fields.rdist.near`](https://rdrr.io/pkg/fields/man/rdist.html) for
  details.

- delta:

  For determining what area or subarea points are nearest to if they do
  not directly fall into an area. See
  [`fields.rdist.near`](https://rdrr.io/pkg/fields/man/rdist.html) for
  details.

- return.popp.tables:

  If TRUE, poppa and poppsub will be calculated based on the generated
  population integration matrix and input area/subarea map data

- set.na.to.zero:

  If TRUE, sets NA populations to 0.

- fix.zero.pop.density.subareas:

  If TRUE, if population density in a subarea is estimated to be zero,
  but the total population in the subarea is nonzero, population is
  filled into the area uniformly

- extract.method:

  Either 'bilinear' or 'simple'. see `method` from
  [`extract`](https://rspatial.github.io/terra/reference/extract.html)

- pop.mat:

  Pixellated grid data frame with variables `area` and `pop` such as
  that generated by `makePopIntegrationTab`

- poppa.target:

  Target population per area stratified by urban rural. Same format as
  poppa

- adjust.by:

  Whether to adjust population density by the `area` or `subarea` level

## Details

Functions for generating pixellated population information and
population frames at the `area` and `subarea` levels. The `area` and
`subarea` levels can be thought of as big regions and little regions,
where areas can be partitioned into unique sets of subareas. For
example, Admin-1 and Admin-2 areas might be areas and subareas
respectively. The population totals are either tabulated at the area x
urban/rural level, the subarea x urban/rural level, or at the pixel
level of a specified resolution. Totals are calculated using population
density information, shapefiles, and, possibly, preexisting population
frames at different areal levels. Note that area names should each be
unique, and similarly for subarea names.

## Functions

- `makePopIntegrationTab()`: Generate pixellated `grid` of coordinates
  (both longitude/latitude and east/north) over spatial domain of the
  given resolution with associated population totals, areas, subareas,
  and urban/rural levels. For very small areas that might not otherwise
  have a grid point in them, a custom integration point is added at
  their centroid. Sets urbanicity classifications by thresholding input
  population density raster using area and subarea population tables,
  and generates area and subarea population tables from population
  density information if not already given. Can be used for integrating
  predictions from the given coordinates to area and subarea levels
  using population weights.

- `getPoppsub()`: Generate table of estimates of population totals per
  subarea x urban/rural combination based on population density raster
  at `kmres` resolution "grid", including custom integration points for
  any subarea too small to include grid points at their centroids.

- `adjustPopMat()`: Adjust population densities in grid based on a
  population frame.

## See also

[`setThresholdsByRegion`](https://richardli.github.io/SUMMER/reference/setThresholdsByRegion.md),
[`poppRegionFromPopMat`](https://richardli.github.io/SUMMER/reference/poppRegionFromPopMat.md),
[`simPopSPDE`](https://richardli.github.io/SUMMER/reference/simPop.md),
[`simPopCustom`](https://richardli.github.io/SUMMER/reference/simPop.md)

## Author

John Paige

## Examples

``` r
if (FALSE) { # \dontrun{

library(sp)
library(sf)
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
require(fields)
# NOTE: the following function will typically take ~15-20 minutes. Can speed up 
#       by setting km.res to be higher, but we recommend fine resolution for 
#       this step, since it only needs to be done once. Instead of running 
#       the code in the following if(FALSE) section, 
#       you can simply run data(kenyaPopulationData)
if(FALSE){
  system.time(poppsubKenya <- getPoppsub(
    km.res=1, pop=pop, domain.map.dat=adm0,
    east.lim=east.lim, north.lim=north.lim, map.projection=projKenya,
    poppa = poppaKenya, areapa=areapaKenya, areapsub=areapsubKenya, 
    area.map.dat=adm1, subarea.map.dat=adm2, 
    areaNameVar = "NAME_1", subareaNameVar="NAME_2"))
}
data(kenyaPopulationData)

# Now generate a general population integration table at 5km resolution, 
# based on subarea (Admin-2) x urban/rural population totals. This takes 
# ~1 minute
system.time(pop.matKenya <- makePopIntegrationTab(
  km.res=5, pop=pop, domain.map.dat=adm0,
  east.lim=east.lim, north.lim=north.lim, map.projection=projKenya,
  poppa = poppaKenya, poppsub=poppsubKenya, 
  area.map.dat = adm1, subarea.map.dat = adm2,
  areaNameVar = "NAME_1", subareaNameVar="NAME_2"))

## Adjust pop.mat to be target (neonatal) rather than general population density. First
## create the target population frame
## (these numbers are based on IPUMS microcensus data)
mothersPerHouseholdUrb = 0.3497151
childrenPerMotherUrb = 1.295917
mothersPerHouseholdRur = 0.4787696
childrenPerMotherRur = 1.455222
targetPopPerStratumUrban = easpaKenya$HHUrb * mothersPerHouseholdUrb * childrenPerMotherUrb
targetPopPerStratumRural = easpaKenya$HHRur * mothersPerHouseholdRur * childrenPerMotherRur
easpaKenyaNeonatal = easpaKenya
easpaKenyaNeonatal$popUrb = targetPopPerStratumUrban
easpaKenyaNeonatal$popRur = targetPopPerStratumRural
easpaKenyaNeonatal$popTotal = easpaKenyaNeonatal$popUrb + easpaKenyaNeonatal$popRur
easpaKenyaNeonatal$pctUrb = 100 * easpaKenyaNeonatal$popUrb / easpaKenyaNeonatal$popTotal
easpaKenyaNeonatal$pctTotal = 
  100 * easpaKenyaNeonatal$popTotal / sum(easpaKenyaNeonatal$popTotal)

# Generate the target population density by scaling the current population density grid 
# at the Admin1 x urban/rural level
pop.matKenyaNeonatal = adjustPopMat(pop.matKenya, easpaKenyaNeonatal)

# Generate neonatal population table from the neonatal population integration matrix.
# This is technically not necessary for population simulation purposes, but is here 
# for illustrative purposes
poppsubKenyaNeonatal = poppRegionFromPopMat(pop.matKenyaNeonatal, pop.matKenyaNeonatal$subarea)
poppsubKenyaNeonatal = cbind(subarea=poppsubKenyaNeonatal$region, 
                             area=adm2@data$NAME_1[match(poppsubKenyaNeonatal$region, 
                               adm2@data$NAME_2)], 
                             poppsubKenyaNeonatal[,-1])
print(head(poppsubKenyaNeonatal))
} # }
```
