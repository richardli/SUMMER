#' Map GPS points to polygon regions 
#' 
#' @param data point data with two columns of GPS locations.
#' @param geo SpatialPolygonsDataFrame of the map
#' @param long column name for longitudinal coordinate in the data
#' @param lat column name for latitude coordinate in the data
#' @param names  character vector of region ids to be added to the neighbours list
#' 
#' @return Spatial djacency matrix.
#' @examples
#' data(DemoMap) 
#' dat <- data.frame(ID = c(1,2,3), lon = c(32.2, 33.7, 33), lat = c(0.1, 0.9, 2.8))
#' dat2 <- mapPoints(dat, DemoMap$geo, long = "lon", lat = "lat", names = "REGNAME")
#' dat2
#'  
#' @export


mapPoints <- function(data, geo, long, lat, names){
 	 gridP = data.frame(Longitude = data[, long],
                    	 Latitude = data[, lat])
  
 	 sp::coordinates(gridP) = ~ Longitude + Latitude
	 sp::proj4string(gridP) = sp::proj4string(geo)
	 new = data.frame(sp::over(gridP, geo))
	 new <- new[, names, drop=FALSE]
	 colnames(new) <- names
	 data <- cbind(data, new)

	 return(data)
}