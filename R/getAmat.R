#' Extract adjacency matrix from the map 
#' 
#' @param geo SpatialPolygonsDataFrame of the map
#' @param names  character vector of region ids to be added to the neighbours list
#' 
#' @return Spatial djacency matrix.
#' @examples
#' data(DemoMap) 
#' mat <- getAmat(geo = DemoMap$geo, names = DemoMap$geo$REGNAME)
#' mat
#' DemoMap$Amat
#'  
#' @export


getAmat <- function(geo, names){
	nb.r <- spdep::poly2nb(geo, queen=F, row.names = names)
	mat <- spdep::nb2mat(nb.r, style="B",zero.policy=TRUE)
	regions <- colnames(mat) <- rownames(mat) 
	mat <- as.matrix(mat[1:dim(mat)[1], 1:dim(mat)[1]])
	return(mat)
}