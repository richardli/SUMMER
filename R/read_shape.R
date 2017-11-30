#' Function to read shape files.
#' 
#'
#' @param filepath file path for .shp files
#' @param regionnames vector of strings of final region names
#' @param data optional country summary data, for checking
#' 
#' @return A list including shape files and the adjacency matrix.
#' 
#' @examples 
#' \dontrun{
#' my_region_names <- c("central","eastern","northern","western")
#' my_fp <- "myExampleFilepath/sdr_subnational_boundaries.shp"
#' my_map <- read_shape(filepath = my_fp, regionnames = my_region_names)
#' }
#' @export
read_shape <- function(filepath, regionnames, data = NULL) {
    if (file.exists(filepath)) {
        geo <- maptools::readShapePoly(filepath)
        
        geo$NAME_final <- regionnames
        
        nb.r <- spdep::poly2nb(geo, queen = FALSE, row.names = geo$NAME_final)
        # Construct adj matrix
        mat <- spdep::nb2mat(nb.r, style = "B", zero.policy = TRUE)
        colnames(mat) <- rownames(mat)
        mat <- as.matrix(mat[1:dim(mat)[1], 1:dim(mat)[1]])
        
        # check map and data contains the same region names
        if (!is.null(data)) {
            regions_in_data <- as.character(unique(data$region))
            regions_in_map <- regionnames
            if (sum(1 - regions_in_map %in% regions_in_data) > 0) {
                stop("Exist regions in Map file, but not in Data")
            } else if (sum(1 - regions_in_data %in% c("All", regions_in_map)) > 0) {
                stop("Exist regions in Data, but not in Map file")
            }
        }
    } else {
        stop("Map file do not exist!")
    }
    return(list(geo = geo, Amat = mat))
}
