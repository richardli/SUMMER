#' Function to map region names to a common set.
#'
#' @param data Preprocessed data
#' @param Bmat Matrix of changes. Each row corresponds to a region name possibly in the data files, and each column corresponds to a region after mapping. The values in the matrix are binary. The row names and column names need to be specified to the region names. 
#' @param regionVar String indicating the region variable. Defaults to 'region'.
#' 
#' @return Data after changing region names
#' @examples
#' 
#' # Construct a small test data
#' testdata <- data.frame(region = c("north", "south", "east",
#'  "south", "east"), index = c(1:5))
#' 
#' # Construct a changing rule: combining south and east
#' Bmat <- matrix(c(1, 0, 0, 0, 1, 1), 3, 2)
#' colnames(Bmat) <- c("north", "south and east")
#' rownames(Bmat) <- c("north", "south", "east")
#' print(Bmat)
#' 
#' # New data after transformation
#' test <- ChangeRegion(testdata, Bmat, "region")
#' print(test)
#' @export
ChangeRegion <- function(data, Bmat, regionVar = "region") {
    final_names <- colnames(Bmat)
    current_names <- rownames(Bmat)
    nc <- length(current_names)
    current_region <- as.character(data[, regionVar])
    
    # check if there are regions not contained in the names given
    missing <- which(current_region %in% current_names == FALSE)
    if (length(missing) > 0) {
        missingregion <- unique(current_region[missing])
        warning(paste("Name for regions are inconsistent. Found the following regions in data:", missingregion), immediate. = TRUE)
    }
    
    # count changes
    nchange <- ncount <- 0
    
    for (i in 1:nc) {
        tmp <- Bmat[i, ]
        if (sum(tmp) > 1) {
            stop(paste("more than one region to map to: ", current_names[i]))
        } else if (sum(tmp) == 1) {
            which <- which(current_region == current_names[i])
            to <- final_names[which(tmp == 1)]
            if (current_names[i] != to && length(which) > 0) {
                nchange <- nchange + 1
                ncount <- ncount + length(which)
            }
            current_region[which] <- to
        } else if (sum(tmp) == 0) {
            # nchange <- nchange + 1 ncount <- ncount + length(which(current_region == current_names[i]))
        }
    }
    cat(paste(nchange, "names changed, in total", ncount, "rows in data changed\n"))
    data[, regionVar] = current_region
    return(data)
}
