#' Aggregate person-month data into counts and totals by groups.
#' 
#' @param data dataset in person-month format
#' @param variables a character vector of the variables to aggregate
#' @param by a character vector of columns that specifies which groups to aggregate by.
#' @param ignore list of conditions not to impute 0. If left unspecified, any group levels not in the data will be imputed to have 0 counts. 
#' @param addtotal logical indicator of whether to add a column of group total counts.
#' @param drop logical indicator of whether to drop all rows with total = 0.
#' 
#' @return data.frame of the ggregated counts. 
#' @author Zehang Richard Li
#' @examples
#' 
#'  
#' # a toy dataset with 4 time periods but one missing in data
#' timelist <- factor(1:4)
#' data = data.frame(died = c(0,0,0,1,1,0,0), 
#' 					area = c(rep(c("A", "B"), 3), "A"), 
#' 					time = timelist[c(1,1,2,3,3,3,3)])
#' data
#' # without ignore argument, all levels will be imputed
#' getCounts(data, variables = "died", by = c("area", "time"))
#'
#' # ignoring time = 4, the ignored level will not be imputed (but still in the output)
#' getCounts(data, variables = "died", by = c("area", "time"), 
#' 			ignore = list("time"=c(4)) )
#'
#'  
#' @export


getCounts <- function(data, variables, by, ignore = NULL, addtotal = TRUE, drop=TRUE){
	if(addtotal){
		data$total <- 1
		variables <- c(variables, "total")
	}
	data <- data[, c(variables, by)]
	formula <- as.formula(paste0(".~", paste(by, collapse = " + ")))
	out <- aggregate(formula, data = data, FUN = sum, drop = drop)
	for(v in variables){
		out[is.na(out[, v]), v] <- 0
	}
	if(!is.null(ignore)){
		for(name in names(ignore)){
			out[out[, name] %in% ignore[[name]], variables] <- NA
		}
	}
	return(out)
}