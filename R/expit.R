#' Expit transformation
#' 
#' @param x data
#' 
#' @return expit of x
#' 
#' @examples 
#' x <- .5
#' expit(x)
#' 
#' @export
expit <- function(x) {
    return(exp(x)/(1 + exp(x)))
}

#' Logit transformation 
#' 
#' @param x data
#' 
#' @return logit of x
#' 
#' @examples
#' x <- .5
#' logit(x)
#' 
#' @export
logit <- function(x) {
    log(x/(1 - x))
}
