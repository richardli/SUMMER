#' Simulate spatial and temporal random effects
#' 
#' This function simulates spatial and temporal random effects with mean zero. The method is described in Algorithm 3.1 of Rue & Held 2015.
#' 
#' @param n sample size
#' @param type type of random effects: temporal (t), spatial (s), or spatial-temporal (st)
#' @param type.s type of spatial random effect, currently only ICAR is available
#' @param type.t type of temporal random effect, currently only RW1 and RW2 are available
#' @param Amat adjacency matrix for the spatial regions
#' @param n.t number of time points for the temporal random effect
#' @param scale.model logical indicator of whether to scale the random effects to have unit generalized variance. See Sørbye 2013 for more details
#' 
#' @return a matrix (for spatial or temporal) or a three-dimensional array (for spatial-temporal) of the random effects.
#' @references Rue, H., & Held, L. (2005). \emph{Gaussian Markov random fields: theory and applications}. CRC press.
#' @references Sørbye, S. H. (2013). \emph{Tutorial: Scaling IGMRF-models in R-INLA}. Department of Mathematics and Statistics, University of Tromsø.
#' @examples
#' \dontrun{
#' data(DemoMap)
#' ## Spatial random effects 
#' out <- rst(n=10000, type = "s", Amat = DemoMap$Amat)
#' # To verify the mean under the conditional specification
#' mean(out[,1] - apply(out[,c(2,3,4)], 1, mean))  
#' mean(out[,2] - apply(out[,c(1,3)], 1, mean)) 
#' mean(out[,3] - apply(out[,c(1,2,4)], 1, mean))  
#' mean(out[,4] - apply(out[,c(1,3)], 1, mean)) 
#' 
#' ## Temporal random effects (RW1)
#' out <- rst(n=1, type = "t", type.t = "RW1", n.t = 200)
#' par(mfrow = c(1,2))
#' plot(1:dim(out)[2], out, col = 1, type = "l", xlab = "Time", ylab = "Random effects")
#' # verify the first order difference is normally distributed
#' first_diff <- diff(as.numeric(out[1,]))
#' qqnorm(first_diff )	
#' abline(c(0,1))
#' 
#' ## Temporal random effects (RW1)
#' out <- rst(n=1, type = "t", type.t = "RW2", n.t = 200)
#' par(mfrow = c(1,2))
#' plot(1:dim(out)[2], out, col = 1, type = "l", xlab = "Time", ylab = "Random effects")
#' # verify the second order difference is normally distributed
#' first_diff <- diff(as.numeric(out[1,]))
#' second_diff <- diff(first_diff)
#' qqnorm(second_diff)	
#' abline(c(0,1))
#' 
#' ## Spacial-temporal random effects
#' out <- rst(n=1, type = "st", type.t = "RW1", Amat = DemoMap$Amat, n.t = 50)
#' dimnames(out)
#' par(mfrow = c(1,1))
#' plot(1:dim(out)[3], out[1,1,], col = 1,
#'  type = "l", ylim = range(out), xlab = "Time", ylab = "Random effects")
#' for(i in 2:4) lines(1:dim(out)[3], out[1,i,], col = i)
#' legend("bottomright", colnames(DemoMap$Amat), col = c(1:4), lty = rep(1,4))
#' }
#' 
#' 
rst <- function(n = 1, type = c("s", "t", "st")[1], type.s = "ICAR", type.t = c("RW1", "RW2")[2], Amat = NULL, n.t=NULL, scale.model=TRUE){

	envir = environment(sys.call()[[1]]) 
    inla.rw = utils::getFromNamespace("inla.rw", "INLA")

	if(is.null(type)){
		stop("Need to specify type of random effects.")
	}
	type <- tolower(type)

	if(is.null(Amat) && type %in% c("s", "st")){
		stop("No spatial adjacency matrix is given.")
	}


	sim.Q <- function(Q){
	  eigenQ <- eigen(Q)
	  rankQ <- qr(Q)$rank
	  sim <- as.vector(eigenQ$vectors[,1:rankQ] %*% 
	           matrix(
	             stats::rnorm(rep(1, rankQ), rep(0, rankQ), 1/sqrt(eigenQ$values[1:rankQ])),
	           ncol = 1))
	  sim
	}

	if(type == "s"){
		if(ncol(Amat) != nrow(Amat)){
			stop("Amat does not have the same number of rows and columns.")
		}
		if(sum(Amat %in% c(0, 1)) < length(Amat)){
			stop("Amat contains values other than 0 and 1.")
		}
		Q <- Amat * -1
		diag(Q) <- 0
		diag(Q) <- -1 * apply(Q, 2, sum)
		if(scale.model)	Q <- as.matrix(INLA::inla.scale.model(Q, constr = list(A=matrix(1,1,dim(Q)[1]), e=0)))

		rsample <- matrix(NA, n, dim(Amat)[1])
		for(i in 1:n) rsample[i, ] <- sim.Q(Q) 
		if(is.null(colnames(Amat))){
			id <-  paste0("s", 1:dim(Amat)[1])
		}else{
			id <- colnames(Amat)
		}
		out <- rsample
		colnames(out) <- id
 	}else if(type == "t"){
		if(is.null(n.t)){
			stop("Need to specify the number of time points.")
		}
		if(type.t == "RW1"){
			order <- 1
		}else if(type.t == "RW2"){
			order <- 2
		}else{
			stop("Need to specify the type of temporal random effects.")
		}
		if(scale.model)	Q <- inla.rw(n.t, order = order, sparse=FALSE)
		rsample <- matrix(NA, n, n.t)
		for(i in 1:n) rsample[i, ] <- sim.Q(Q) 
		id <- 1:n.t
		out <- rsample
		colnames(out) <- id
 	}else if(type == "st"){
 		if(ncol(Amat) != nrow(Amat)){
			stop("Amat does not have the same number of rows and columns.")
		}
		if(sum(Amat %in% c(0, 1)) < length(Amat)){
			stop("Amat contains values other than 0 and 1.")
		}
		Q1 <- Amat * -1
		diag(Q1) <- 0
		diag(Q1) <- -1 * apply(Q1, 2, sum)
		if(scale.model)	Q1 <- as.matrix(INLA::inla.scale.model(Q1, constr = list(A=matrix(1,1,dim(Q1)[1]), e=0)))
		if(is.null(n.t)){
			stop("Need to specify the number of time points.")
		}
		if(type.t == "RW1"){
			order <- 1
		}else if(type.t == "RW2"){
			order <- 2
		}else{
			stop("Need to specify the type of temporal random effects.")
		}
		Q2 <- inla.rw(n.t, order = order, sparse=FALSE, scale.model=scale.model)
		Q <- Q1 %x% Q2
		rsample <- array(NA, dim = c(n, dim(Amat)[1], n.t))
		for(i in 1:n){
			tmp <- sim.Q(Q)  # first space then time
			rsample[i, , ] <- tmp
		}
		idt <- 1:n.t
		if(is.null(colnames(Amat))){
			ids <- paste0("s", 1:dim(Amat)[1])
		}else{
			ids <- colnames(Amat)
		}
		out <- rsample
		dimnames(out)[[2]] <- ids
		dimnames(out)[[3]] <- idt
 	}else{
 		stop("Unknown type. Need to be one of the following: s, t, and st")
 	}
 	return(out)
}