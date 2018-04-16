#' Function to simulate hyperpriors from an adjacency matrix.
#' 
#'
#' @param R Desired prior odds ratio. Default to 2, i.e., a 95\% prior interval for the residual odds ratios lies in the interval (R, 1/R).
#' @param nsamp Sample to simulate for scaling factor
#' @param nsamp.check Sample to simulate for checking range
#' @param Amat Adjacency matrix of the areas in the data.
#' @param nperiod numerical value of how many time periods in the data
#' @param only.iid Indicator for whether or not only IID hyperpriors are simulated
#' @references Wakefield, J. Multi-level modelling, the ecologic fallacy, and hybrid study designs. \emph{International Journal of Epidemiology}, 2009, vol. 38 (pg. 330-336).
#' @examples
#' \dontrun{
#' data(DemoMap)
#' mat <- DemoMap$Amat
#' priors <- simhyper(R = 2, nsamp = 1e+05, nsamp.check = 5000, Amat = mat)
#' }
#' @export
simhyper <- function(R = 2, nsamp = 1e+05, nsamp.check = 5000, Amat, nperiod = 6, only.iid = TRUE) {
    #################################################################### (R,1/R) is the range of the residual odds ratios gives a=d/2 where d = degrees of freedom of marginal Student’s t
    d <- 1
    a <- d/2
    p <- 0.025
    b <- (log(R))^2 * d/(2 * stats::qt(p, df = d)^2)
    a.iid <- a
    b.iid <- b
    if (!only.iid) {
    ##################################################################### Check range using simulation #
    tausamp <- stats::rgamma(nsamp, a, b)
    Usamp <- stats::rnorm(nsamp, mean = 0, sd = 1/sqrt(tausamp))
    # The adjacency matrix for 6 years
    m2 = c(1, rep(2, nperiod - 2), 1)
    # create the adjacency #
    before <- 1:(nperiod - 2)
    after <- 3:(nperiod)
    alternate <- c(rbind(before, after))
    adj2 <- c(2, alternate, nperiod - 1)
    
    make.Q <- function(num.neighbors, neighbors, omega.sq = 1) {
        n <- length(num.neighbors)
        mat <- matrix(0, ncol = n, nrow = n)
        diag(mat) <- num.neighbors
        mat[cbind(rep(1:n, num.neighbors), neighbors)] <- -1
        mat/omega.sq
    }
    vars.Q <- function(eigenvalues, eigenvectors) {
        margsum <- 0
        nloop <- length(eigenvalues) - 1
        for (i in 1:nloop) {
            ev <- eigenvectors[, i]
            margsum <- margsum + ev %*% t(ev)/eigenvalues[i]
        }
        margvars <- diag(margsum)
        margvars
    }
    # 
    sim.Q <- function(Q) {
        eigenQ <- eigen(Q)
        rankQ <- qr(Q)$rank
        sim <- as.vector(eigenQ$vectors[, 1:rankQ] %*% matrix(stats::rnorm(rep(1, rankQ), rep(0, rankQ), 1/sqrt(eigenQ$values[1:rankQ])), 
            ncol = 1))
        sim
    }
    # 
    Q <- make.Q(m2, adj2, 1)
    eigentemp <- eigen(Q)
    eigenvaluesQ <- eigentemp$values
    eigenvectorsQ <- eigentemp$vectors
    rankQ <- qr(Q)$rank  # 5
    margy <- mean(vars.Q(eigenvaluesQ, eigenvectorsQ))
    a.rw1 = a
    b.rw1 = b/margy
    taustarsamp <- stats::rgamma(nsamp.check, a.rw1, b.rw1)
    Ustarsamp <- matrix(nrow = nsamp.check, ncol = nperiod)
    for (i in 1:nsamp.check) {
        Qstar <- Q * taustarsamp[i]
        Ustarsamp[i, ] <- sim.Q(Qstar)
    }
    check.rw1 = stats::quantile(exp(Ustarsamp), p = c(0.025, 0.5, 0.975))
    
    
    # tausamp <- rgamma(nsamp,a,b) Usamp <- rnorm(nsamp,mean=0,sd=1/sqrt(tausamp)) quantile(exp(Usamp),p=c(0.025,0.5,0.975)) m2
    # <- c(2, 3, rep(4, nperiod - 4), 3, 2) # create the adjacency # before2 <- 1:(nperiod - 4) before1 <- 2:(nperiod - 3) after1
    # <- 4:(nperiod - 1) after2 <- 5:(nperiod) alternate <- c(rbind(before1, before2, after1, after2))
    
    # # adj2<-c(2, 3, # 1, 3, 4, # alternate, # nperiod - 3, nperiod - 3, nperiod, # nperiod - 2, nperiod - 1) # Q <- make.Q(m2,
    # adj2, 1) The adjacency matrix for 6 years
    
    if (nperiod > 4) {
        Q <- matrix(0, nrow = nperiod, ncol = nperiod)
        Q[1, 1:3] <- c(1, -2, 1)
        Q[2, 1:4] <- c(-2, 5, -4, 1)
        for (j in 3:(nperiod - 2)) {
            Q[j, (j - 2):(j + 2)] <- c(1, -4, 6, -4, 1)
        }
        Q[nperiod, (nperiod - 2):nperiod] <- c(1, -2, 1)
        Q[nperiod - 1, (nperiod - 3):nperiod] <- c(1, -4, 5, -2)
    } else {
        stop("RW2 prior not specified for n < 5")
    }
    vars.Q2 <- function(eigenvalues, eigenvectors, rankQ) {
        margsum <- 0
        nloop <- rankQ
        for (i in 1:nloop) {
            ev <- eigenvectors[, i]
            margsum <- margsum + ev %*% t(ev)/eigenvalues[i]
        }
        margvars <- diag(margsum)
        margvars
    }
    sim.Q <- function(Q) {
        eigenQ <- eigen(Q)
        rankQ <- qr(Q)$rank
        sim <- as.vector(eigenQ$vectors[, 1:rankQ] %*% matrix(stats::rnorm(rep(1, rankQ), rep(0, rankQ), 1/sqrt(eigenQ$values[1:rankQ])), 
            ncol = 1))
        sim
    }
    eigentemp <- eigen(Q)
    eigenvaluesQ <- eigentemp$values
    eigenvectorsQ <- eigentemp$vectors
    rankQ <- qr(Q)$rank  # 4
    margy <- mean(vars.Q2(eigenvaluesQ, eigenvectorsQ, rankQ))
    a.rw2 = a
    b.rw2 = b/margy
    taustarsamp <- stats::rgamma(nsamp.check, a.rw2, b.rw2)
    Ustarsamp <- matrix(nrow = nsamp.check, ncol = nperiod)
    for (i in 1:nsamp.check) {
        Qstar <- Q * taustarsamp[i]
        Ustarsamp[i, ] <- sim.Q(Qstar)
    }
    check.rw2 = stats::quantile(exp(Ustarsamp), p = c(0.025, 0.5, 0.975))
    
    
    #################################################################### 
    tausamp <- stats::rgamma(nsamp, a, b)
    Usamp <- stats::rnorm(nsamp, mean = 0, sd = 1/sqrt(tausamp))
    m2 <- apply(Amat, 1, sum)
    # create the adjacency list #
    nums <- c(1:dim(Amat)[1])
    adj2 <- NULL
    for (i in 1:dim(Amat)[1]) {
        adj2 <- c(adj2, nums[as.numeric(Amat[i, ]) == 1])
    }
    make.Q <- function(num.neighbors, neighbors, omega.sq = 1) {
        n <- length(num.neighbors)
        mat <- matrix(0, ncol = n, nrow = n)
        diag(mat) <- num.neighbors
        mat[cbind(rep(1:n, num.neighbors), neighbors)] <- -1
        mat/omega.sq
    }
    vars.Q3 <- function(eigenvalues, eigenvectors, rankQ) {
        margsum <- 0
        nloop <- rankQ
        for (i in 1:nloop) {
            ev <- eigenvectors[, i]
            margsum <- margsum + ev %*% t(ev)/eigenvalues[i]
        }
        margvars <- diag(margsum)
        margvars
    }
    # 
    sim.Q <- function(Q) {
        eigenQ <- eigen(Q)
        rankQ <- qr(Q)$rank
        sim <- as.vector(eigenQ$vectors[, 1:rankQ] %*% matrix(stats::rnorm(rep(1, rankQ), rep(0, rankQ), 1/sqrt(eigenQ$values[1:rankQ])), 
            ncol = 1))
        sim
    }
    # 
    Q <- make.Q(m2, adj2, 1)
    eigentemp <- eigen(Q)
    eigenvaluesQ <- eigentemp$values
    eigenvectorsQ <- eigentemp$vectors
    rankQ <- qr(Q)$rank  # 20
    margy <- mean(vars.Q3(eigenvaluesQ, eigenvectorsQ, rankQ))
    a.icar = a
    b.icar = b/margy
    taustarsamp <- stats::rgamma(nsamp.check, a.icar, b.icar)
    Ustarsamp <- matrix(nrow = nsamp.check, ncol = dim(Amat)[1])
    for (i in 1:nsamp.check) {
        Qstar <- Q * taustarsamp[i]
        Ustarsamp[i, ] <- sim.Q(Qstar)
    }
    check.icar = stats::quantile(exp(Ustarsamp), p = c(0.025, 0.5, 0.975))

    }else{
        a.rw1 <- a.rw2 <- a.icar <- a.iid
        b.rw1 <- b.rw2 <- b.icar <- b.iid
        check.rw1 <- NULL
        check.rw2 <- NULL
        check.icar <- NULL
    }


    return(list(a.iid = a.iid, b.iid = b.iid, a.rw1 = a.rw1, b.rw1 = b.rw1, check.rw1 = check.rw1, a.rw2 = a.rw2, b.rw2 = b.rw2, 
        check.rw2 = check.rw2, a.icar = a.icar, b.icar = b.icar, check.icar = check.icar))
}
