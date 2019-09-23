#' New random walk 1 and 2 models for m-year to period random effects
#' 
#' @param cmd list of model components
#' @param theta log precision
#' 
rw.new.pc = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
    ## assume 'tau', 'order', 'n' and 'm' 'n' is the dim of RW and 'm' is the aggregated length,
    ## averaging over n/m variables, non-overlapping
    
    ## the environment of this function which holds the variables and we can store 'my.cache'
    ## there.
    envir = environment(sys.call()[[1]]) 
    if(is.null(envir)) envir <- environment()
     n = eval(n, environment()) 
     m = eval(m, environment()) 
     order = eval(order, environment()) 
     tau = eval(tau, environment()) 
     alpha0 = eval(alpha0, environment()) 
     u0 = eval(u0, environment())

    inla.rw = utils::getFromNamespace("inla.rw", "INLA")
    inla.ginv = utils::getFromNamespace("inla.ginv", "INLA")
    
    if (!exists("my.cache", envir = envir, mode = "list")) {
      nn = n %/% m
      stopifnot (nn == as.integer(n/m))
      R = inla.rw(n, order = order,  scale.model=TRUE, sparse=TRUE)
      A = matrix(0, nn, n)
      j = 1
      for(i in 1:nn) {
        A[i, j:(j+m-1)] = 1/m
        j = j + m
      }
      A =  INLA::inla.as.sparse(A)
      D = Matrix::Diagonal(nn, x=1)
      assign("my.cache", list(R=R, A=A, D=D, nn=nn), envir = envir)
    }else{
      my.cache <- eval(my.cache, envir = envir)
    } 
    
    interpret.theta = function() {
      return(list(kappa = exp(theta[1L])))
    }
    
    graph = function() {
      return (Q())
    }
    
    Q = function() {
      QQ = rbind(cbind(p$kappa * my.cache$R + tau * t(my.cache$A) %*% my.cache$A,
                       -tau * t(my.cache$A)),
                 cbind(-tau * my.cache$A, tau * my.cache$D))
      return(QQ)
    }
    
    mu = function() {
      return(numeric(0))
    }
    
    log.norm.const = function() {
      val = (n-order) * (-0.5 * log(2 * pi) + 0.5 * log(p$kappa)) +
        (my.cache$nn * (-0.5 * log(2 * pi) + 0.5 * log(tau)))
      return(val)
    }
    
    log.prior = function() {
      # val = dgamma(p$kappa, shape = shape0, rate = rate0, log = TRUE) + theta[1]
      val = INLA::inla.pc.dprec(p$kappa, u = u0, alpha = alpha0, log = TRUE) + theta[1]
      return(val)
    }
    
    initial = function() {
      return(4)
    }
    
    quit = function() {
      return(invisible())
    }
    
    ## as some calls to this function does not define 'theta',  its convenient to have to
    ## defined still (like in the graph-function)
    if (is.null(theta))
      theta = initial()
    
    p = interpret.theta()
    val = do.call(match.arg(cmd), args = list())
    return(val)
   }   
  
#' New random IID models for m-year to period random effects
#' 
#' @param cmd list of model components
#' @param theta log precision
#' 
    iid.new.pc = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
    
    envir = environment(sys.call()[[1]]) 
    if(is.null(envir)) envir <- environment()
     n = eval(n, environment()) 
     m = eval(m, environment()) 
     order = eval(order, environment()) 
     tau = eval(tau, environment()) 
     alpha0 = eval(alpha0, environment()) 
     u0 = eval(u0, environment())


    if (!exists("my.cache", envir = envir, mode = "list")) {
      nn = n %/% m
      stopifnot (nn == as.integer(n/m))
      R = Matrix::Diagonal(n, x = rep(1, n))
      A = matrix(0, nn, n)
      j = 1
      for(i in 1:nn) {
        A[i, j:(j+m-1)] = 1/m
        j = j + m
      }
      A = INLA::inla.as.sparse(A)
      D = Matrix::Diagonal(nn, x=1)
      assign("my.cache", list(R=R, A=A, D=D, nn=nn), envir = envir)
    }else{
        my.cache <- eval(my.cache, envir = envir)
    } 
    
    interpret.theta = function() {
      return(list(kappa = exp(theta[1L])))
    }
    
    graph = function() {
      return (Q())
    }
    
    Q = function() {
      QQ = rbind(cbind(p$kappa * my.cache$R + tau * t(my.cache$A) %*% my.cache$A,
                       -tau * t(my.cache$A)),
                 cbind(-tau * my.cache$A, tau * my.cache$D))
      return(QQ)
    }
    
    mu = function() {
      return(numeric(0))
    }
    
    log.norm.const = function() {
      val = (n * (-0.5 * log(2 * pi) + 0.5 * log(p$kappa)) +
        (my.cache$nn * (-0.5 * log(2 * pi) + 0.5 * log(tau))))
      return(val)
    }
    
    log.prior = function() {
      # val = dgamma(p$kappa, shape = shape0, rate = rate0, log = TRUE) + theta[1]
      val = INLA::inla.pc.dprec(p$kappa, u = u0, alpha = alpha0, log = TRUE) + theta[1]
      return(val)
    }
    
    initial = function() {
      return(4)
    }
    
    quit = function() {
      return(invisible())
    }
    
    ## as some calls to this function does not define 'theta',  its convenient to have to
    ## defined still (like in the graph-function)
    if (is.null(theta))
      theta = initial()
    
    p = interpret.theta()
    val = do.call(match.arg(cmd), args = list())
    return(val)
  }  

  
#' New Type I to IV space time interaction models for m-year to period random effects
#' 
#' @param cmd list of model components
#' @param theta log precision
#' 
  st.new.pc = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
  
  envir = environment(sys.call()[[1]]) 
        if(is.null(envir)) envir <- environment()
  n = eval(n, environment()) 
  m = eval(m, environment()) 
  order = eval(order, environment()) 
  tau = eval(tau, environment()) 
  alpha0 = eval(alpha0, environment()) 
  u0 = eval(u0, environment())
  S = eval(S, environment())
  Amat = eval(Amat, environment())
  type = eval(type, environment())

  inla.rw = utils::getFromNamespace("inla.rw", "INLA")
  inla.ginv = utils::getFromNamespace("inla.ginv", "INLA")
  # The new structure takes the following order
  # (x_11, ..., x_1T, ..., x_S1, ..., x_ST, xx_11, ..., xx_1t, ..., xx_S1, ..., xx_St)
  #  x_ij : random effect of region i, year j 
  # xx_ik : random effect of region i, period k

  if (!exists("my.cache", envir = envir, mode = "list")) {
    nn = n %/% m
    stopifnot (nn == as.integer(n/m))
    R1 = Matrix::Diagonal(n, x = rep(1, n))
    R2 = inla.rw(n, order = order, scale.model=TRUE, sparse=TRUE)
    R3 = Matrix::Diagonal(S, x = rep(1, S))
    R4 = Amat
    diag(R4) <- 0
    diag <- apply(R4, 1, sum)
    R4[R4 != 0] <- -1
    diag(R4) <- diag

    R4 <- INLA::inla.scale.model(R4, constr = list(A=matrix(1,1,dim(R4)[1]), e=0))
    # both independent
    if(type == 1){
        R <- R3 %x% R1
    # AR * independent    
    }else if(type == 2){
        R <- R3 %x% R2
    # independent * besag    
    }else if(type == 3){
        R <- R4 %x% R1
    # AR * besag
    }else if(type == 4){
        R <- R4 %x% R2
    }

    A = matrix(0, nn*S, n*S)
    j = 1
    for(i in 1:(nn*S)) {
      A[i, j:(j+m-1)] = 1/m
      j = j + m
    }
    A = INLA::inla.as.sparse(A)
    D = Matrix::Diagonal(nn*S, x=1)
    assign("my.cache", list(R=INLA::inla.as.sparse(R), A=A, D=D, nn=nn), envir = envir)
  }else{
      my.cache <- eval(my.cache, envir = envir)
  } 
  
  interpret.theta = function() {
    return(list(kappa = exp(theta[1L])))
  }
  
  graph = function() {
    return (Q())
  }
  
  Q = function() {
    QQ = rbind(cbind(p$kappa * my.cache$R + tau * t(my.cache$A) %*% my.cache$A,
                       -tau * t(my.cache$A)),
                 cbind(-tau * my.cache$A, tau * my.cache$D))
    return(QQ)
  }
  
  mu = function() {
    return(numeric(0))
  }
  ## Type I   : S * n
  ## Type II  : S * (n - order)
  ## Type III : (S-1) * n 
  ## Type IV  : (S-1) * (n - order)
  log.norm.const = function() {
    df <- S * n
    if(type == 2){
      df <- S * (n - order)
    }else if(type == 3){
      df <- (S-1) * n
    }else if(type == 4){
      df <- (S-1) * (n - order)
    }
    val = (df * (-0.5 * log(2 * pi) + 0.5 * log(p$kappa)) +
      (S * my.cache$nn * (-0.5 * log(2 * pi) + 0.5 * log(tau))))
    return(val)
  }
  
  log.prior = function() {
    # val = dgamma(p$kappa, shape = shape0, rate = rate0, log = TRUE) + theta[1]
      val = INLA::inla.pc.dprec(p$kappa, u = u0, alpha = alpha0, log = TRUE) + theta[1]
    return(val)
  }
  
  initial = function() {
    return(4)
  }
  
  quit = function() {
    return(invisible())
  }
  
  ## as some calls to this function does not define 'theta',  its convenient to have to
  ## defined still (like in the graph-function)
  if (is.null(theta))
    theta = initial()
  
  p = interpret.theta()
  val = do.call(match.arg(cmd), args = list())
  return(val)
}  
