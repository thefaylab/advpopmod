# ---------- Part 2:  Kalman Filter Function -----------------------
kf <- function(obs,A,B,Q,R,n0) {
  #-- Kalman filter estimates of state vectors
  
  # Output: 
  #  Kalman filter estimates of state vector for static SSM
  #  Also returns matrix of 1 step ahead predicted states and corresponding
  #  array of 1 step ahead predicted covariance matrices
  
  # Input:
  # obs is a p x [T+1] matrix
  # A is the state transition matrix
  # Q is the state covariance matrix, n[t] = A x n[t-1] + "Q"
  # B is the observation "linkage" matrix
  # R is the obs'n covariance matrix, y[t] = B x n[t] + "R"
  # n0 is the initial state vector
  y.vec <- obs
  
  num.times <- dim(obs)[2]-1  # start with time = 0
  dim.state <- dim(Q)[1]
  filter.n <- matrix(data=0,nrow=dim.state,ncol=num.times+1)
  pred.n <- update.n <- matrix(0,nrow=dim.state,ncol=1)
  cov.pred.n <- cov.update.n <- diag(x=rep(0,dim.state),nrow=dim.state)
  
  #-- the following are used for likelihood evaluation
  pred.n.matrix <- matrix(data=0,nrow=dim.state,ncol=num.times+1)
  cov.pred.array <- array(data=0,dim=c(dim.state,dim.state,num.times+1))
  cov.filter.array <- array(data=0,dim=c(dim.state,dim.state,num.times+1))
  
  pred.n.matrix[,1] <- n0
  
  # the iterations
  update.n <- n0
  filter.n[,1] <- update.n
  identity.mat <- diag(x=1,nrow=dim.state)
  for(i in 1:num.times) {
    pred.n <- A %*% update.n
    cov.pred.n <- A %*% cov.update.n %*% t(A) + Q
    
    Kalman.gain <- cov.pred.n %*% t(B) %*% solve(B %*% cov.pred.n %*% t(B) + R)
    
    update.n <- pred.n + Kalman.gain %*% (y.vec[,i+1,drop=FALSE]-B%*%pred.n)
    cov.update.n <- (diag(x=1,nrow=dim.state)-Kalman.gain %*% B) %*% cov.pred.n
    
    filter.n[,i+1] <- update.n
    
    pred.n.matrix[,i+1] <- pred.n
    cov.pred.array[,,i+1] <- cov.pred.n
    cov.filter.array[,,i+1] <- cov.update.n
  }
  out <- list(filter.n=filter.n,pred.n.matrix=pred.n.matrix,
              cov.pred.array=cov.pred.array,cov.filter.array=cov.filter.array)
  return(out)
}


# ---------- Part 3: Negative Log-Likelihood Function  ------ 
#  For state matrix A, with independent state components.
#  Need to have mvtnorm package installed
library(package=mvtnorm)
like.fn.2state.SSM <- function(par.vals,obs,gamma) 
  {   
  n0 <- exp(par.vals[1])
  lamda <- exp(par.vals[2])
  tau <- exp(par.vals[3])
  sigma <- exp(par.vals[4])
  A=matrix(lamda,nrow=1,ncol=1)
  Q=matrix(tau^2,nrow=1,ncol=1)
  B=matrix(gamma,nrow=1,ncol=1)
  R=matrix(sigma^2,nrow=1,ncol=1)
  #obs <- matrix(c(-1,y),nrow=1,ncol=1+length(y))
  
  num.times  <- dim(obs)[2]-1  # time starts at 0
  num.states <- dim(Q)[1] 
  #A <- matrix(data=theta,nrow=num.states,ncol=num.states,byrow=FALSE)
  out <- kf(obs=obs,A=A,B=B,Q=Q,R=R,n0=n0)
  pred.n.matrix  <- out$pred.n.matrix
  cov.pred.array <- out$cov.pred.array
  
  log.l <- 0
  for(i in 1:num.times) { 
    cond.mean.vector <- B %*% pred.n.matrix[,i]
    cond.cov.matrix  <- B %*% cov.pred.array[,,i] %*% t(B) + R
    log.l.comp <- dmvnorm(x=obs[,i+1], mean=cond.mean.vector, 
                          sigma=cond.cov.matrix, log=TRUE)
    log.l      <- log.l + log.l.comp
  }
  out <- -log.l 
  return(out)
}

