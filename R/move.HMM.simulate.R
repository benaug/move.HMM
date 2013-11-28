#'Simulate data from a hidden markov model
#'
#'This function simulates data from a user-specified hidden markov model 
#'
#'@param dists A length ndist vector of distributions from the following list:
#'weibull, gamma, exponential, normal, lognormal, lnorm3, posnorm,
#'invgamma, rayleigh, f, ncf, dagum, frechet, beta, binom, poisson, nbinom,
#'zapois, wrpcauchy, wrpnorm
#'@param params a list of length ndist+1 containing matrices of starting parameter
#'values.  The first element of the list must be the starting values for the
#'transition matrix.  If any distributions only have 1 parameter, the list entry
#'must be a nstates x 1 matrix.
#'@param n The length of the hidden markov model to be simulated
#'@param delta an optional vector of starting state probabilities.  If no vector
#'is supplied, the stationary distribution is used.
#'@include Distributions.R
#'@return The observations from the simulated HMM and the true state sequence
#'@export
#'@examples \dontrun{
#'#2 states, 2 dist-lognorm, wrapped normal
#'lmean=c(-3,-1) #meanlog parameters
#'sd=c(1,1) #sdlog parameters
#'rho<-c(0.2,0.3) # wrapped normal concentration parameters
#'mu<-c(pi,0) # wrapped normal mean parameters
#'gamma0=matrix(c(0.6,0.4,0.2,0.8),byrow=T,nrow=2)
#'
#'dists=c("lognormal","wrpnorm")
#'nstates=2
#'turn=c(1,2)
#'params=vector("list",3)
#'params[[1]]=gamma0
#'params[[2]]=cbind(lmean,sd)
#'params[[3]]=cbind(mu,rho)
#'obs=move.HMM.simulate(dists,params,1000)
#'#see more examples in move.HMM.mle
#'}
#'@export
#'
move.HMM.simulate <- function(dists,params,n,delta=NULL){
  if(!all(unlist(lapply(params,is.matrix))))stop("argument 'params' must contain nstate x nparam matrices")
  if(any((rowSums(params[[1]])-1)>1e10))stop("Transition matrix rows should sum to 1")
  if(!all(nrow(params[[1]])-unlist(lapply(params,nrow))==0))stop("All parameter matrices must have the same number of rows")
  nstates <- ncol(params[[1]])
  if(!is.null(delta)){
    if(!is.null(dim(delta)))stop("User-specified delta must be a vector")
    if((sum(delta)!=1))stop("User-specified delta must sum to 1")
    if(length(delta)!=nstates)stop("User-specified delta must be of length 'nstates'")
  }
  out <- Distributions(dists,nstates)
  generate <- out[[5]]
  ndists <- length(dists)
  obs <- matrix(NA,nrow=n,ncol=ndists)
  #If no initial distribution given, use stationary
  if(is.null(delta))delta <-solve(t(diag(nstates)-params[[1]] +1),rep(1,nstates))
  statevect <- 1:nstates
  state <- numeric(n)
  state [1] <- sample(statevect ,1,prob=delta)
  for (i in 2:n){
    state[i] <-sample(statevect ,1,prob=params[[1]][state[i-1] ,])
  }
  for(i in 1:ndists){
    nparam <- max(1,ncol(params[[i+1]]))
    if(nparam==2){
      obs[,i] <- generate[[i]](n,params[[i+1]][state,1],params[[i+1]][state,2])
    }else if(nparam==1){
      obs[,i] <- generate[[i]](n,params[[i+1]][state,1])
    }else if(nparam==3){
      obs[,i] <- generate[[i]](n,params[[i+1]][state,1],params[[i+1]][state,2],params[[i+1]][state,3])
    }
  }
  list(obs=obs,states=state)
}
