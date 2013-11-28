#'Simulate data from a hidden semi markov model
#'
#'This function simulates data from a user-specified hidden semi markov model 
#'
#'@param dists A vector of distributions of length ndist+1, with the first distribution chosen from: shiftpois, shiftnegbin  The following distribution(s)
#'are  from the following list: weibull, gamma, exponential, normal, lognormal, lnorm3, posnorm,
#'invgamma, rayleigh, f, ncf, dagum, frechet, beta, binom, poisson, nbinom,
#'zapois, wrpcauchy, wrpnorm.
#'@param params a list of starting parameter values.  If nstates=2, params is of length ndist+1
#'The first element of the list must be the starting values for the
#'dwell time distribution.  If any distributions only have 1 parameter, the list entry
#'must be a nstates x 1 matrix. If nstates>2, params is of lenght ndist+2, with the first element
#'being the starting values for the t.p.m., the second element being the starting values for the
#'dwell time distribution, and the following element(s) being the starting values for the observation
#'distributions
#'@param n The length of the hidden markov model to be simulated
#'@param nstates The number of hidden states
#'@param delta an optional vector of starting state probabilities.  If no vector
#'is supplied, the stationary distribution is used.
#'@include Distributions.R
#'@return The observations from the simulated HMM and the true state sequence
#'@export
move.HSMM.simulate=function(dists,params,n,nstates,delta=NULL){
  if(!all(unlist(lapply(params,is.matrix))))stop("argument 'params' must contain nstate x nparam matrices")
  if(!all(nrow(params[[1]])-unlist(lapply(params,nrow))==0))stop("All parameter matrices must have the same number of rows")
  #make sure first matrix is valid tpm for nstates>2
   if(nstates>2){
     if(!all(diag(params[[1]])==0))stop("Diagonal elements of t.p.m. must be 0")
     if(!all((rowSums(params[[1]])-1)<1e-12))stop("Rows of t.p.m. must sum to 1")
   }
  if(nstates==2){
    #make sure no tpm is entered for nstates=2
    if(all(diag(params[[1]])==0))stop("Do not enter t.p.m if nstates=2")
  }
  if(!is.null(delta)){
    if(!is.null(dim(delta)))stop("User-specified delta must be a vector")
    if((sum(delta)!=1))stop("User-specified delta must sum to 1")
    if(length(delta)!=nstates)stop("User-specified delta must be of length 'nstates'")
  }
  out=Distributions(dists,nstates)
  generate=out[[5]]
  ndists=length(dists)-1
  
  #Make tpm if states = 2
  if(nstates==2){
    gamma=matrix(c(0,1,1,0),byrow=T,nrow=2)
    params=c(list(gamma),params)
  }
  gamma=params[[1]]  
  obs=matrix(NA,nrow=n,ncol=ndists)
  #Initial distribution is stationary distribution
  if(is.null(delta)){
    delta=solve(t(diag(nrow(gamma))-gamma+1),rep(1,nrow(gamma)))
  }
  nparam=max(1,ncol(params[[2]]))
  statevect <- 1:nstates
  state <- numeric(n)
  
  #Simulate states
  if (length(dim(delta))==2) delta <- delta[,1]
  start=sample(statevect ,1,prob=delta)
  if(nparam==2){
    dwell=generate[[1]](1,params[[2]][start,1],params[[2]][start,2])
  }else if(nparam==1){
    dwell=generate[[1]](1,params[[2]][start,1])
  }else if(nparam==3){
    dwell=generate[[1]](1,params[[2]][start,1],params[[2]][start,2],params[[2]][start,2])
  }
  state[1:dwell]=rep(start,dwell)
  idx=dwell
  while(idx<n){
    #update state
    current=sample(statevect ,1,prob=gamma[state[idx],])
    #how long do we dwell?
    if(nparam==2){
      dwell=generate[[1]](1,params[[2]][current,1],params[[2]][current,2])
    }else if(nparam==1){
      dwell=generate[[1]](1,params[[2]][current,1])
    }else if(nparam==3){
      dwell=generate[[1]](1,params[[2]][current,1],params[[2]][current,2],params[[2]][current,2])
    }
    idx2=idx+dwell
    state[(idx+1):(idx2)]=rep(current,dwell)
    idx=idx2
  }
  state=state[1:n]
  
  #Generate observations given state
  for(i in 1:ndists){
    nparam=max(1,ncol(params[[i+2]]))
    if(nparam==2){
      obs[,i]=generate[[i+1]](n,params[[i+2]][state,1],params[[i+2]][state,2])
    }else if(nparam==1){
      obs[,i]=generate[[i+1]](n,params[[i+2]][state,1])
    }else if(nparam==3){
      obs[,i]=generate[[i+1]](n,params[[i+2]][state,1],params[[i+2]][state,2],params[[i+2]][state,3])
    }
  }
  list(obs=obs,states=state)
}
