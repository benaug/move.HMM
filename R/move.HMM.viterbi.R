#'Assign states using the Viterbi algorithm
#'
#'This function, modified from Zucchini and MacDonald (2009), assigns states
#'to observations using the Viterbi algorithm.  It takes as input a 
#'move.HMM object and an optional vector containing the starting state 
#'probabilities.
#'
#'@param move.HMM A move.HMM object containing a fitted HMM model.
#'@param delta An optional vector of starting state probabilities.  If no vector
#'is supplied, the stationary distribution is used.
#'@include Distributions.R
#'@return A vector of state assignments.
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
#'obs=move.HMM.simulate(dists,params,1000)$obs
#'turn=c(1,2)
#'move.HMM=move.HMM.mle(obs,dists,params,stepm=35,iterlim=100,turn=turn)
#'#get Viterbi state assignments
#'move.HMM.viterbi(move.HMM)
#'}
#'@export
#'
move.HMM.viterbi <-function(move.HMM,delta=NULL){
  obs <- move.HMM$obs
  params <- move.HMM$params
  nstates <- move.HMM$nstates
  dists <- move.HMM$dists
  out <- Distributions(dists,nstates)
  PDFs <- out[[3]]
  #If no initial distribution given, use stationary distribution
  ## (first column = parameter estimates only)
  if(is.null(delta))delta=move.HMM$delta[,1]
  n <- nrow(obs)
  allprobs <- matrix(rep(1,nstates*n),nrow=n)#f(y_t|s_t=k)
  for (k in 1:n){
    for (j in 1:nstates){
      for(i in 1:length(PDFs)){
        nparam=max(1,ncol(params[[i+1]]))
        argList <- c(list(obs[k,i]),
                          lapply(1:nparam,
                                 function(m) params[[i+1]][j,m]))
        pdfvals <- do.call(PDFs[[i]],argList)
        allprobs[k,j] <- allprobs[k,j]*ifelse(is.na(obs[k,i]),1,
                                              pdfvals)
      } #i index
    } # j index
  } # k index
  xi <- matrix(0,n,nstates)
  foo <- delta*allprobs [1,]
  xi[1,] <- foo/sum(foo)
  for (i in 2:n){
    foo <- apply(xi[i-1,]*params[[1]] ,2,max)*allprobs[i,]
    xi[i,] <- foo/sum(foo)
  }
  iv <-numeric(n)
  iv[n] <-which.max(xi[n,])
  for (i in (n-1) :1)
    iv[i] <- which.max(params[[1]][,iv[i+1]]* xi[i,])
  iv
}
