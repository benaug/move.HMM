#'Compute conditional and posterior state probabilities
#'
#'This function, modified from Zucchini and MacDonald (2009), computes the
#'conditional state  probabilities and posterior state probabilities
#'(Patterson et al. 2009). It takes as input a move.HSMM object.
#'
#'@param move.HSMM A move.HSMM object containing a fitted HSMM model.
#'@include move.HSMM.lalphabeta.R
#'@return A list of conditional and posterior state probabilities.
#'#'@examples \dontrun{
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
#'obs=move.HSMM.simulate(dists,params,1000)
#'turn=c(1,2)
#'move.HSMM=move.HSMM.mle(obs,dists,params,stepm=35,iterlim=100,turn=turn)
#'#get conditional and posterior state probabilities
#'move.HSMM.state_probs(move.HSMM)
#'}
#'@export
#'
move.HSMM.state_probs <-  function(move.HSMM){
  obs=move.HSMM$obs
  nstates=move.HSMM$nstates
  params=move.HSMM$params
  turn=move.HSMM$turn
  dists=move.HSMM$dists
  m=move.HSMM$m
  out=Distributions(dists,nstates,turn)
  PDFs=out[[3]]
  CDFs=out[[4]]
  Gamma <- gen.Gamma(m,params,PDFs,CDFs)
  #forward probabilities  
  delta <- solve(t(diag(sum(m))-Gamma+1),rep(1,sum(m)))
  n <- nrow(obs)
  fb <- move.HSMM.lalphabeta(move.HSMM)
  la <- fb$la
  lb <- fb$lb
  c <- max(la[n,])
  llk <- c+log(sum(exp(la[n,]-c)))
  stateprobs <- matrix(rep(1,(sum(m))*n),nrow=n)
  for (i in 1:n){
    stateprobs[i,]<-exp(la[i,]+lb[i,]-llk)
  }

  #backwards probabilities
  backprobs <- matrix(rep(1,(sum(m))*n),nrow=n)
  backprobs[n,]=stateprobs[n,]
  for(i in (n-1):1){
    backprobs[i,]=stateprobs[i,]*(Gamma%*%t(backprobs[i+1,]/(stateprobs[i,]%*%Gamma)))
  }
  
  #Combine state aggregates
  stateprobs2=matrix(NA,nrow=n,ncol=nstates)
  backprobs2=matrix(NA,nrow=n,ncol=nstates)
  mstart=c(1,cumsum(m)+1)
  mstart=mstart[-length(mstart)]
  mstop=cumsum(m)
  for(i in 1:nstates){
    stateprobs2[,i]=rowSums(stateprobs[,mstart[i]:mstop[i]])
    backprobs2[,i]=rowSums(backprobs[,mstart[i]:mstop[i]])
  }
  return(list(forward=stateprobs2,backward=backprobs2))
}