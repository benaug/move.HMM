#'Assign states using the Viterbi algorithm
#'
#'This function, modified from Zucchini and MacDonald (2009), assigns states
#'to observations using the Viterbi algorithm.  It takes as input a 
#'move.HSMM object and an optional vector containing the starting state 
#'probabilities.
#'
#'@param move.HSMM A move.HSMM object containing a fitted HSMM model.
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
#'obs=move.HSMM.simulate(dists,params,1000)$obs
#'turn=c(1,2)
#'move.HSMM=move.HSMM.mle(obs,dists,params,stepm=35,iterlim=100,turn=turn)
#'#get Viterbi state assignments
#'move.HSMM.viterbi(move.HSMM)
#'}
#'@export
#'
move.HSMM.viterbi <-function(move.HSMM,delta=NULL){
  obs=move.HSMM$obs
  params=move.HSMM$params
  nstates=move.HSMM$nstates
  dists=move.HSMM$dists
  m=move.HSMM$m1
  out=Distributions(dists,nstates)
  PDFs=out[[3]]
  CDFs=out[[4]]
  Gamma <- gen.Gamma(m,params,PDFs,CDFs)
  if(nstates>2){
    params[[1]]=NULL
  }
  if(is.null(delta))delta=solve(t(diag(sum(m))-Gamma+1),rep(1,sum(m)))
  n <- nrow(obs)
  allprobs <- matrix(rep(1,(sum(m))*n),nrow=n)
  mstart=c(1,cumsum(m)+1)
  mstart=mstart[-length(mstart)]
  mstop=cumsum(m)
  for (k in 1:n){
    for (j in 1:nstates){
      for(i in 2:length(PDFs)){
        nparam=max(1,ncol(params[[i]]))
        if(nparam==2){
          #for 2 parameter distributions
          allprobs[k,mstart[j]:mstop[j]]=allprobs[k,mstart[j]:mstop[j]]*rep(ifelse(is.na(obs[k,i-1]),1,PDFs[[i]](obs[k,i-1],params[[i]][j,1],params[[i]][j,2])),m[j])
        }else if(nparam==1){
          #for 1 parameter distributions. 
          allprobs[k,mstart[j]:mstop[j]]=allprobs[k,mstart[j]:mstop[j]]*rep(ifelse(is.na(obs[k,i-1]),1,PDFs[[i]](obs[k,i-1],params[[i]][j])),m[j])
        }else if(nparam==3){
          #for 3 parameter distributions
          allprobs[k,mstart[j]:mstop[j]]=allprobs[k,mstart[j]:mstop[j]]*rep(ifelse(is.na(obs[k,i-1]),1,PDFs[[i]](obs[k,i-1],params[[i]][j,1],params[[i]][j,2],params[[i]][j,3])),m[j])
        }
      } #i index
    } # j index
  } # k index
  xi <- matrix(0,n,sum(m))
  foo <- delta*allprobs [1,]
  xi[1,] <- foo/sum(foo)
  for (i in 2:n)
  {
    foo <- apply(xi[i-1,]*Gamma ,2,max)*allprobs[i,]
    xi[i,] <- foo/sum(foo)
  }
  iv <-numeric(n)
  iv[n] <-which.max(xi[n,])
  for (i in (n-1) :1)
    iv[i] <- which.max(Gamma[,iv[i+1]]* xi[i,])
  statepreds=numeric(n)
  mstart=c(1,cumsum(m)+1)
  mstart=mstart[-length(mstart)]
  mstop=cumsum(m)
  for(i in 1:nstates){
    statepreds[((iv>=mstart[i])&(iv<=mstop[i]))]=i
  }
  statepreds
}