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
  m1=move.HSMM$m1
  out=Distributions(dists,nstates)
  PDFs=out[[3]]
  CDFs=out[[4]]
  Gamma <- gen.Gamma(m1,params,PDFs,CDFs)
  if(nstates>2){
    params[[1]]=NULL
  }
  if(is.null(delta))delta=solve(t(diag(sum(m1))-Gamma+1),rep(1,sum(m1)))
  n <- nrow(obs)
  allprobs <- matrix(rep(1,(sum(m1))*n),nrow=n)
  mstart=c(1,cumsum(m1)+1)
  mstart=mstart[-length(mstart)]
  mstop=cumsum(m1)
  nparam=unlist(lapply(params,ncol))
  ndists=length(PDFs)
  #make index for NAs
  use=!is.na(obs)*1
  for(i in 2:ndists){
    if(nparam[i]==2){
      #for 2 parameter distributions
      for (j in 1:nstates){
        allprobs[use[,i-1],mstart[j]:mstop[j]] <- allprobs[use[,i-1],mstart[j]:mstop[j]]*matrix(rep(PDFs[[i]](obs[use[,i-1],i-1],params[[i]][j,1],params[[i]][j,2]),m1[j]),ncol=m1[j])
      }
    }else if(nparam[i]==1){
      #for 1 parameter distributions. 
      for (j in 1:nstates){
        allprobs[use[,i-1],mstart[j]:mstop[j]] <- allprobs[use[,i-1],mstart[j]:mstop[j]]*matrix(rep(PDFs[[i]](obs[use[,i-1],i-1],params[[i]][j]),m1[j]),ncol=m1[j])
      }
    }else if(nparam[i]==3){
      #for 3 parameter distributions
      for (j in 1:nstates){
        allprobs[use[,i-1],mstart[j]:mstop[j]] <- allprobs[use[,i-1],mstart[j]:mstop[j]]*matrix(rep(PDFs[[i]](obs[use[,i-1],i-1],params[[i]][j,1],params[[i]][j,2],params[[i]][j,3]),m1[j]),ncol=m1[j])
      }
    }
  }
  xi <- matrix(0,n,sum(m1))
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
  mstart=c(1,cumsum(m1)+1)
  mstart=mstart[-length(mstart)]
  mstop=cumsum(m1)
  for(i in 1:nstates){
    statepreds[((iv>=mstart[i])&(iv<=mstop[i]))]=i
  }
  statepreds
}