#'Compute negative log likelihood of HSMM using all parameters
#'
#'This function computes the negative log likelihood of the hidden Markov model. 
#'using all parameters, untransformed.  It is used to get the covariance matrix
#'of the fitted model.
#'
#'@param parvect The vector of parameters to be estimated
#'@param obs A n x ndist matrix of data.  If ndist=1, obs must be a n x 1 matrix.
#'@param PDFs A list of PDFs for the dwell time and ndist observation distributions.
#'@param CDFs A list of CDFs for the dwell time and ndist observation distributions.
#'@param skeleton A list with the original parameter structure used to reassemble
#'parvect
#'@param nstates Number of hidden states
#'@param m1 a vector of length nstates that specifies how many states will be used to approximate each
#'state of the HSMM (see Langrock and Zuchinni 2011)
#'@param ini numeric value that specifies how the initial state distribution is calculated. 0 sets the
#'initial distribution to the stationary distribution.  If this matrix is not invertible, 1 sets
#'the initial distribution for each state within each state agreggate to 1/m(state).
#'@param useRcpp Logical indicating whether or not to use Rcpp.
#'@return The negative log likelihood of the hidden markov model.
#'@include gen.Gamma.R
#'@export
## function that computes the negative log-likelihood
move.HSMM.mllk.full <- function(parvect,obs,PDFs,CDFs,skeleton,nstates,m1,ini,useRcpp=FALSE){
  n=nrow(obs)
  if(nstates>2){
    #Put 0's from dwell times back in and put tpm back in correct order
    idx=seq(1,nstates^2,nstates+1)
    for(i in 1:nstates){
      parvect=append(parvect,0,idx[i]-1)
    }
    params=relist(parvect,skeleton)
    params[[1]]=t(params[[1]])
  }else{
    params=relist(parvect,skeleton)
  }
  Gamma <- gen.Gamma(m=m1,params,PDFs,CDFs)
  if(nstates>2){
    #Remove t.p.m.
    params[[1]]=NULL
  }
  nparam=unlist(lapply(params,ncol))
  sm1=sum(m1)
  delta <- solve(t(diag(sm1)-Gamma+1),rep(1,sm1))
  allprobs <- matrix(rep(1,(sm1)*n),nrow=n)
  cm1=cumsum(m1)
  mstart=c(1,cm1+1)
  mstart=mstart[-length(mstart)]
  mstop=cm1
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
  foo <- delta 
  if(class(useRcpp)=="CFunc"){
    foo=matrix(foo,ncol=sm1)
    mllk=useRcpp(Gamma,allprobs,foo)
  }else{
    lscale <- 0
    for (i in 1:n){
      foo <- foo%*%Gamma*allprobs[i,]  
      sumfoo <- sum(foo) #f_t+1,t
      lscale <- lscale+log(sumfoo) #adding log likelihood contributions
      foo <- foo/sumfoo
    }
    mllk <- -lscale
  }
  mllk
}
