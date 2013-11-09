#'Compute negative log likelihood of HMM using all parameters
#'
#'This function computes the negative log likelihood of the hidden Markov model. 
#'using all parameters, untransformed.  It is used to get the covariance matrix
#'of the fitted model.
#'
#'@param parvect The vector of parameters to be estimated
#'@param obs A n x ndist matrix of data.  If ndist=1, obs must be a n x 1 matrix.
#'@param PDFs A list of PDFs for the ndist distributions.
#'@param skeleton A list with the original parameter structure used to reassemble
#'parvect
#'@param nstates Number of hidden states
#'@return The negative log likelihood of the hidden markov model.
#'@export
## function that computes the negative log-likelihood
move.HMM.mllk.full <- function(parvect,obs,PDFs,skeleton,nstates){
  n=nrow(obs)
  if(nstates==1){
    parvect=c(1,parvect)
  }
  params=relist(parvect,skeleton)
  params[[1]]=t(params[[1]])
  delta=solve(t(diag(nrow(params[[1]]))-params[[1]]+1),rep(1,nrow(params[[1]])))
  allprobs <- matrix(rep(1,nstates*n),nrow=n)#f(y_t|s_t=k)
  if(nstates>1){
    nparam=unlist(lapply(params,ncol))[-1]
    }else{
      nparam=unlist(lapply(params,ncol))
    }
  ndists=length(PDFs)
  #make index for NAs
  use=!is.na(obs)*1
  for(i in 1:ndists){
    if(nparam[i]==2){
      #for 2 parameter distributions
      for (j in 1:nstates){
        allprobs[use[,i],j] <- allprobs[use[,i],j]*PDFs[[i]](obs[use[,i],i],params[[i+1]][j,1],params[[i+1]][j,2])
      }
    }else if(nparam[i]==1){
      #for 1 parameter distributions. 
      for (j in 1:nstates){
        allprobs[use[,i],j] <- allprobs[use[,i],j]*PDFs[[i]](obs[use[,i],i],params[[i+1]][j])
      }
    }else if(nparam[i]==3){
      #for 3 parameter distributions
      for (j in 1:nstates){
        allprobs[use[,i],j] <- allprobs[use[,i],j]*PDFs[[i]](obs[use[,i],i],params[[i+1]][j,1],params[[i+1]][j,2],params[[i+1]][j,3])
      }
    }
  }
  foo <- delta 
  lscale <- 0
  gamma=params[[1]]
  for (i in 1:n){
    # foo*Gamma is Pr(s_t=k|Y_t-1) is transition matrix
    foo <- foo%*%gamma*allprobs[i,]  
    sumfoo <- sum(foo) #f_t+1,t
    lscale <- lscale+log(sumfoo) #adding log likelihood contributions
    foo <- foo/sumfoo
  }
  mllk <- -lscale
  mllk
}
