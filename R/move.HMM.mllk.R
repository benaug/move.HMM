#'Compute negative log likelihood of HMM
#'
#'This function, modified from Langrock et al. (2012), computes the negative
#'log likelihood of the hidden Markov model. 
#'
#'@param parvect The vector of parameters to be estimated
#'@param obs A n x ndist matrix of data.  If ndist=1, obs must be a n x 1 matrix.
#'@param PDFs A list of PDFs for the ndist distributions.
#'@param skeleton A list with the original parameter structure used to reassemble
#'parvect
#'@param inv.transforms A list of inverse transformations used to transform
#'parvect back to the original scale
#'@param nstates Number of hidden states
#'@return The negative log likelihood of the hidden markov model.
#'@export
## function that computes the negative log-likelihood
move.HMM.mllk <- function(parvect,obs,PDFs,skeleton,inv.transforms,nstates){
  n <- dim(obs)[1]
  lpn <- move.HMM.pw2pn(inv.transforms,parvect,skeleton,nstates)
  params <- lpn$params
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
  foo <- lpn$delta 
  lscale <- 0
  gamma=params[[1]]
  for (i in 1:n){
    foo <- foo%*%gamma*allprobs[i,]  
    sumfoo <- sum(foo) #f_t+1,t
    lscale <- lscale+log(sumfoo)
    foo <- foo/sumfoo
  }
  mllk <- -lscale
  mllk
}
