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
move.HMM.mllk.full <- function(parvect,obs,PDFs,skeleton,nstates){
  n=nrow(obs)
  params=relist(parvect,skeleton)
  delta=solve(t(diag(nrow(params[[1]]))-params[[1]]+1),rep(1,nrow(params[[1]])))
  allprobs <- matrix(rep(1,nstates*n),nrow=n)#f(y_t|s_t=k)
  for (k in 1:n){
      for (j in 1:nstates){
        for(i in 1:length(PDFs)){
          nparam=max(1,ncol(params[[i+1]]))
          if(nparam==2){
            #for 2 parameter distributions
            allprobs[k,j]=allprobs[k,j]*ifelse(is.na(obs[k,i]),1,PDFs[[i]](obs[k,i],params[[i+1]][j,1],params[[i+1]][j,2]))
          }else if(nparam==1){
            #for 1 parameter distributions. 
            allprobs[k,j]=allprobs[k,j]*ifelse(is.na(obs[k,i]),1,PDFs[[i]](obs[k,i],params[[i+1]][j]))
          }else if(nparam==3){
            #for 3 parameter distributions
            allprobs[k,j]=allprobs[k,j]*ifelse(is.na(obs[k,i]),1,PDFs[[i]](obs[k,i],params[[i+1]][j,1],params[[i+1]][j,2],params[[i+1]][j,3]))
          }
        } #i index
      } # j index
  } # k index
  foo <- delta 
  lscale <- 0
  for (i in 1:n){
    # foo*Gamma is Pr(s_t=k|Y_t-1) is transition matrix
    foo <- foo%*%params[[1]]*allprobs[i,]  
    sumfoo <- sum(foo) #f_t+1,t
    lscale <- lscale+log(sumfoo) #adding log likelihood contributions
    foo <- foo/sumfoo
  }
  mllk <- -lscale
  mllk
}
