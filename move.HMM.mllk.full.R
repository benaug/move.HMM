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
move.HMM.mllk.full <- function(parvect,obs,PDFs,skeleton,inv.transforms,nstates){
  n <- dim(obs)[1]
  lpn <- move.HMM.pw2pn.full(inv.transforms,parvect,skeleton,nstates)
  params=lpn$params
  allprobs <- matrix(rep(1,nstates*n),nrow=n)#f(y_t|s_t=k)
  for (k in 1:n){
    #if (is.na(obs[k,1])) {
    #  allprobs[k,] <- rep(1,nstates)
    #}  this is redundant
    #if (!is.na(obs[k,1])) {
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
   # } # non-miss check--I don't think this works correctly if data is missing for only 1 obs variable
  } # k index
  foo <- lpn$delta 
  lscale <- 0
  #allprobs[allprobs==0]=1  #rigged for location parameter dists-not correct
  for (i in 1:n){
    # foo*Gamma is Pr(s_t=k|Y_t-1) is transition matrix
    foo <- foo%*%lpn$params[[1]]*allprobs[i,]  
    sumfoo <- sum(foo) #f_t+1,t
    lscale <- lscale+log(sumfoo) #adding log likelihood contributions
    foo <- foo/sumfoo
  }
  mllk <- -lscale
  mllk
}
