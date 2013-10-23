#'Compute negative log likelihood of HSMM
#'
#'This function, modified from Langrock et al. (2012), computes the negative
#'log likelihood of the hidden Markov model. 
#'
#'@param parvect The vector of parameters to be estimated
#'@param obs A n x ndist matrix of data.  If ndist=1, obs must be a n x 1 matrix.
#'@param PDFs A list of PDFs for the dwell time and ndist observation distributions.
#'@param CDFs A list of CDFs for the dwell time and ndist observation distributions.
#'@param skeleton A list with the original parameter structure used to reassemble
#'parvect
#'@param inv.transforms A list of inverse transformations used to transform
#'parvect back to the original scale
#'@param nstates Number of hidden states
#'@param m a vector of length nstates that specifies how many states will be used to approximate each
#'state of the HSMM (see Langrock and Zuchinni 2011)
#'@param ini numeric value that specifies how the initial state distribution is calculated. 0 sets the
#'initial distribution to the stationary distribution.  If this matrix is not invertible, 1 sets
#'the initial distribution for each state within each state agreggate to 1/m(state).
#'@return The negative log likelihood of the hidden markov model.
#'@include gen.Gamma
#'@export
## function that computes the negative log-likelihood
move.HSMM.mllk <- function(parvect,obs,PDFs,CDFs,skeleton,inv.transforms,nstates,m1,ini){
  n <- dim(obs)[1]
  lpn <- move.HSMM.pw2pn(inv.transforms,parvect,skeleton,nstates)
  params=lpn$params
  Gamma <- gen.Gamma(m=m1,params,PDFs,CDFs)
  if(nstates>2){
    #Remove t.p.m.
    params[[1]]=NULL
  }
  if (ini==1) {delta <- rep(1/(sum(m1)),sum(m1))} # if invertibility problems
  if (ini==0) {delta <- solve(t(diag(sum(m1))-Gamma+1),rep(1,sum(m1)))}
  #if(any(delta<0))stop("Set stationary to no or both")
  #allprobs <- matrix(rep(1,nstates*n),nrow=n)#f(y_t|s_t=k)
  allprobs <- matrix(rep(1,(sum(m1))*n),nrow=n)
  mstart=c(1,cumsum(m1)+1)
  mstart=mstart[-length(mstart)]
  mstop=cumsum(m1)
  for (k in 1:n){
    for (j in 1:nstates){
      for(i in 2:length(PDFs)){
        nparam=max(1,ncol(params[[i]]))
        if(nparam==2){
          #for 2 parameter distributions
          allprobs[k,mstart[j]:mstop[j]]=allprobs[k,mstart[j]:mstop[j]]*rep(ifelse(is.na(obs[k,i-1]),1,PDFs[[i]](obs[k,i-1],params[[i]][j,1],params[[i]][j,2])),m1[j])
        }else if(nparam==1){
          #for 1 parameter distributions. 
          allprobs[k,mstart[j]:mstop[j]]=allprobs[k,mstart[j]:mstop[j]]*rep(ifelse(is.na(obs[k,i-1]),1,PDFs[[i]](obs[k,i-1],params[[i]][j])),m1[j])
        }else if(nparam==3){
          #for 3 parameter distributions
          allprobs[k,mstart[j]:mstop[j]]=allprobs[k,mstart[j]:mstop[j]]*rep(ifelse(is.na(obs[k,i-1]),1,PDFs[[i]](obs[k,i-1],params[[i]][j,1],params[[i]][j,2],params[[i]][j,3])),m1[j])
        }
      } #i index
    } # j index
  } # k index
  foo <- delta 
  lscale <- 0
  for (i in 1:n){
    # foo*Gamma is Pr(s_t=k|Y_t-1) is transition matrix
    foo <- foo%*%Gamma*allprobs[i,]  
    sumfoo <- sum(foo) #f_t+1,t
    lscale <- lscale+log(sumfoo) #adding log likelihood contributions
    foo <- foo/sumfoo
  }
  mllk <- -lscale
  mllk
}