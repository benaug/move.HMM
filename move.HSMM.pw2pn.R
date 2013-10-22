#'Transform parameters back to real scale
#'
#'This function, modified from Langrock et al. (2012), transforms parameters
#'back to the real scale.
#'
#'@param parvect The vector of parameters to be estimated
#'@param inv.transforms A list of inverse transformations used to transform
#'parvect back to the original scale
#'@param nstates Number of hidden states
#'@param skeleton A list with the original parameter structure used to reassemble
#'parvect
#'@return A list of parameter values on the real scale with an element
#'for each distribution
#'@export
#'
move.HSMM.pw2pn <- function(inv.transforms,parvect,skeleton,nstates){
  start=1
  params=relist(parvect,skeleton)
  #back-transform t.p.m. if present
   if(nstates>2){    
     design=matrix(rep(1,nstates*nstates),nrow=nstates)
     design[1:(nstates-1),nstates]=0
     design[nstates,nstates-1]=0
     design2=design
     diag(design2)=rep(0,nstates)
     gamma=matrix(rep(0,nstates*nstates),nrow=nstates)
     gamma [design2==1] <- exp(parvect[1:((nstates-2)*nstates)])
     gamma[design!=1]=1
     gamma <- gamma/apply(gamma,1,sum)     
     start=2
     inv.transforms=c(1,inv.transforms)
     inv.transforms$type=c(1,inv.transforms$type)
     params[[1]]=gamma
   }
  skeleton2=skeleton
  skeleton2[[1]]=NULL
  if(nstates>2){
    params=c(list(gamma),relist(parvect[((nstates-2)*nstates+1):length(parvect)],skeleton2))
  }else{
    params=relist(parvect,skeleton) 
  }
  #transform all other parameters
  ndist=length(params)
  for(i in start:ndist){
    #if not circular distribution
    if(inv.transforms$type[i]==1){
      #loop over parameter types
      for(j in 1:(max(ncol(params[[i]]),1))){
        params[[i]][,j]=inv.transforms[[i]][[j]](params[[i]][,j])
      }
    }else{#if circular dist
      #loop over direction parameters
      for(j in 1:(nrow(params[[i]]))){
        #direction parameter
        params[[i]][j,1]=inv.transforms[[i]][[j]](params[[i]][j,1])
      }
      #concentration parameter
      params[[i]][,2]=inv.transforms[[i]][[length(inv.transforms[[i]])]](params[[i]][,2])
    }
  }
  names(params)=names(inv.transforms[-length(inv.transforms)])
  return(list(params=params))
}