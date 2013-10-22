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
## inverse transformation back to the natural parameter space
move.HMM.pw2pn.full <- function(inv.transforms,parvect,skeleton,nstates){
  
  #transform gamma
  if(nstates>1){
    gamma [! gamma] <- (parvect[1:nstates^2])
    delta <- solve(t(diag(nrow(gamma))-gamma+1),rep(1,nrow(gamma)))
  }else{
    delta=1
  }
  skeleton2=skeleton
  skeleton2[[1]]=NULL
  if(nstates>1){
    params=c(list(gamma),relist(parvect[(nstates^2+1):length(parvect)],skeleton2))
  }else{
    params=c(list(1),relist(parvect[(nstates^2+1):length(parvect)],skeleton2))
  }
  #for each other variable
  ndist=length(params)
  for(i in 2:ndist){
    #if not circular distribution
    if(inv.transforms$type[i-1]==1){
      #loop over parameter types
      for(j in 1:(max(ncol(params[[i]]),1))){
        params[[i]][,j]=inv.transforms[[i-1]][[j]](params[[i]][,j])
      }
    }else{#if circular dist
      #loop over direction parameters
      for(j in 1:(nrow(params[[i]]))){
        #direction parameter
        params[[i]][j,1]=inv.transforms[[i-1]][[j]](params[[i]][j,1])
      }
      #concentration parameter
      params[[i]][,2]=inv.transforms[[i-1]][[length(inv.transforms[[i-1]])]](params[[i]][,2])
    }
  }
  
  names(params)=c("tmat",names(inv.transforms)[1:(length(names(inv.transforms))-1)])
  return(list(params=params,delta=delta))
}