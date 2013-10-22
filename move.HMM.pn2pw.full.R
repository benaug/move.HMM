#'Transform parameters to the working scale
#'
#'This function, modified from Langrock et al. (2012), transforms parameters
#'to the working scale for maximization.
#'
#'@param params A list of length d+1 containing matrices of starting parameter
#'values.  The first element of the list must be the starting values for the
#'transition matrix.  If any distributions only have 1 parameter, the list entry
#'must be a nstates x 1 matrix.
#'@param transforms A list of transformations used to transform
#'the parameters to the working scale for maximization.
#'@param nstates Number of hidden states
#'@return A vector of transformed parameters
#'@export

## function that transforms each of the (possibly constrained) parameters to the working scale
move.HMM.pn2pw.full <- function(transforms,params,nstates){
  #transform gamma
  #params[[1]] <- qlogis(params[[1]])
  if(nstates>1){
    params[[1]] <-  as.vector((params[[1]]))
  }
  #for each other variable
  ndist=length(params)
  for(i in 2:ndist){
    #if not circular distribution
    if(transforms$type[i-1]==1){
      #loop over parameter types
      for(j in 1:(max(ncol(params[[i]]),1))){
        params[[i]][,j]=transforms[[i-1]][[j]](params[[i]][,j])
      }
      }else{#if circular dist
        #loop over direction parameters
        for(j in 1:(max(nrow(params[[i]]),1))){
        #for(j in 1:(max(ncol(params[[i]]),1))){
          #direction parameter
          params[[i]][j,1]=transforms[[i-1]][[j]](params[[i]][j,1])
        }
        #concentration parameter
        params[[i]][,2]=transforms[[i-1]][[length(transforms[[i-1]])]](params[[i]][,2])
  }
  }
  params=unlist(params)
  if(nstates==1){
    params=params[-1]
  }
  if(any(params==-Inf)){
    warning("Invalid starting values")
  }
  return(unlist(params))
}