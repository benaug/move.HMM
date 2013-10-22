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
move.HSMM.pn2pw <- function(transforms,params,nstates){
  start=1
  #transform t.p.m. if nstates>2.  Off diagonal elements are all 0.  Model last element
  #of each row, except for the last state where we model the second to last since the 
  #last is 0
  if(nstates>2){
    gamma=params[[1]]
    design=matrix(rep(1,nstates*nstates),nrow=nstates)
    design[1:(nstates-1),nstates]=0
    design[nstates,nstates-1]=0
    fixed=gamma[design==0]
    fixed=c(fixed[2:length(fixed)],fixed[1]) #reshape
    foo <- log(gamma/fixed)
    diag(design)=rep(0,nstates)
    params[[1]] <- as.vector(foo[design==1])
    start=2
    transforms=c(1,transforms)
    transforms$type=c(1,transforms$type)
  }
  ndist=length(params)
  for(i in start:ndist){
    #if not circular distribution
    if(transforms$type[i]==1){
      #loop over parameter types
      for(j in 1:(max(ncol(params[[i]]),1))){
        params[[i]][,j]=transforms[[i]][[j]](params[[i]][,j])
      }
    }else{#if circular dist
      #loop over direction parameters
      #for(j in 1:(max(ncol(params[[i]]),1))){
      for(j in 1:(max(nrow(params[[i]]),1))){ 
        #direction parameter
        params[[i]][j,1]=transforms[[i]][[j]](params[[i]][j,1])
      }
      #concentration parameter
      params[[i]][,2]=transforms[[i]][[length(transforms[[i]])]](params[[i]][,2])
    }
  }
  params=unlist(params)
  if(any(params==-Inf)){
    warning("Invalid starting values")
  }
  return(unlist(params))
}