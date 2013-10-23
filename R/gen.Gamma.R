#'Generate t.p.m. to approximate HSMM with HMM
#'
#'This function generates a t.p.m. of a HMM that approximates a HSMM (see Langrock and Zucchini, 2011) 
#'@param m vector of length nstates indicating the number of states to be in each state aggregate (see Langrock and Zuchinni 2011).
#'@param params A list of length ndist+1 containing matrices of starting parameter
#'values.  If nstates=2, the first element of the list must be the starting values for the
#'dwell time distribution.   If nstates>2, the first element must be the starting values for 
#'t.p.m., which must have 0's on the diagonal and rowSums=1.
#'@param PDFs A list of PDFs for the ndist distributions.
#'@param CDFs A list of CDFs for the ndist distributions.
#'@return A vector of shifted negative binomial pdf values
#'@export
#'

gen.Gamma <- function(m,params,PDFs,CDFs){
  Gamma <- diag(sum(m))*0
  nstates=nrow(params[[1]])
  #Locate dwell time distribution parameters and off-diag transition probs if nstates>2
  if(nstates>2){
    dwell=params[[2]]
    gam=params[[1]]
  }else{
    dwell=params[[1]]
    gam=matrix(rep(1,4),nrow=2)
  }
  nparam=max(1,ncol(dwell))
  within=vector("list",nstates)
  between=vector("list",nstates)
  
  #Build within and between state transition probabilities
  for(j in 1:nstates){
    ## Build within
    within[[j]]=matrix(0,nrow=m[j],ncol=m[j])
    if(nparam==2){
      p <- PDFs[[1]](1:m[j],dwell[j,1],dwell[j,2])
      dd=c(1,1-CDFs[[1]](1:(m[j]-1),dwell[j,1],dwell[j,2]))
    }else if(nparam==1){
      p <- PDFs[[1]](1:m[j],dwell[j,1])
      dd=c(1,1-CDFs[[1]](1:(m[j]-1),dwell[j,1]))
    }else if(nparam==3){
      p <- PDFs[[1]](1:m[j],dwell[j,1],dwell[j,2],dwell[j,3])
      dd=c(1,1-CDFs[[1]](1:(m[j]-1),dwell[j,1],dwell[j,2],dwell[j,3]))
    }
    haz=p/dd
    test=which(dd<1e-12)
    if(length(test)>0){
      haz[test]=1
    }
    diag(within[[j]][1:(nrow(within[[j]])-1),2:ncol(within[[j]])])=1-haz[1:(length(haz)-1)]
    within[[j]][nrow(within[[j]]),ncol(within[[j]])]=1-haz[length(haz)]
    #Build between
    between[[j]]=matrix(0,nrow=m[j],ncol=m[j])
    between[[j]][,1]=haz
  }
  
  #Assemble and multiply by offdiagnoal tpm elements if nstate>2
  mstart=c(1,cumsum(m)+1)
  mstart=mstart[-length(mstart)]
  mstop=cumsum(m)
  idx=seq(1,nstates,1)
  for(j in 1:nstates){
    Gamma[mstart[j]:mstop[j],mstart[j]:mstop[j]]=within[[j]]
    idx2=idx[-j]
    for(i in 1:length(idx2)){
      Gamma[mstart[j]:mstop[j],mstart[idx2[i]]:mstop[idx2[i]]]=between[[j]]*gam[j,idx2[i]]
    }
  }
  return(Gamma)
}