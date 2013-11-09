#'Add confidence intervals to a move.HMM object
#'
#'move.HMM.CI is used to add confidence intervals to a move.HMM object.  Current
#'options are parametric bootstrap percentile CIs or CIs calculated from the finite differences Hessian.
#'Bootstrapping utilizes parallel processing on a local host.  For each bootstrap sample, NAs in the original
#'data are inserted in the simulated data in the same positions.  CIs for the stationary distribution are
#'obtained using a numerical derivative and the delta method following Patterson et al. (2009).  These intervals
#'appear to be substantially too narrow compared to bootstrap CIs and should be treated with caution.  The
#'CI's based on the finite differences Hessian should generally be treated skeptically unless you've
#'done a coverage analysis for your model with the amount of data you have.  Bootstrap samples are stored in move.HMM$store boot so they can be inspected or combined with
#'more bootstrap samples.
#'
#'@param move.HMM A fitted move.HMM object.
#'@param CI A character determining which type of CI is to be calculated.  Current options are
#'"FD" for the finitie differences Hessian and "boot" for parametric bootstrapping and percentile CIs.
#'@param alpha Type I error rate for CIs.  alpha=0.05 for 95 percent CIs
#'@param B Number of bootstrap resamples
#'@param cores Number of cores to be used in parallell bootstrapping
#'@param stepm a positive scalar which gives the maximum allowable scaled step
#'length. stepm is used to prevent steps which would cause the optimization
#'function to overflow, to prevent the algorithm from leaving the area of
#'interest in parameter space, or to detect divergence in the algorithm.
#'stepm would be chosen small enough to prevent the first two of these
#'occurrences, but should be larger than any anticipated reasonable step.
#'@param iterlim a positive integer specifying the maximum number of iterations to be performed before the nlm is terminated.
#'@return A list containing the lower and upper confidence bounds
#'@include Distributions.R
#'@include move.HMM.pw2pn.R
#'@include move.HMM.mllk.R
#'@include stationaryDist.R
#'@include move.HMM.mllk.full.R
#'@export
#'
move.HMM.CI=function(move.HMM,CI="boot",alpha=0.05,B=100,cores=2,stepm,iterlim){
  nstates=move.HMM$nstates
  dists=move.HMM$dists
  params=move.HMM$params
  obs=move.HMM$obs
  n=nrow(obs)
  turn=move.HMM$turn
  if(nstates>1){
    nparams=length(unlist(params))+nstates
  }else{
    nparams=length(unlist(params))+nstates-2
  }
  if(CI=="boot"){
    cat('\n Calculating bootstrap CIs')
    require(foreach)
    require(snow)
    require(doSNOW)
    #Find locations of missing data
    idx.na=which(is.na(obs))
    any.missing=any(idx.na==TRUE)
    cl.tmp = makeCluster(rep("localhost",cores), type="SOCK")
    registerDoSNOW(cl.tmp)
    out=foreach(k=1:B, .packages=c("move.HMM","pscl","VGAM","psych","CircStats")) %dopar% {
      obs=move.HMM.simulate(dists,params,n)$obs
      #Add in missing values to data
      if(any.missing){
        obs[idx.na]=NA
      }
      move.HMM=move.HMM.mle(obs,dists,params,stepm=stepm,iterlim=iterlim,turn=turn,CI=F)
      move.HMM$parout[,1]
    }
    stopCluster(cl.tmp)
    store=matrix(unlist(out),ncol=nparams,byrow=T)
    storeboot=matrix(unlist(out),ncol=nparams,byrow=T)
    colnames(storeboot)=names(out[[1]])
    CIs=t(apply(store,2,quantile,probs=c(alpha/2,1-alpha/2)))
    lower=CIs[,1]
    upper=CIs[,2]
    move.HMM$storeboot=storeboot
  }
  if(CI=="FD"){
    cat('\n Calculating CIs using finite differences Hessian')    
    out=Distributions(dists,nstates,turn)
    transforms=out[[1]]
    inv.transforms=out[[2]]
    PDFs=out[[3]]
    skeleton=params
    delta=move.HMM$delta
    if(nstates==1){
      parvect=move.HMM.pn2pw(transforms,params,nstates) 
      H <- hessian(move.HMM.mllk,parvect,obs=obs,PDFs=PDFs,skeleton=skeleton,nstates=nstates,inv.transforms=inv.transforms)
      vars=diag(solve(H))
      se=sqrt(vars)
      est=parvect
      lower=est-qnorm(alpha/2,0,1)*se
      upper=est+qnorm(alpha/2,0,1)*se
      est=move.HMM.pw2pn(inv.transforms,parvect,skeleton,nstates)
      lower=move.HMM.pw2pn(inv.transforms,lower,skeleton,nstates)
      upper=move.HMM.pw2pn(inv.transforms,upper,skeleton,nstates)
      est=unlist(est)[2:(nparams+1)]
      lower=unlist(lower)[2:(nparams+1)]
      upper=unlist(upper)[2:(nparams+1)]
    }else{
      #transform tpm so that it unlists in right order
      params$tmat=t(params$tmat)
      parvect=unlist(params)
      H <- hessian(move.HMM.mllk.full,parvect,obs=obs,PDFs=PDFs,skeleton=skeleton,nstates=nstates)  
      #build constraint matrix
      K=matrix(0,ncol=length(parvect),nrow=nstates)
      st=1
      for(i in 1:nstates){
        K[i,st:(st+nstates-1)]=1
        st=st+nstates
      }
      D=H+t(K)%*%K
      Dinv=solve(D)
      KDinv=K%*%Dinv
      C=Dinv-Dinv%*%t(K)%*%solve(KDinv%*%t(K))%*%KDinv
      vars=diag(C)
      se=sqrt(vars)
      est=unlist(params)
      upper=est+qnorm(alpha/2,0,1)*se
      lower=est-qnorm(alpha/2,0,1)*se
      #Get stationary derivative
      gamvect=unlist(t(params$tmat))[1:nstates^2]
      G <- jacobian(stationaryDist,gamvect)
      vardelta=G%*%C[1:nstates^2,1:nstates^2]%*%t(G)
      sedelta=sqrt(diag(vardelta))
      deltlow=delta-qnorm(alpha/2,0,1)*sedelta
      deltup=delta+qnorm(alpha/2,0,1)*sedelta
      #Add CIs for delta
      upper=c(upper,deltup)
      lower=c(lower,deltlow)
    }
    #Flip any uppers<lowers
    flip=upper<lower
    both=cbind(lower,upper)
    if(any(flip)){
      for(i in 1:length(flip)){
        if(flip[i]){
          both[i,1:2]=both[i,2:1]
        }
      }
    }
    lower=both[,1]
    upper=both[,2]
  }
#    if(nstates==1){
#      lower=c(1,lower,1)
#      upper=c(1,upper,1)
#    }
  move.HMM$parout[,2:3]=cbind(lower,upper)
  return(move.HMM)
}