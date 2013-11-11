#'Add confidence intervals to a move.HSMM object
#'
#'move.HSMM.CI is used to add confidence intervals to a move.HSMM object.  Current
#'options are parametric bootstrap percentile CIs or CIs calculated from the finite differences Hessian.
#'Bootstrapping utilizes parallel processing on a local host.  For each bootstrap sample, NAs in the original
#'data are inserted in the simulated data in the same positions.  CI's based on the finite differences Hessian should generally be treated skeptically unless you've
#'done a coverage analysis for your model with the amount of data you have.  Bootstrap samples are stored in move.HSMM$store boot so they can be inspected or combined with
#'more bootstrap samples.
#'
#'@param move.HSMM A fitted move.HSMM object.
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
#'@param useRcpp Logical indicating whether or not to use Rcpp.
#'@return A list containing the lower and upper confidence bounds
#'@include Distributions.R
#'@include move.HSMM.pw2pn.R
#'@include move.HSMM.mllk.R
#'@include move.HSMM.mllk.full.R
#'@export
#'
move.HSMM.CI=function(move.HSMM,CI="boot",alpha=0.05,B=100,cores=2,stepm,iterlim,useRcpp=useRcpp){
  nstates=move.HSMM$nstates
  dists=move.HSMM$dists
  params=move.HSMM$params
  obs=move.HSMM$obs
  n=nrow(obs)
  m1=move.HSMM$m1
  turn=move.HSMM$turn
  nparams=length(unlist(params))
  
  if(CI=="boot"){
    cat('\n Calculating bootstrap CIs (This will take a while when # cores is small and/or B is large)')
    suppressMessages(require(foreach))
    suppressMessages(require(snow))
    suppressMessages(require(doSNOW))
    #add in space for stationary dist
    nparams=nparams+nstates
    #Find locations of missing data
    idx.na=which(is.na(obs))
    any.missing=any(idx.na==TRUE)
    cl.tmp = makeCluster(rep("localhost",cores), type="SOCK")
    registerDoSNOW(cl.tmp)
    if(useRcpp){
      out=foreach(k=1:B, .packages=c("move.HMM","pscl","VGAM","psych","CircStats","Rcpp","RcppArmadillo","inline")) %dopar% {
        obs=move.HSMM.simulate(dists,params,n,nstates)$obs
        #Add in missing values to data
        if(any.missing){
          obs[idx.na]=NA
        }
        move.HSMM=move.HSMM.mle(obs,dists,params,stepm=stepm,CI=F,iterlim=iterlim,turn=turn,m1=m1,useRcpp=TRUE)
        c(move.HSMM$parout[,1],move.HSMM$delta[,1])
      }
    }else{
      out=foreach(k=1:B, .packages=c("move.HMM","pscl","VGAM","psych","CircStats")) %dopar% {
        obs=move.HSMM.simulate(dists,params,n,nstates)$obs
        #Add in missing values to data
        if(any.missing){
          obs[idx.na]=NA
        }
        move.HSMM=move.HSMM.mle(obs,dists,params,stepm=stepm,CI=F,iterlim=iterlim,turn=turn,m1=m1)
        c(move.HSMM$parout[,1],move.HSMM$delta[,1])
      }
    }
    stopCluster(cl.tmp)
    storeboot=matrix(unlist(out),ncol=nparams,byrow=T)
    colnames(storeboot)=names(out[[1]])
    CIs=t(apply(storeboot,2,quantile,probs=c(alpha/2,1-alpha/2)))
    lower=CIs[1:(nparams-nstates),1]
    upper=CIs[1:(nparams-nstates),2]
    deltlow=CIs[(nparams-nstates+1):nparams,1]
    deltup=CIs[(nparams-nstates+1):nparams,2]
    move.HSMM$storeboot=storeboot
  }
  if(CI=="FD"){
    cat('\n Calculating CIs using finite differences Hessian')
    if(useRcpp){
      suppressMessages(require(Rcpp))
      suppressMessages(require(inline))
      suppressMessages(require(RcppArmadillo))
      code <- '
   arma::mat gam = Rcpp::as<arma::mat>(Gamma);
   arma::mat foo2 = Rcpp::as<arma::mat>(foo);
   arma::mat probs = Rcpp::as<arma::mat>(allprobs);
   int n = probs.n_rows;
   arma::mat lscale(1,1);
   arma::mat sumfoo(1,1);
   for (int i=0; i<n; i++) {
     foo2=foo2*gam%probs.row(i);
     sumfoo=sum(foo2,1);
     lscale=lscale+log(sumfoo);
     double sf =  arma::as_scalar(sumfoo);
     foo2=foo2/sf;

   }
   return Rcpp::wrap(-lscale);
 '
      useRcpp <- cxxfunction(signature(Gamma="numeric",allprobs="numeric",foo="numeric"),
                             code,plugin="RcppArmadillo")
    }
    out=Distributions(dists,nstates,turn)
    transforms=out[[1]]
    inv.transforms=out[[2]]
    PDFs=out[[3]]
    CDFs=out[[4]]
    skeleton=params
    delta=move.HSMM$delta
    #transform so params unlist in correct order
    if(nstates>2){
      params[[1]]=t(params[[1]])
      parvect=unlist(params)
      #Remove 0's from dwell times
      parvect=parvect[-seq(1,nstates^2,nstates+1)]
    }else{
      parvect=unlist(params)
    }
    H <- hessian(move.HSMM.mllk.full,parvect,obs=obs,PDFs=PDFs,CDFs=CDFs,skeleton=skeleton,nstates=nstates,m1=m1,ini=0,useRcpp=useRcpp)
    if(nstates>2){
      #build constraint matrix
      K=matrix(0,ncol=length(parvect),nrow=nstates)
      st=1
      for(i in 1:nstates){
        K[i,st:(st+nstates-2)]=1
        st=st+nstates-1
      }
      D=H+t(K)%*%K
      Dinv=solve(D)
      KDinv=K%*%Dinv
      C=Dinv-Dinv%*%t(K)%*%solve(KDinv%*%t(K))%*%KDinv
      vars=diag(C)
      #calculate CIs
      se=sqrt(vars)
      upper=parvect+qnorm(alpha/2,0,1)*se
      lower=parvect-qnorm(alpha/2,0,1)*se
      #Put 0's from dwell times back in
      idx=seq(1,nstates^2,nstates+1)
      for(i in 1:nstates){
        upper=append(upper,0,idx[i]-1)
        lower=append(lower,0,idx[i]-1)        
      }
    }else{
      vars=diag(solve(H))
      se=sqrt(vars)
      upper=parvect+qnorm(alpha/2,0,1)*se
      lower=parvect-qnorm(alpha/2,0,1)*se
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
    deltlow=rep(NA,nstates)
    deltup=rep(NA,nstates)
  }
  move.HSMM$parout[,2:3]=cbind(lower,upper)
  move.HSMM$delta[,2:3]=cbind(deltlow,deltup)
  return(move.HSMM)
}