#'Fit a Hidden Markov Model (HMM) via maximum likelihood
#'
#'move.HMM.mle is used to fit HMMs, allowing for multiple observation 
#'variables with different distributions (ndist=number of distributions).
#'Maximization is performed in nlm.
#'
#'@param obs A n x ndist matrix of data.  If ndist=1, obs must be a n x 1 matrix. It
#'is recommended that movement path distances are modeled at the kilometer scale
#'rather than at the meter scale.  Turning angle observations must be in
#'the interval [-pi,pi].
#'@param dists A length d vector of distributions from the following list:
#'weibull, gamma, exponential, normal, lognormal, lnorm3, posnorm,
#'invgamma, rayleigh, f, ncf, dagum, frechet, beta, binom, poisson, nbinom,
#'zapois, wrpcauchy, wrpnorm
#'@param params A list of length ndist+1 containing matrices of starting parameter
#'values.  The first element of the list must be the starting values for the
#'transition matrix.  If modeling a single behavioral state, the transition matrix
#'must be a 1 x 1 matrix with value 1 (see example).  If any distributions only
#'have 1 parameter, the list entry must be a nstates x 1 matrix.  Users should use
#'reasonable starting values.  One method of finding good starting values is to 
#'plot randomly generated data from distributions with known parameter values
#'and compare them to histograms of the data.  More complex models
#'generally require better starting values.  One strategy is to fit a one distribution
#'model, then fit a two distribution model using the MLEs for the first distribution
#'as starting values.  Signs that the model has converged to a local maximum are
#'poor residual plots, poor plot.HMM plots, transition probabilities very close to 1
#'(especially for the shifted lognormal for which only a local maximum exists),
#'and otherwise unreasonable parameter estimates.  Analysts should try multiple combinations of starting
#'values to increase confidence that a global maximum has been achieved.
#'@param stepm a positive scalar which gives the maximum allowable scaled step
#'length. stepm is used to prevent steps which would cause the optimization
#'function to overflow, to prevent the algorithm from leaving the area of
#'interest in parameter space, or to detect divergence in the algorithm.
#'stepm would be chosen small enough to prevent the first two of these
#'occurrences, but should be larger than any anticipated reasonable step.
#'@param CI Logical indicating if confidence intervals should be produced.  CIs not yet implemented.
#'@param iterlim a positive integer specifying the maximum number of iterations to be performed before the nlm is terminated.
#'@param turn Parameters determining the transformation for circular distributions.
#'turn=1 leads to support on (0,2pi) and turn=2 leads to support on (-pi,pi).  For
#'animal movement models, the "encamped" state should use turn=1 and the "traveling"
#'state should use turn=2.
#'@return A list containing model parameters, the stationary distribution, and
#'the AICc
#'@include Distributions.R
#'@include move.HMM.pw2pn.R
#'@include move.HMM.mllk.R
#'@examples \dontrun{
#'#2 states, 2 dist-lognorm, wrapped normal
#'lmean=c(-3,-1) #meanlog parameters
#'sd=c(1,1) #sdlog parameters
#'rho<-c(0.2,0.3) # wrapped normal concentration parameters
#'mu<-c(pi,0) # wrapped normal mean parameters
#'gamma0=matrix(c(0.6,0.4,0.2,0.8),byrow=T,nrow=2)
#'
#'dists=c("lognormal","wrpnorm")
#'turn=c(1,2)
#'params=vector("list",3)
#'params[[1]]=gamma0
#'params[[2]]=cbind(lmean,sd)
#'params[[3]]=cbind(mu,rho)
#'obs=move.HMM.simulate(dists,params,1000)$obs
#'move.HMM=move.HMM.mle(obs,dists,params,stepm=35,iterlim=100,turn=turn)
#'#Assess fit
#'xlim=matrix(c(0.001,-pi,2,pi),ncol=2)
#'breaks=c(200,20)
#'HMM.plot(move.HMM,xlim,breaks=breaks)
#'move.HMM.psresid(move.HMM)
#'move.HMM.Altman(move.HMM)
#'move.HMM.dwell.plot(move.HMM)
#'move.HMM.ACF(move.HMM)
#'#Get CIs
#'params=move.HMM$params
#'move.HMM=move.HMM.mle(obs,dists,params,stepm=35,iterlim=150,turn=turn,CI=T)
#'
#'#2 states, 1 dist-shifted lognormal
#'mlog=c(-4.2,-2.2) #meanlog parameters
#'sdlog=c(1,1) #sdlog parameters
#'shift=c(0.0000001,0.03)
#'gamma0=matrix(c(0.6,0.4,0.2,0.8),byrow=T,nrow=2)
#'dists=c("lnorm3")
#'params=vector("list",2)
#'params[[1]]=gamma0
#'params[[2]]=cbind(mlog,sdlog,shift)
#'obs=move.HMM.simulate(dists,params,1500)$obs
#'move.HMM=move.HMM.mle(obs,dists,params,stepm=35,iterlim=100)
#'#Assess fit
#'xlim=matrix(c(0.001,0.8),ncol=2)
#'breaks=c(200)
#'HMM.plot(move.HMM,xlim,breaks=breaks)
#'move.HMM.psresid(move.HMM)
#'move.HMM.Altman(move.HMM)
#'move.HMM.dwell.plot(move.HMM)
#'move.HMM.ACF(move.HMM)
#'#Get CIs
#'params=move.HMM$params
#'move.HMM=move.HMM.mle(obs,dists,params,stepm=35,iterlim=150,CI=T)
#'
#'#1 state, 1 dist-gamma
#'shape=c(1) #shape parameters
#'rate=c(10) #rate parameters
#'gamma0=matrix(1,byrow=T,nrow=1)
#'dists=c("gamma")
#'params=vector("list",2)
#'params[[1]]=gamma0
#'params[[2]]=cbind(shape,rate)
#'obs=move.HMM.simulate(dists,params,1000)$obs
#'move.HMM=move.HMM.mle(obs,dists,params,stepm=35,iterlim=100)
#'#Assess fit
#'xlim=matrix(c(0.001,1),ncol=2)
#'breaks=c(20)
#'HMM.plot(move.HMM,xlim,breaks=breaks)
#'move.HMM.psresid(move.HMM)
#'move.HMM.Altman(move.HMM)
#'move.HMM.dwell.plot(move.HMM)
#'move.HMM.ACF(move.HMM)
#'#Get CIs
#'params=move.HMM$params
#'move.HMM=move.HMM.mle(obs,dists,params,stepm=35,iterlim=150,CI=T)
#'
#'#1 state, 2 dist- Weibull, Wrapped cauchy
#'wshape=c(0.9) #Weibull shape parameter
#'wscale=c(0.3) #Weibull scale parameter
#'rho<-c(0.2) # wrapped cauchy concentration parameter
#'mu<-c(pi) # wrapped cauchy mean parameter
#'gamma0=matrix(1,byrow=T,nrow=1)
#'dists=c("weibull", "wrpcauchy")
#'turn=1
#'params=vector("list",3)
#'params[[1]]=gamma0
#'params[[2]]=cbind(wshape,wscale)
#'params[[3]]=cbind(mu, rho)
#'obs=move.HMM.simulate(dists,params,1000)$obs
#'move.HMM=move.HMM.mle(obs,dists,params,stepm=35,iterlim=100, turn=turn)
#'xlim=matrix(c(0.001,-pi,0,2,pi,40),ncol=2)
#'by=c(0.001,0.001,1)
#'breaks=c(200,20,20)
#'HMM.plot(move.HMM,xlim,breaks=breaks)
#'move.HMM.psresid(move.HMM)
#'move.HMM.Altman(move.HMM)
#'move.HMM.dwell.plot(move.HMM)
#'move.HMM.ACF(move.HMM)
#'#Get CIs
#'params=move.HMM$params
#'move.HMM=move.HMM.mle(obs,dists,params,stepm=35,iterlim=150,turn=turn,CI=T)
#'
#'#3 states, 1 dist lognormal
#'lmean=c(-4,-2,-0.5) #shape parameters
#'sd=c(1,1,1) #rate parameters
#'gamma0=matrix(c(0.3,0.3,0.4,0.3,0.3,0.4,0.3,0.3,0.4),byrow=T,nrow=3)
#'dists=c("lognormal")
#'params=vector("list",2)
#'params[[1]]=gamma0
#'params[[2]]=cbind(lmean,sd)
#'obs=move.HMM.simulate(dists,params,1500)$obs
#'move.HMM=move.HMM.mle(obs,dists,params,stepm=35,iterlim=200)
#'#Assess fit
#'xlim=matrix(c(0.001,1),ncol=2)
#'breaks=c(200)
#'HMM.plot(move.HMM,xlim,breaks=breaks)
#'move.HMM.psresid(move.HMM)
#'move.HMM.Altman(move.HMM)
#'move.HMM.dwell.plot(move.HMM)
#'move.HMM.ACF(move.HMM)
#'#Get CIs
#'params=move.HMM$params
#'move.HMM=move.HMM.mle(obs,dists,params,stepm=35,iterlim=150,CI=T)
#'
#'#2 states, 3 dist-lognorm, wrapped cauchy, poisson
#'#For example, this could be movement path lengths, turning angles,
#'#and accelerometer counts from a GPS-collared animal.
#'lmean=c(-3,-1) #meanlog parameters
#'sd=c(1,1) #sdlog parameters
#'rho<-c(0.2,0.3) # wrapped Cauchy concentration parameters
#'mu<-c(pi,0) # wrapped Cauchy mean parameters
#'lambda=c(2,20)
#'gamma0=matrix(c(0.6,0.4,0.2,0.8),byrow=T,nrow=2)
#'dists=c("lognormal","wrpcauchy","poisson")
#'turn=c(1,2)
#'params=vector("list",4)
#'params[[1]]=gamma0
#'params[[2]]=cbind(lmean,sd)
#'params[[3]]=cbind(mu,rho)
#'params[[4]]=matrix(lambda,ncol=1)
#'obs=move.HMM.simulate(dists,params,1500)$obs
#'move.HMM=move.HMM.mle(obs,dists,params,stepm=35,iterlim=150,turn=turn)
#'#Assess fit - note, not great--need more data.
#'xlim=matrix(c(0.001,-pi,0,2,pi,40),ncol=2)
#'by=c(0.001,0.001,1)
#'breaks=c(200,20,20)
#'HMM.plot(move.HMM,xlim,breaks=breaks)
#'move.HMM.psresid(move.HMM)
#'move.HMM.Altman(move.HMM)
#'move.HMM.dwell.plot(move.HMM)
#'move.HMM.ACF(move.HMM)
#'#Get CIs
#'params=move.HMM$params
#'move.HMM=move.HMM.mle(obs,dists,params,stepm=35,iterlim=150,turn=turn,CI=T)
#'}
#'@export
move.HMM.mle <- function(obs,dists,params,stepm=35,CI=F,iterlim=150,turn=NULL){
  #check input
  if(is.matrix(obs)==F&is.data.frame(obs)==F)stop("argument 'obs' must be a ndist x n matrix or data frame")
  if(!is.null(dim(params[[1]]))){
    if(!all(unlist(lapply(params,is.matrix))))stop("argument 'params' must contain nstate x nparam matrices")
    if(!all(nrow(params[[1]])-unlist(lapply(params,nrow))==0))stop("All parameter matrices must have the same number of rows")
  }
  if(!all(nrow(params[[1]])-unlist(lapply(params,nrow))==0))stop("All parameter matrices must have the same number of rows")
  nstates=nrow(params[[1]])
  if(any(is.element(dists,c("wrpnorm","wrpcauchy")))){
    if(is.null(turn))stop("Must input turn")
    if(length(turn)!=nstates)stop("Number of turn elements must = number of hidden states")
  }
  if(!(any(is.element(dists,c("wrpnorm","wrpcauchy"))))&(!is.null(turn)))stop("No turn argument needed--no circular distribution.")
  #Get appropriate linearizing transformations and PDFs 
  out=Distributions(dists,nstates,turn)
  if(!all(unlist(lapply(params,ncol))[2:length(params)]==out[[7]]))stop("Incorrect number of parameters supplied for at least 1 distribution.")
  transforms=out[[1]]
  inv.transforms=out[[2]]
  PDFs=out[[3]]
  skeleton=params
  parvect <- move.HMM.pn2pw(transforms,params,nstates)  
  mod <- nlm(move.HMM.mllk,parvect,obs,print.level=2,stepmax=stepm,PDFs=PDFs,skeleton=skeleton,inv.transforms=inv.transforms,nstates=nstates,iterlim=iterlim)
  mllk <- -mod$minimum
  pn <- move.HMM.pw2pn(inv.transforms,mod$estimate,skeleton,nstates)
  #t.p.m must be matrix
  if(nstates==1){
    pn$params$tmat=as.matrix(pn$params$tmat)
  }
  params=pn$params
  npar=length(parvect)
  AIC=2*npar-2*mllk
  AICc=AIC+(2*npar*(npar+1))/(nrow(obs)-npar-1)
  
  #Get CIs
  if(CI==T){
    #Get SEs from hessian
    cat("Calculating CIs")
    #should add code checking for singularity of hessian
    H <- hessian(move.HMM.mllk,mod$estimate,obs=obs,PDFs=PDFs,skeleton=skeleton,inv.transforms=inv.transforms,nstates=nstates)
    
    ##Code attempting to get CIs for transition probabilities - it would work if I had the correct hessian
    
    # parvect.full=move.HMM.pn2pw.full(transforms,pn$params,nstates)
    #H <- hessian(move.HMM.mllk.full,parvect.full,obs=obs,PDFs=PDFs,skeleton=skeleton,inv.transforms=inv.transforms,nstates=nstates)
    #build constraint matrix
    #K=matrix(0,ncol=length(parvect.full),nrow=nstates)
    #       st=1
    #       for(i in 1:nstates){
    #         K[i,st:(st+nstates-1)]=1
    #         st=st+nstates
    #       }
    #       D=H+t(K)%*%K
    #       Dinv=solve(D)
    #       KDinv=K%*%Dinv
    #       C=Dinv-Dinv%*%t(K)%*%solve(KDinv%*%t(K))%*%KDinv
    
    vars=diag(solve(H))
    se=rep(NA,length(vars))
    se[vars>0]=sqrt(vars[vars>0])
    #calculate CIs on transformed scale and back transform
    upper=mod$estimate+1.96*se
    lower=mod$estimate-1.96*se
    upper=move.HMM.pw2pn(inv.transforms,upper,skeleton,nstates)
    lower=move.HMM.pw2pn(inv.transforms,lower,skeleton,nstates)
    #Remove CIs for delta and t.p.m. - not correct
    upper$delta=rep(NA,nstates)
    lower$delta=rep(NA,nstates)
    upper$params$tmat=matrix(NA,nrow=nstates,ncol=nstates)
    lower$params$tmat=matrix(NA,nrow=nstates,ncol=nstates)
    #transpose t.p.m. for presentation of results
    pn$params$tmat=t(pn$params$tmat)
    upper$params$tmat=t(upper$params$tmat)
    lower$params$tmat=t(lower$params$tmat)
  }else{
    if(nstates==1){
      upper=lower=rep(NA,length(mod$estimate)+2)
      
    }else{
      upper=lower=rep(NA,length(unlist(pn)))
    }
  }
  
  #Make est. ci structure
  parout=cbind(unlist(pn),unlist(lower),unlist(upper))
  colnames(parout)=c("est.","95% lower","95% upper")
  par=1
  if(nstates==1){
    parout=parout[-c(1,nrow(parout)),]
  }
  if(nstates>1){
    for(i in 1:nrow(pn$params[[1]])){
      for(j in 1:nrow(pn$params[[1]])){
        if(i==j){
          rownames(parout)[par]=paste("P(",i,"|",j,")*")
          parout[par,2:3]=parout[par,3:2]
        }else{
          rownames(parout)[par]=paste("P(",i,"|",j,")")
        }
        par=par+1
      }
    }
  }
  for(k in 2:length(pn$params)){
    for(j in 1:ncol(pn$params[[k]])){
      for(i in 1:nrow(pn$params[[k]])){
        rownames(parout)[par]=paste(dists[k-1],colnames(pn$params[[k]])[j],i)
        if(CI==T){
          if(!is.na(parout[par,2])){
            if(parout[par,2]>parout[par,3]){
              parout[par,2:3]=parout[par,3:2]
            }
          }
        }
        par=par+1
      }
    }
  }
  if(nstates>1){
    for(i in 1:length(pn$delta)){
      rownames(parout)[par]=paste("Stationary",i,"*")
      par=par+1
    }
  }
  out=list(dists=dists,nstates=nstates,params=params,delta=pn$delta,parout=parout,mllk=mllk,npar=npar,AICc=AICc,turn=turn,obs=obs)
  class(out)="move.HMM"
  out
}