#'Fit a Hidden Semi Markov Model (HSMM) via maximum likelihood
#'
#'move.HSMM.mle is used to fit HSMMs, allowing for multiple observation 
#'variables with different distributions (ndist=number of distributions).  It approximates
#'a HSMM with a HMM as outlined in Langrock et al. (2012).
#'Maximization is performed in nlm. Currently only 2 hidden states are supported
#'but this will be extended to an arbitrany number of states >1.
#'
#'@param obs A n x ndist matrix of data.  If ndist=1, obs must be a n x 1 matrix. It
#'is recommended that movement path distances are modeled at the kilometer scale
#'rather than at the meter scale.
#'@param dists A length d vector of distributions.  The first distribution must be
#'the dwell time distribution chosen from the following list: shiftpois,shifnegbin. The
#'subsequent distributions are for observation variables and must be chosen from the following list:
#'weibull, gamma, exponential, normal, lognormal, lnorm3, posnorm,
#'invgamma, rayleigh, f, ncf, dagum, frechet, beta, binom, poisson, nbinom,
#'zapois, wrpcauchy, wrpnorm.  Note wrpnorm is much slower to evaluate than wrpcauchy.
#'Differences in the amount of time taken to maximize can be substantial.
#'@param params A list containing matrices of starting parameter
#'values.  The structure of this list differs between 2 state models and models with
#'greater than 2 states.  For 2 state models, the first element of the list must be the starting values for the
#'dwell time distribution.  For 3+ state models, the first element of the list must be the starting values
#'for the t.p.m.  Since dwell times are modeled explicitly, the diagonal of the t.p.m must be all zeros.  Further, 
#'rows should still sum to 1.  The second element of the list will now be the dwell time distribution
#'parameter starting values.  See the example below for a demonstration.  As before, if any distributions only
#'have 1 parameter, the list entry must be a nstates x 1 matrix.  Users should use
#'reasonable starting values.  One method of finding good starting values is to 
#'plot randomly-generated data from distributions with known values and compare
#' them to histograms of the data.  More complex models
#'generally require better starting values.  One strategy for fitting HSMMs is to first fit
#'a HMM and use the MLEs as starting values for the HSMM.  Signs that the model has converged
#'to a local maximum are poor residual plots, poor plot.HMM plots, transition probabilities very close to 1
#'(especially for the shifted lognormal for which only a local maximum exists),
#'and otherwise unreasonable parameter estimates.  Analysts should try multiple combinations of starting
#'values to increase confidence that a global maximum has been achieved.
#'@param stepm a positive scalar which gives the maximum allowable scaled step
#'length. stepm is used to prevent steps which would cause the optimization
#'function to overflow, to prevent the algorithm from leaving the area of
#'interest in parameter space, or to detect divergence in the algorithm.
#'stepm would be chosen small enough to prevent the first two of these
#'occurrences, but should be larger than any anticipated reasonable step.
#'@param iterlim a positive integer specifying the maximum number of iterations to be performed before the nlm is terminated.
#'@param turn Parameters determining the transformation for circular distributions.
#'turn=1 leads to support on (0,2pi) and turn=2 leads to support on (-pi,pi).  For
#'animal movement models, the "encamped" state should use turn=1 and the "traveling"
#'state should use turn=2.
#'@param CI A character determining which type of CI is to be calculated.  Current options are
#'"FD" for the finitie differences Hessian and "boot" for parametric bootstrapping and percentile CIs.
#'@param alpha Type I error rate for CIs.  alpha=0.05 for 95 percent CIs
#'@param B Number of bootstrap resamples
#'@param cores Number of cores to be used in parallell bootstrapping
#'@param m1 vector of length nstates indicating the number of states to be in each state aggregate (see Langrock and Zuchinni 2011).
#'@return A list containing model parameters, the stationary distribution, and
#'the AICc
#'@include Distributions.R
#'@include move.HSMM.pw2pn.R
#'@include move.HSMM.mllk.R
#'@include move.HSMM.mllk.full.R
#'@examples \dontrun{
#'######2 state 2 distribution with Poisson dwell time distribution
#'lmean=c(-3,-1) #meanlog parameters
#'sd=c(1,1) #sdlog parameters
#'rho<-c(0.2,0.3) # wrapped normal concentration parameters
#'mu<-c(pi,0) # wrapped normal mean parameters
#'dists=c("shiftpois","lognormal","wrpcauchy")
#'turn=c(1,2)
#'params=vector("list",3)
#'params[[1]]=matrix(c(2,9),nrow=2)
#'params[[2]]=cbind(lmean,sd)
#'params[[3]]=cbind(mu,rho)
#'nstates=2
#'obs=move.HSMM.simulate(dists,params,1000,nstates)$obs
#'turn=c(1,2)
#'move.HSMM=move.HSMM.mle(obs,dists,params,stepm=35,CI=F,iterlim=100,turn=turn,m1=c(30,30))
#'#Assess fit
#'xlim=matrix(c(0.001,-pi,2,pi),ncol=2)
#'breaks=c(200,20)
#'HSMM.plot(move.HSMM,xlim,breaks=breaks)
#'move.HSMM.psresid(move.HSMM)
#'move.HSMM.Altman(move.HSMM)
#'move.HSMM.dwell.plot(move.HSMM)
#'move.HSMM.ACF(move.HSMM,simlength=10000)
#'#Get bootstrap CIs
#'move.HSMM=move.HSMM.CI(move.HSMM,CI="boot",alpha=0.05,B=100,cores=4,stepm=5,iterlim=100)
#'
#'#2 state, 1 distribution shifted lognormal with Negative Binomial dwell time distribution 
#'mlog=c(-4.2,-2.2) #meanlog parameters
#'sdlog=c(1,1) #sdlog parameters
#'size=c(3,1)
#'prob=c(0.8,0.2)
#'shift=c(0.0000001,0.03)
#'dists=c("shiftnegbin","lnorm3")
#'stationary="yes"
#'params=vector("list",2)
#'params[[1]]=cbind(size,prob)
#'params[[2]]=cbind(mlog,sdlog,shift)
#'nstates=2
#'obs=move.HSMM.simulate(dists,params,5000,nstates)$obs
#'move.HSMM=move.HSMM.mle(obs,dists,params,stepm=35,iterlim=200,m1=c(30,30),)
#'#Assess fit
#'xlim=matrix(c(0.001,0.8),ncol=2)
#'breaks=c(200)
#'HSMM.plot(move.HSMM,xlim,breaks=breaks)
#'move.HSMM.psresid(move.HSMM)
#'move.HSMM.Altman(move.HSMM)
#'move.HSMM.dwell.plot(move.HSMM)
#'move.HSMM.ACF(move.HSMM,simlength=10000)
#'#Get bootstrap CIs
#'move.HSMM=move.HSMM.CI(move.HSMM,CI="boot",alpha=0.05,B=100,cores=4,stepm=5,iterlim=100)
#'
#'#2 states, 3 dist-lognorm, wrapped cauchy, poisson
#'#For example, this could be movement path lengths, turning angles,
#'#and accelerometer counts from a GPS-collared animal.
#'lmean=c(-3,-1) #meanlog parameters
#'sd=c(1,1) #sdlog parameters
#'rho<-c(0.2,0.3) # wrapped Cauchy concentration parameters
#'mu<-c(pi,0) # wrapped Cauchy mean parameters
#'lambda=c(2,20)
#'dists=c("shiftpois","lognormal","wrpcauchy","poisson")
#'turn=c(1,2)
#'params=vector("list",4)
#'params[[1]]=matrix(c(2,9),nrow=2)
#'params[[2]]=cbind(lmean,sd)
#'params[[3]]=cbind(mu,rho)
#'params[[4]]=matrix(lambda,ncol=1)
#'nstates=2
#'stationary="yes"
#'obs=move.HSMM.simulate(dists,params,5000,nstates)$obs
#'move.HSMM=move.HSMM.mle(obs,dists,params,stepm=35,iterlim=200,m1=c(30,30),turn=turn)
#'#Assess fit.
#'xlim=matrix(c(0.001,-pi,0,2,pi,40),ncol=2)
#'by=c(0.001,0.001,1)
#'breaks=c(200,20,20)
#'HSMM.plot(move.HSMM,xlim,breaks=breaks)
#'move.HSMM.psresid(move.HSMM)
#'move.HSMM.Altman(move.HSMM)
#'move.HSMM.dwell.plot(move.HSMM)
#'move.HSMM.ACF(move.HSMM,simlength=10000)
#'#Get bootstrap CIs
#'move.HSMM=move.HSMM.CI(move.HSMM,CI="boot",alpha=0.05,B=100,cores=4,stepm=5,iterlim=100)
#'
#'######3 state 2 distribution with Poisson dwell time distribution
#'lmean=c(-3,-2,-1) #meanlog parameters
#'sd=c(1,1,1) #sdlog parameters
#'rho<-c(0.2,0.3,0.4) # wrapped normal concentration parameters
#'mu<-c(pi,0,0) # wrapped normal mean parameters
#'gamma0=matrix(c(0,0.2,0.8,0.6,0,0.4,0.5,0.5,0),byrow=T,nrow=3)
#'dists=c("shiftpois","lognormal","wrpcauchy")
#'nstates=3
#'turn=c(1,2,2)
#'stationary="yes"
#'params=vector("list",4)
#'params[[1]]=gamma0
#'params[[2]]=matrix(c(2,4,9),nrow=3)
#'params[[3]]=cbind(lmean,sd)
#'params[[4]]=cbind(mu,rho)
#'obs=move.HSMM.simulate(dists,params,5000,nstates)$obs
#'turn=c(1,2,2)
#'move.HSMM=move.HSMM.mle(obs,dists,params,stepm=35,CI=F,iterlim=100,turn=turn,m1=c(30,30,30))
#'#Assess fit
#'xlim=matrix(c(0.001,-pi,2,pi),ncol=2)
#'breaks=c(200,20)
#'HSMM.plot(move.HSMM,xlim,breaks=breaks)
#'move.HSMM.psresid(move.HSMM)
#'move.HSMM.Altman(move.HSMM)
#'move.HSMM.dwell.plot(move.HSMM)
#'move.HSMM.ACF(move.HSMM,simlength=10000)
#'#Get bootstrap CIs
#'move.HSMM=move.HSMM.CI(move.HSMM,CI="boot",alpha=0.05,B=100,cores=4,stepm=5,iterlim=100)
#'}
#'@export
move.HSMM.mle <- function(obs,dists,params,stepm=5,CI=F,iterlim=150,turn=NULL,m1,alpha=0.05,B=100,cores=4){
  #check input
  nstates=nrow(params[[length(params)]])
  if(length(m1)!=nstates)stop("length(m1) must equal nstates")
  if(is.matrix(obs)==F&is.data.frame(obs)==F)stop("argument 'obs' must be a ndist x n matrix or data frame")
  if(!all(unlist(lapply(params,is.matrix))))stop("argument 'params' must contain nstate x nparam matrices")
  #if(any(rowSums(params[[1]])!=1))stop("Transition matrix rows should sum to 1")
  if(!all(nrow(params[[1]])-unlist(lapply(params,nrow))==0))stop("All parameter matrices must have the same number of rows")
  dwelldists=c("pospois","posnegbin","posgeom","logarithmic","shiftnegbin","shiftpois")
  if(sum(dists[1]==dwelldists)==0)stop("The first distribution must be the dwell time distribution")
  nstates=nrow(params[[1]])
  ndists=length(dists)
  if(nstates==1)stop("A 1 state HSMM does not make sense")
  if((nstates>2)&(length(params)==(ndists)))stop("Must include tpm in params when nstates>2")
  if((nstates==2)&(length(params)==(ndists+1)))stop("Don't include tpm in params when nstate=2")
  if(ncol(obs)!=(length(dists)-1))stop("Number of columns in obs much match number of observation distributions")
  out=Distributions(dists,nstates,turn)
  if(nstates==2){
    if(!all(unlist(lapply(params,ncol))==out[[7]]))stop("Incorrect number of parameters supplied for at least 1 distribution.")
  }else if(nstates>2){
    params2=params
    params2[[1]]=NULL
    if(!all(unlist(lapply(params2,ncol))==out[[7]]))stop("Incorrect number of parameters supplied for at least 1 distribution.")
  }
  if(any(is.element(dists,c("wrpnorm","wrpcauchy")))){
    if(is.null(turn))stop("Must input turn")
    if(length(turn)!=nstates)stop("Number of turn elements must = number of hidden states")
  }
  if(!(any(is.element(dists,c("wrpnorm","wrpcauchy"))))&(!is.null(turn)))stop("No turn argument needed--no circular distribution.")
  #Get appropriate linearizing transformations and PDFs 
  transforms=out[[1]]
  inv.transforms=out[[2]]
  PDFs=out[[3]]
  CDFs=out[[4]]
  skeleton=params
  #transform parameters
  parvect <- move.HSMM.pn2pw(transforms,params,nstates)
  #maximize likelihood.
  #try starting with stationary distribution, may have problems inverting t.p.m to get stationary dist.
  mod <- try(nlm(move.HSMM.mllk,parvect,obs,print.level=2,stepmax=stepm,PDFs=PDFs,CDFs=CDFs,skeleton=skeleton,inv.transforms=inv.transforms,nstates=nstates,iterlim=iterlim,m1=m1,ini=0),silent=T)
  if(!is.null(attributes(mod)$condition)){
    cat('\n Cannot invert t.p.m.  Maximizing with equal state probabilities at time 1.')
    #If that doesn't work, start with 1/nstates for all states
    mod <- nlm(move.HSMM.mllk,parvect,obs,print.level=2,stepmax=stepm,PDFs=PDFs,CDFs=CDFs,skeleton=skeleton,inv.transforms=inv.transforms,nstates=nstates,iterlim=iterlim,m1=m1,ini=1)
    #Then use these MLEs as starting values and start with stationary distribution
    parvect=mod$estimate
    cat('\n Now maximizing starting with the stationary distribution.')
    mod <- nlm(move.HSMM.mllk,parvect,obs,print.level=2,stepmax=stepm,PDFs=PDFs,CDFs=CDFs,skeleton=skeleton,inv.transforms=inv.transforms,nstates=nstates,iterlim=iterlim,m1=m1,ini=0)  
  }
#   
#   if((stationary=="no")|(stationary=="both")){
#     mod <- nlm(move.HSMM.mllk,parvect,obs,print.level=2,stepmax=stepm,PDFs=PDFs,CDFs=CDFs,skeleton=skeleton,inv.transforms=inv.transforms,nstates=nstates,iterlim=iterlim,m1=m1,ini=1)
#   }
#   if((stationary=="both")|(stationary=="yes")){
#     if(stationary=="both"){
#       parvect=mod$estimate
#     }
#     mod <- nlm(move.HSMM.mllk,parvect,obs,print.level=2,stepmax=stepm,PDFs=PDFs,CDFs=CDFs,skeleton=skeleton,inv.transforms=inv.transforms,nstates=nstates,iterlim=iterlim,m1=m1,ini=0)  
#   }
  mllk <- -mod$minimum
  pn <- move.HSMM.pw2pn(inv.transforms,mod$estimate,skeleton,nstates)
  params=pn$params
  gamma=gen.Gamma(m1,params,PDFs,CDFs)
  delta <- solve(t(diag(sum(m1))-gamma+1),rep(1,sum(m1)))
  npar=length(parvect)
  AIC=2*npar-2*mllk
  AICc=AIC+(2*npar*(npar+1))/(nrow(obs)-npar-1)
  #Calculate approximate stationary distribution
  mstart=c(1,cumsum(m1)+1)
  mstart=mstart[-length(mstart)]
  mstop=cumsum(m1)
  delta2=rep(NA,length(m1))
  for(i in 1:length(m1)){
    delta2[i]=sum(delta[mstart[i]:mstop[i]])
  }
  delta2=cbind(delta2,rep(NA,nstates),rep(NA,nstates))
  state=seq(1,nstates,1)
  rownames(delta2)=paste("state",state)
  level=100*(1-alpha)
  colnames(delta2)=c("est.",paste(level,"% lower",sep=""),paste(level,"% upper",sep=""))
  #Get CIs
  if(CI!=FALSE){
    parout=cbind(unlist(pn),rep(NA,length(unlist(pn))),rep(NA,length(unlist(pn))))
    move=list(dists=dists,nstates=nstates,params=pn$params,parout=parout,delta=delta2,npar=npar,mllk=mllk,AICc=AICc,turn=turn,m1=m1,obs=obs)
    class(move)="move.HSMM"
    out=move.HSMM.CI(move,CI=CI,alpha=alpha,B=B,cores=cores,stepm=stepm,iterlim=iterlim)
    lower=out$parout[,2]
    upper=out$parout[,3]
    delta2[,2:3]=out$delta[,2:3]
  }else{
    upper=lower=rep(NA,length(unlist(pn)))
  }
  #build structure for parameter estimates and confidence intervals
  parout=cbind(unlist(pn),unlist(lower),unlist(upper))
  colnames(parout)=c("est.",paste(level,"% lower",sep=""),paste(level,"% upper",sep=""))
  par=1
  #Check nstates>2 code
  if(nstates>2){
    #which parameters are derived?
    design=matrix(rep(0,nstates*nstates),nrow=nstates)
    design[1:(nstates-1),nstates]=1
    design[nstates,nstates-1]=1
    fixed=which(as.numeric(t(design))==1)
    for(i in 1:nrow(params[[1]])){
      for(j in 1:nrow(params[[1]])){
        if(i==j){
          rownames(parout)[par]=paste("P(",i,"|",j,")**")
        }else{
          if(par%in%fixed){
            rownames(parout)[par]=paste("P(",i,"|",j,")*")
          }else{
            rownames(parout)[par]=paste("P(",i,"|",j,")")
          }
        }
        par=par+1
      }
    }
    params[[1]]=NULL
  }
  for(k in 1:ndists){
    for(j in 1:ncol(params[[k]])){
      for(i in 1:nrow(params[[k]])){
        rownames(parout)[par]=paste(dists[k],colnames(params[[k]])[j],i)
        if(CI==T){
          if(parout[par,2]>parout[par,3]){
            parout[par,2:3]=parout[par,3:2]
          }
        }
        par=par+1
      }
    }
  }
  out=list(dists=dists,nstates=nstates,params=pn$params,parout=parout,delta=delta2,npar=npar,mllk=mllk,AICc=AICc,turn=turn,m1=m1,obs=obs)
  class(out)="move.HSMM"
  out
}