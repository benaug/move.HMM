\name{move.HMM.mle}
\alias{move.HMM.mle}
\title{Fit a Hidden Markov Model (HMM) via maximum likelihood}
\usage{
  move.HMM.mle(obs, dists, params, stepm = 35, CI = FALSE,
    iterlim = 150, turn = NULL, alpha = 0.05, B = 100,
    cores = 4, useRcpp = FALSE)
}
\arguments{
  \item{obs}{A n x ndist matrix of data.  If ndist=1, obs
  must be a n x 1 matrix. It is recommended that movement
  path distances are modeled at the kilometer scale rather
  than at the meter scale.  Turning angle observations must
  be in the interval [-pi,pi].}

  \item{dists}{A length d vector of distributions from the
  following list: weibull, gamma, exponential, normal,
  lognormal, lnorm3, posnorm, invgamma, rayleigh, f, ncf,
  dagum, frechet, beta, binom, poisson, nbinom, zapois,
  wrpcauchy, wrpnorm.  Note wrpnorm is much slower to
  evaluate than wrpcauchy. Differences in the amount of
  time taken to maximize can be substantial.}

  \item{params}{A list of length ndist+1 containing
  matrices of starting parameter values.  The first element
  of the list must be the starting values for the
  transition matrix.  If modeling a single behavioral
  state, the transition matrix must be a 1 x 1 matrix with
  value 1 (see example).  If any distributions only have 1
  parameter, the list entry must be a nstates x 1 matrix.
  Users should use reasonable starting values.  One method
  of finding good starting values is to plot randomly
  generated data from distributions with known parameter
  values and compare them to histograms of the data.  More
  complex models generally require better starting values.
  One strategy is to fit a one distribution model, then fit
  a two distribution model using the MLEs for the first
  distribution as starting values.  Signs that the model
  has converged to a local maximum are poor residual plots,
  poor plot.HMM plots, transition probabilities very close
  to 1 (especially for the shifted lognormal for which only
  a local maximum exists), and otherwise unreasonable
  parameter estimates.  Analysts should try multiple
  combinations of starting values to increase confidence
  that a global maximum has been achieved.}

  \item{stepm}{a positive scalar which gives the maximum
  allowable scaled step length. stepm is used to prevent
  steps which would cause the optimization function to
  overflow, to prevent the algorithm from leaving the area
  of interest in parameter space, or to detect divergence
  in the algorithm. stepm would be chosen small enough to
  prevent the first two of these occurrences, but should be
  larger than any anticipated reasonable step.  If
  maximization is failing due to the parameter falling
  outside of it's support, decrease stepm.}

  \item{CI}{A logical or character determining which type
  of CI is to be calculated.  If CI=FALSE, no CIs are
  calculated.  Otherwise, current options are "FD" for the
  finitie differences Hessian and "boot" for parametric
  bootstrapping and percentile CIs.}

  \item{stepm}{a positive scalar which gives the maximum
  allowable scaled step}

  \item{iterlim}{a positive integer specifying the maximum
  number of iterations to be performed before the nlm is
  terminated.}

  \item{turn}{Parameters determining the transformation for
  circular distributions. turn=1 leads to support on
  (0,2pi) and turn=2 leads to support on (-pi,pi).  For
  animal movement models, the "encamped" state should use
  turn=1 and the "traveling" state should use turn=2.}

  \item{alpha}{Type I error rate for CIs.  alpha=0.05 for
  95 percent CIs}

  \item{B}{Number of bootstrap resamples}

  \item{cores}{Number of cores to be used in parallel
  bootstrapping}

  \item{useRcpp}{Logical indicating whether or not to use
  Rcpp.  Doing so leads to significant speedups in model
  fitting and obtaining CIs for longer time series, say of
  length 3000+. See this site for getting Rcpp working on
  windows:
  http://www.r-bloggers.com/installing-rcpp-on-windows-7-for-r-and-c-integration/}
}
\value{
  A list containing model parameters, the stationary
  distribution, and the AICc
}
\description{
  move.HMM.mle is used to fit HMMs, allowing for multiple
  observation variables with different distributions
  (ndist=number of distributions). Maximization is
  performed in nlm.
}
\examples{
\dontrun{
#2 states, 2 dist-lognorm, wrapped normal
lmean=c(-3,-1) #meanlog parameters
sd=c(1,1) #sdlog parameters
rho<-c(0.2,0.3) # wrapped normal concentration parameters
mu<-c(pi,0) # wrapped normal mean parameters
gamma0=matrix(c(0.6,0.4,0.2,0.8),byrow=T,nrow=2)

dists=c("lognormal","wrpcauchy")
turn=c(1,2)
params=vector("list",3)
params[[1]]=gamma0
params[[2]]=cbind(lmean,sd)
params[[3]]=cbind(mu,rho)
obs=move.HMM.simulate(dists,params,1000)$obs
move.HMM=move.HMM.mle(obs,dists,params,stepm=35,iterlim=100,turn=turn)
#Assess fit
xlim=matrix(c(0.001,-pi,2,pi),ncol=2)
breaks=c(200,20)
HMM.plot(move.HMM,xlim,breaks=breaks)
move.HMM.psresid(move.HMM)
move.HMM.Altman(move.HMM)
move.HMM.dwell.plot(move.HMM)
move.HMM.ACF(move.HMM)
#Get bootstrap CIs (should usually use B>100)
move.HMM=move.HMM.CI(move.HMM,CI="boot",alpha=0.05,B=100,cores=4,stepm=35,iterlim=100)

#2 states, 1 dist-shifted lognormal
mlog=c(-4.2,-2.2) #meanlog parameters
sdlog=c(1,1) #sdlog parameters
shift=c(0.0000001,0.03)
gamma0=matrix(c(0.6,0.4,0.2,0.8),byrow=T,nrow=2)
dists=c("lnorm3")
params=vector("list",2)
params[[1]]=gamma0
params[[2]]=cbind(mlog,sdlog,shift)
obs=move.HMM.simulate(dists,params,1500)$obs
move.HMM=move.HMM.mle(obs,dists,params,stepm=35,iterlim=100)
#Assess fit
xlim=matrix(c(0.001,0.8),ncol=2)
breaks=c(200)
HMM.plot(move.HMM,xlim,breaks=breaks)
move.HMM.psresid(move.HMM)
move.HMM.Altman(move.HMM)
move.HMM.dwell.plot(move.HMM)
move.HMM.ACF(move.HMM)
#Get bootstrap CIs (should usually use B>100)
move.HMM=move.HMM.CI(move.HMM,CI="boot",alpha=0.05,B=100,cores=4,stepm=35,iterlim=100)

#1 state, 1 dist-gamma
shape=c(1) #shape parameters
rate=c(10) #rate parameters
gamma0=matrix(1,byrow=T,nrow=1)
dists=c("gamma")
params=vector("list",2)
params[[1]]=gamma0
params[[2]]=cbind(shape,rate)
obs=move.HMM.simulate(dists,params,1000)$obs
move.HMM=move.HMM.mle(obs,dists,params,stepm=35,iterlim=100)
#Assess fit
xlim=matrix(c(0.001,1),ncol=2)
breaks=c(20)
HMM.plot(move.HMM,xlim,breaks=breaks)
move.HMM.psresid(move.HMM)
move.HMM.Altman(move.HMM)
move.HMM.dwell.plot(move.HMM)
move.HMM.ACF(move.HMM)
#Get bootstrap CIs (should usually use B>100)
move.HMM=move.HMM.CI(move.HMM,CI="boot",alpha=0.05,B=100,cores=4,stepm=35,iterlim=100)

#1 state, 2 dist- Weibull, Wrapped cauchy
wshape=c(0.9) #Weibull shape parameter
wscale=c(0.3) #Weibull scale parameter
rho<-c(0.2) # wrapped cauchy concentration parameter
mu<-c(pi) # wrapped cauchy mean parameter
gamma0=matrix(1,byrow=T,nrow=1)
dists=c("weibull", "wrpcauchy")
turn=1
params=vector("list",3)
params[[1]]=gamma0
params[[2]]=cbind(wshape,wscale)
params[[3]]=cbind(mu, rho)
obs=move.HMM.simulate(dists,params,1000)$obs
move.HMM=move.HMM.mle(obs,dists,params,stepm=35,iterlim=100, turn=turn)
xlim=matrix(c(0.001,-pi,0,2,pi,40),ncol=2)
by=c(0.001,0.001,1)
breaks=c(200,20,20)
HMM.plot(move.HMM,xlim,breaks=breaks)
move.HMM.psresid(move.HMM)
move.HMM.Altman(move.HMM)
move.HMM.dwell.plot(move.HMM)
move.HMM.ACF(move.HMM)
#Get bootstrap CIs (should usually use B>100)
move.HMM=move.HMM.CI(move.HMM,CI="boot",alpha=0.05,B=100,cores=4,stepm=35,iterlim=100)

#3 states, 1 dist lognormal
lmean=c(-4,-2,-0.5) #shape parameters
sd=c(1,1,1) #rate parameters
gamma0=matrix(c(0.3,0.3,0.4,0.3,0.3,0.4,0.3,0.3,0.4),byrow=T,nrow=3)
dists=c("lognormal")
params=vector("list",2)
params[[1]]=gamma0
params[[2]]=cbind(lmean,sd)
obs=move.HMM.simulate(dists,params,10000)$obs
move.HMM=move.HMM.mle(obs,dists,params,stepm=35,iterlim=200)
#Assess fit
xlim=matrix(c(0.001,1),ncol=2)
breaks=c(200)
HMM.plot(move.HMM,xlim,breaks=breaks)
move.HMM.psresid(move.HMM)
move.HMM.Altman(move.HMM)
move.HMM.dwell.plot(move.HMM)
move.HMM.ACF(move.HMM)
#Get bootstrap CIs (should usually use B>100)
move.HMM=move.HMM.CI(move.HMM,CI="boot",alpha=0.05,B=100,cores=4,stepm=35,iterlim=100)

#2 states, 3 dist-lognorm, wrapped cauchy, poisson
#For example, this could be movement path lengths, turning angles,
#and accelerometer counts from a GPS-collared animal.
lmean=c(-3,-1) #meanlog parameters
sd=c(1,1) #sdlog parameters
rho<-c(0.2,0.3) # wrapped Cauchy concentration parameters
mu<-c(pi,0) # wrapped Cauchy mean parameters
lambda=c(2,20)
gamma0=matrix(c(0.6,0.4,0.2,0.8),byrow=T,nrow=2)
dists=c("lognormal","wrpcauchy","poisson")
turn=c(1,2)
params=vector("list",4)
params[[1]]=gamma0
params[[2]]=cbind(lmean,sd)
params[[3]]=cbind(mu,rho)
params[[4]]=matrix(lambda,ncol=1)
obs=move.HMM.simulate(dists,params,1500)$obs
move.HMM=move.HMM.mle(obs,dists,params,stepm=35,iterlim=150,turn=turn)
#Assess fit - note, not great--need more data.
xlim=matrix(c(0.001,-pi,0,2,pi,40),ncol=2)
by=c(0.001,0.001,1)
breaks=c(200,20,20)
HMM.plot(move.HMM,xlim,breaks=breaks)
move.HMM.psresid(move.HMM)
move.HMM.Altman(move.HMM)
move.HMM.dwell.plot(move.HMM)
move.HMM.ACF(move.HMM)
#Get bootstrap CIs (should usually use B>100)
move.HMM=move.HMM.CI(move.HMM,CI="boot",alpha=0.05,B=100,cores=4,stepm=35,iterlim=100)
}
}

