\name{move.HSMM.mle}
\alias{move.HSMM.mle}
\title{Fit a Hidden Semi Markov Model (HSMM) via maximum likelihood}
\usage{
  move.HSMM.mle(obs, dists, params, stepm = 5, CI = F,
    iterlim = 150, turn = NULL, m1, alpha = 0.05, B = 100,
    cores = 4, useRcpp = FALSE)
}
\arguments{
  \item{obs}{A n x ndist matrix of data.  If ndist=1, obs
  must be a n x 1 matrix. It is recommended that movement
  path distances are modeled at the kilometer scale rather
  than at the meter scale.}

  \item{dists}{A length d vector of distributions.  The
  first distribution must be the dwell time distribution
  chosen from the following list: shiftpois,shifnegbin. The
  subsequent distributions are for observation variables
  and must be chosen from the following list: weibull,
  gamma, exponential, normal, lognormal, lnorm3, posnorm,
  invgamma, rayleigh, f, ncf, dagum, frechet, beta, binom,
  poisson, nbinom, zapois, wrpcauchy, wrpnorm.  Note
  wrpnorm is much slower to evaluate than wrpcauchy.
  Differences in the amount of time taken to maximize can
  be substantial.}

  \item{params}{A list containing matrices of starting
  parameter values.  The structure of this list differs
  between 2 state models and models with greater than 2
  states.  For 2 state models, the first element of the
  list must be the starting values for the dwell time
  distribution.  For 3+ state models, the first element of
  the list must be the starting values for the t.p.m.
  Since dwell times are modeled explicitly, the diagonal of
  the t.p.m must be all zeros.  Further, rows should still
  sum to 1.  The second element of the list will now be the
  dwell time distribution parameter starting values.  See
  the example below for a demonstration.  As before, if any
  distributions only have 1 parameter, the list entry must
  be a nstates x 1 matrix.  Users should use reasonable
  starting values.  One method of finding good starting
  values is to plot randomly-generated data from
  distributions with known values and compare them to
  histograms of the data.  More complex models generally
  require better starting values.  One strategy for fitting
  HSMMs is to first fit a HMM and use the MLEs as starting
  values for the HSMM.  Signs that the model has converged
  to a local maximum are poor residual plots, poor plot.HMM
  plots, transition probabilities very close to 1
  (especially for the shifted lognormal for which only a
  local maximum exists), and otherwise unreasonable
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
  larger than any anticipated reasonable step.}

  \item{iterlim}{a positive integer specifying the maximum
  number of iterations to be performed before the nlm is
  terminated.}

  \item{turn}{Parameters determining the transformation for
  circular distributions. turn=1 leads to support on
  (0,2pi) and turn=2 leads to support on (-pi,pi).  For
  animal movement models, the "encamped" state should use
  turn=1 and the "traveling" state should use turn=2.}

  \item{CI}{A character determining which type of CI is to
  be calculated.  Current options are "FD" for the finitie
  differences Hessian and "boot" for parametric
  bootstrapping and percentile CIs.}

  \item{alpha}{Type I error rate for CIs.  alpha=0.05 for
  95 percent CIs}

  \item{B}{Number of bootstrap resamples}

  \item{cores}{Number of cores to be used in parallell
  bootstrapping}

  \item{m1}{vector of length nstates indicating the number
  of states to be in each state aggregate (see Langrock and
  Zuchinni 2011).}

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
  move.HSMM.mle is used to fit HSMMs, allowing for multiple
  observation variables with different distributions
  (ndist=number of distributions).  It approximates a HSMM
  with a HMM as outlined in Langrock et al. (2012).
  Maximization is performed in nlm. Currently only 2 hidden
  states are supported but this will be extended to an
  arbitrany number of states >1.
}
\examples{
\dontrun{
######2 state 2 distribution with Poisson dwell time distribution
lmean=c(-3,-1) #meanlog parameters
sd=c(1,1) #sdlog parameters
rho<-c(0.2,0.3) # wrapped normal concentration parameters
mu<-c(pi,0) # wrapped normal mean parameters
dists=c("shiftpois","lognormal","wrpcauchy")
turn=c(1,2)
params=vector("list",3)
params[[1]]=matrix(c(2,9),nrow=2)
params[[2]]=cbind(lmean,sd)
params[[3]]=cbind(mu,rho)
nstates=2
obs=move.HSMM.simulate(dists,params,1000,nstates)$obs
turn=c(1,2)
move.HSMM=move.HSMM.mle(obs,dists,params,stepm=35,CI=F,iterlim=100,turn=turn,m1=c(30,30))
#Assess fit
xlim=matrix(c(0.001,-pi,2,pi),ncol=2)
breaks=c(200,20)
HSMM.plot(move.HSMM,xlim,breaks=breaks)
move.HSMM.psresid(move.HSMM)
move.HSMM.Altman(move.HSMM)
move.HSMM.dwell.plot(move.HSMM)
move.HSMM.ACF(move.HSMM,simlength=10000)
#Get bootstrap CIs
move.HSMM=move.HSMM.CI(move.HSMM,CI="boot",alpha=0.05,B=100,cores=4,stepm=5,iterlim=100)

#2 state, 1 distribution shifted lognormal with Negative Binomial dwell time distribution
mlog=c(-4.2,-2.2) #meanlog parameters
sdlog=c(1,1) #sdlog parameters
size=c(3,1)
prob=c(0.8,0.2)
shift=c(0.0000001,0.03)
dists=c("shiftnegbin","lnorm3")
stationary="yes"
params=vector("list",2)
params[[1]]=cbind(size,prob)
params[[2]]=cbind(mlog,sdlog,shift)
nstates=2
obs=move.HSMM.simulate(dists,params,5000,nstates)$obs
move.HSMM=move.HSMM.mle(obs,dists,params,stepm=35,iterlim=200,m1=c(30,30),)
#Assess fit
xlim=matrix(c(0.001,0.8),ncol=2)
breaks=c(200)
HSMM.plot(move.HSMM,xlim,breaks=breaks)
move.HSMM.psresid(move.HSMM)
move.HSMM.Altman(move.HSMM)
move.HSMM.dwell.plot(move.HSMM)
move.HSMM.ACF(move.HSMM,simlength=10000)
#Get bootstrap CIs
move.HSMM=move.HSMM.CI(move.HSMM,CI="boot",alpha=0.05,B=100,cores=4,stepm=5,iterlim=100)

#2 states, 3 dist-lognorm, wrapped cauchy, poisson
#For example, this could be movement path lengths, turning angles,
#and accelerometer counts from a GPS-collared animal.
lmean=c(-3,-1) #meanlog parameters
sd=c(1,1) #sdlog parameters
rho<-c(0.2,0.3) # wrapped Cauchy concentration parameters
mu<-c(pi,0) # wrapped Cauchy mean parameters
lambda=c(2,20)
dists=c("shiftpois","lognormal","wrpcauchy","poisson")
turn=c(1,2)
params=vector("list",4)
params[[1]]=matrix(c(2,9),nrow=2)
params[[2]]=cbind(lmean,sd)
params[[3]]=cbind(mu,rho)
params[[4]]=matrix(lambda,ncol=1)
nstates=2
stationary="yes"
obs=move.HSMM.simulate(dists,params,5000,nstates)$obs
move.HSMM=move.HSMM.mle(obs,dists,params,stepm=35,iterlim=200,m1=c(30,30),turn=turn)
#Assess fit.
xlim=matrix(c(0.001,-pi,0,2,pi,40),ncol=2)
by=c(0.001,0.001,1)
breaks=c(200,20,20)
HSMM.plot(move.HSMM,xlim,breaks=breaks)
move.HSMM.psresid(move.HSMM)
move.HSMM.Altman(move.HSMM)
move.HSMM.dwell.plot(move.HSMM)
move.HSMM.ACF(move.HSMM,simlength=10000)
#Get bootstrap CIs
move.HSMM=move.HSMM.CI(move.HSMM,CI="boot",alpha=0.05,B=100,cores=4,stepm=5,iterlim=100)

######3 state 2 distribution with Poisson dwell time distribution
lmean=c(-3,-2,-1) #meanlog parameters
sd=c(1,1,1) #sdlog parameters
rho<-c(0.2,0.3,0.4) # wrapped normal concentration parameters
mu<-c(pi,0,0) # wrapped normal mean parameters
gamma0=matrix(c(0,0.2,0.8,0.6,0,0.4,0.5,0.5,0),byrow=T,nrow=3)
dists=c("shiftpois","lognormal","wrpcauchy")
nstates=3
turn=c(1,2,2)
stationary="yes"
params=vector("list",4)
params[[1]]=gamma0
params[[2]]=matrix(c(2,4,9),nrow=3)
params[[3]]=cbind(lmean,sd)
params[[4]]=cbind(mu,rho)
obs=move.HSMM.simulate(dists,params,5000,nstates)$obs
turn=c(1,2,2)
move.HSMM=move.HSMM.mle(obs,dists,params,stepm=35,CI=F,iterlim=100,turn=turn,m1=c(30,30,30))
#Assess fit
xlim=matrix(c(0.001,-pi,2,pi),ncol=2)
breaks=c(200,20)
HSMM.plot(move.HSMM,xlim,breaks=breaks)
move.HSMM.psresid(move.HSMM)
move.HSMM.Altman(move.HSMM)
move.HSMM.dwell.plot(move.HSMM)
move.HSMM.ACF(move.HSMM,simlength=10000)
#Get bootstrap CIs
move.HSMM=move.HSMM.CI(move.HSMM,CI="boot",alpha=0.05,B=100,cores=4,stepm=5,iterlim=100)
}
}

