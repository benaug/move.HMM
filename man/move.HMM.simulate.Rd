\name{move.HMM.simulate}
\alias{move.HMM.simulate}
\title{Simulate data from a hidden markov model}
\usage{
  move.HMM.simulate(dists, params, n, delta = NULL)
}
\arguments{
  \item{dists}{A length ndist vector of distributions from
  the following list: weibull, gamma, exponential, normal,
  lognormal, lnorm3, posnorm, invgamma, rayleigh, f, ncf,
  dagum, frechet, beta, binom, poisson, nbinom, zapois,
  wrpcauchy, wrpnorm}

  \item{params}{a list of length ndist+1 containing
  matrices of starting parameter values.  The first element
  of the list must be the starting values for the
  transition matrix.  If any distributions only have 1
  parameter, the list entry must be a nstates x 1 matrix.}

  \item{n}{The length of the hidden markov model to be
  simulated}

  \item{delta}{an optional vector of starting state
  probabilities.  If no vector is supplied, the stationary
  distribution is used.}
}
\value{
  The observations from the simulated HMM and the true
  state sequence
}
\description{
  This function simulates data from a user-specified hidden
  markov model
}
\examples{
\dontrun{
#2 states, 2 dist-lognorm, wrapped normal
lmean=c(-3,-1) #meanlog parameters
sd=c(1,1) #sdlog parameters
rho<-c(0.2,0.3) # wrapped normal concentration parameters
mu<-c(pi,0) # wrapped normal mean parameters
gamma0=matrix(c(0.6,0.4,0.2,0.8),byrow=T,nrow=2)

dists=c("lognormal","wrpnorm")
nstates=2
turn=c(1,2)
params=vector("list",3)
params[[1]]=gamma0
params[[2]]=cbind(lmean,sd)
params[[3]]=cbind(mu,rho)
obs=move.HMM.simulate(dists,params,1000)
#see more examples in move.HMM.mle
}
}

