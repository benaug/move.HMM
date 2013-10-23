\name{move.HMM.viterbi}
\alias{move.HMM.viterbi}
\title{Assign states using the Viterbi algorithm}
\usage{
  move.HMM.viterbi(move.HMM, delta = NULL)
}
\arguments{
  \item{move.HMM}{A move.HMM object containing a fitted HMM
  model.}

  \item{delta}{An optional vector of starting state
  probabilities.  If no vector is supplied, the stationary
  distribution is used.}
}
\value{
  A vector of state assignments.
}
\description{
  This function, modified from Zucchini and MacDonald
  (2009), assigns states to observations using the Viterbi
  algorithm.  It takes as input a move.HMM object and an
  optional vector containing the starting state
  probabilities.
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
obs=move.HMM.simulate(dists,params,1000)$obs
turn=c(1,2)
move.HMM=move.HMM.mle(obs,dists,params,stepm=35,iterlim=100,turn=turn)
#get Viterbi state assignments
move.HMM.viterbi(move.HMM)
}
}

