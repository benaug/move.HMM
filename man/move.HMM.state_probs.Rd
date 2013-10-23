\name{move.HMM.state_probs}
\alias{move.HMM.state_probs}
\title{Compute conditional and posterior state probabilities}
\usage{
  move.HMM.state_probs(move.HMM)
}
\arguments{
  \item{move.HMM}{A move.HMM object containing a fitted HMM
  model.}
}
\value{
  A list of conditional and posterior state probabilities.
  #'
}
\description{
  This function, modified from Zucchini and MacDonald
  (2009), computes the conditional state probabilities and
  posterior state probabilities (Patterson et al. 2009)
  using the stationary distribution as the initial
  distribution. It takes as input a move.HMM object.
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
turn=c(1,2)
move.HMM=move.HMM.mle(obs,dists,params,stepm=35,iterlim=100,turn=turn)
#get conditional and posterior state probabilities
move.HMM.state_probs(move.HMM)
}
}

