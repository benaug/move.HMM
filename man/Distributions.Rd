\name{Distributions}
\alias{Distributions}
\title{Parameter transformations and distribution functions}
\usage{
  Distributions(dists, nstates, turn = NULL)
}
\arguments{
  \item{dists}{A length ndist vector of distributions from
  the following list: weibull, gamma, exponential, normal,
  lognormal, lnorm3, posnorm, invgamma, rayleigh, f, ncf,
  dagum, frechet, beta, geometric, logarithmic, binom,
  poisson, nbinom, zapois, wrpcauchy,
  wrpnorm,shiftpois,shiftnegbin}

  \item{nstates}{Number of hidden states}

  \item{turn}{Parameters determining the transformation for
  circular distributions. turn=1 yields support on (0,2pi)
  and turn=2 yields support on (-pi,pi).  For animal
  movement models, the "encamped" state should use turn=1
  and the "traveling" state should use turn=2.}
}
\value{
  A list containing transformations, inverse
  transformations, PDFs, CDFs, and functions to generate
  random deviates for the input distributions.
}
\description{
  This function contains supported PDFs, CDFs, and random
  deviate generating functions along with appropriate
  parameter transformations and inverse transformations
  needed for maximum likelihood estimation.  It takes a
  vector of distributions and the number of hidden states
  as input and outputs the appropriate functions for those
  distributions.  This function is for internal use.
  Additional distributions can be added by modifying this
  function.
}

