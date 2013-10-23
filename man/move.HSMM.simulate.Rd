\name{move.HSMM.simulate}
\alias{move.HSMM.simulate}
\title{Simulate data from a hidden semi markov model}
\usage{
  move.HSMM.simulate(dists, params, n, nstates,
    delta = NULL)
}
\arguments{
  \item{dists}{A vector of distributions of length ndist+1,
  with the first distribution chosen from: shiftpois,
  shiftnegbin The following distribution(s) are from the
  following list: weibull, gamma, exponential, normal,
  lognormal, lnorm3, posnorm, invgamma, rayleigh, f, ncf,
  dagum, frechet, beta, binom, poisson, nbinom, zapois,
  wrpcauchy, wrpnorm.}

  \item{params}{a list of starting parameter values.  If
  nstates=2, params is of length ndist+1 The first element
  of the list must be the starting values for the dwell
  time distribution.  If any distributions only have 1
  parameter, the list entry must be a nstates x 1 matrix.
  If nstates>2, params is of lenght ndist+2, with the first
  element being the starting values for the t.p.m., the
  second element being the starting values for the dwell
  time distribution, and the following element(s) being the
  starting values for the observation distributions}

  \item{n}{The length of the hidden markov model to be
  simulated}

  \item{nstates}{The number of hidden states}

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
  semi markov model
}

