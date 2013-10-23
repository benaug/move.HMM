\name{move.HMM.pn2pw}
\alias{move.HMM.pn2pw}
\title{Transform parameters to the working scale}
\usage{
  move.HMM.pn2pw(transforms, params, nstates)
}
\arguments{
  \item{params}{A list of length d+1 containing matrices of
  starting parameter values.  The first element of the list
  must be the starting values for the transition matrix.
  If any distributions only have 1 parameter, the list
  entry must be a nstates x 1 matrix.}

  \item{transforms}{A list of transformations used to
  transform the parameters to the working scale for
  maximization.}

  \item{nstates}{Number of hidden states}
}
\value{
  A vector of transformed parameters
}
\description{
  This function, modified from Langrock et al. (2012),
  transforms parameters to the working scale for
  maximization.
}

