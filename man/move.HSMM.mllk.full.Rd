\name{move.HSMM.mllk.full}
\alias{move.HSMM.mllk.full}
\title{Compute negative log likelihood of HSMM using all parameters}
\usage{
  move.HSMM.mllk.full(parvect, obs, PDFs, CDFs, skeleton,
    nstates, m1, ini)
}
\arguments{
  \item{parvect}{The vector of parameters to be estimated}

  \item{obs}{A n x ndist matrix of data.  If ndist=1, obs
  must be a n x 1 matrix.}

  \item{PDFs}{A list of PDFs for the dwell time and ndist
  observation distributions.}

  \item{CDFs}{A list of CDFs for the dwell time and ndist
  observation distributions.}

  \item{skeleton}{A list with the original parameter
  structure used to reassemble parvect}

  \item{nstates}{Number of hidden states}

  \item{m1}{a vector of length nstates that specifies how
  many states will be used to approximate each state of the
  HSMM (see Langrock and Zuchinni 2011)}

  \item{ini}{numeric value that specifies how the initial
  state distribution is calculated. 0 sets the initial
  distribution to the stationary distribution.  If this
  matrix is not invertible, 1 sets the initial distribution
  for each state within each state agreggate to
  1/m(state).}
}
\value{
  The negative log likelihood of the hidden markov model.
}
\description{
  This function computes the negative log likelihood of the
  hidden Markov model. using all parameters, untransformed.
  It is used to get the covariance matrix of the fitted
  model.
}

