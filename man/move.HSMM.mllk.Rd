\name{move.HSMM.mllk}
\alias{move.HSMM.mllk}
\title{Compute negative log likelihood of HSMM}
\usage{
  move.HSMM.mllk(parvect, obs, PDFs, CDFs, skeleton,
    inv.transforms, nstates, m1, ini)
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

  \item{inv.transforms}{A list of inverse transformations
  used to transform parvect back to the original scale}

  \item{nstates}{Number of hidden states}

  \item{m1}{a vector of length nstates that specifies how
  many states will be used to approximate each state of the
  HSMM (see Langrock and Zucchini 2011)}

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
  This function, modified from Langrock et al. (2012),
  computes the negative log likelihood of the hidden Markov
  model.
}

