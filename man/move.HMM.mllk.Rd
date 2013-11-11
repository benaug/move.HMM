\name{move.HMM.mllk}
\alias{move.HMM.mllk}
\title{Compute negative log likelihood of HMM}
\usage{
  move.HMM.mllk(parvect, obs, PDFs, skeleton,
    inv.transforms, nstates, useRcpp = FALSE)
}
\arguments{
  \item{parvect}{The vector of parameters to be estimated}

  \item{obs}{A n x ndist matrix of data.  If ndist=1, obs
  must be a n x 1 matrix.}

  \item{PDFs}{A list of PDFs for the ndist distributions.}

  \item{skeleton}{A list with the original parameter
  structure used to reassemble parvect}

  \item{inv.transforms}{A list of inverse transformations
  used to transform parvect back to the original scale}

  \item{nstates}{Number of hidden states}

  \item{useRcpp}{Logical indicating whether or not to use
  Rcpp.}
}
\value{
  The negative log likelihood of the hidden markov model.
}
\description{
  This function, modified from Langrock et al. (2012),
  computes the negative log likelihood of the hidden Markov
  model.
}

