\name{move.HMM.mllk.full}
\alias{move.HMM.mllk.full}
\title{Compute negative log likelihood of HMM using all parameters}
\usage{
  move.HMM.mllk.full(parvect, obs, PDFs, skeleton, nstates,
    useRcpp = FALSE)
}
\arguments{
  \item{parvect}{The vector of parameters to be estimated}

  \item{obs}{A n x ndist matrix of data.  If ndist=1, obs
  must be a n x 1 matrix.}

  \item{PDFs}{A list of PDFs for the ndist distributions.}

  \item{skeleton}{A list with the original parameter
  structure used to reassemble parvect}

  \item{nstates}{Number of hidden states}

  \item{useRcpp}{Logical indicating whether or not to use
  Rcpp.}
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

