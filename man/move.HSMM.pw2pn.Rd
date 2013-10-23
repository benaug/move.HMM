\name{move.HSMM.pw2pn}
\alias{move.HSMM.pw2pn}
\title{Transform parameters back to real scale}
\usage{
  move.HSMM.pw2pn(inv.transforms, parvect, skeleton,
    nstates)
}
\arguments{
  \item{parvect}{The vector of parameters to be estimated}

  \item{inv.transforms}{A list of inverse transformations
  used to transform parvect back to the original scale}

  \item{nstates}{Number of hidden states}

  \item{skeleton}{A list with the original parameter
  structure used to reassemble parvect}
}
\value{
  A list of parameter values on the real scale with an
  element for each distribution
}
\description{
  This function, modified from Langrock et al. (2012),
  transforms parameters back to the real scale.
}

