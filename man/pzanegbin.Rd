\name{pzanegbin}
\alias{pzanegbin}
\title{Zero-altered Negative Binomial CDF}
\usage{
  pzanegbin(q, size, prob = NULL, pobs0 = 0)
}
\arguments{
  \item{q}{a vector of values where the cdf is to be
  evaluated}

  \item{size}{a value for the ordinary Negative Binomial
  size parameter}

  \item{prob}{a value for the ordinary Negative Binomial
  prob parameter}

  \item{pobs0}{Probability of zero.}
}
\value{
  A vector of zero-altered negative binomial pdf values
}
\description{
  This function evaluates zero-altered negative binomial
  cdf.  This function is modified from the pzanegbin VGAM
  package so that munb is not a parameter input option.
}

