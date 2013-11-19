\name{dzanegbin}
\alias{dzanegbin}
\title{Zero-altered Negative Binomial PDF}
\usage{
  dzanegbin(x, size, prob = NULL, pobs0 = 0, log = FALSE)
}
\arguments{
  \item{x}{a vector of values where the pdf is to be
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
  pdf.  This function is modified from the dzanegbin VGAM
  package so that munb is not a parameter input option.
}

