\name{rwrpnorm}
\alias{rwrpnorm}
\title{Wrapped normal random number generation}
\usage{
  rwrpnorm(n, mu, rho)
}
\arguments{
  \item{n}{The number of random numbers to generate}

  \item{mu}{A value for the mu parameter}

  \item{rho}{A value for the concentration parameter in the
  interval [0,1]}
}
\value{
  A vector of random numbers drawn from a wrapped normal
  distribution
}
\description{
  This function generates random numbers from the wrapped
  normal distribution. It is modified from the functoin in
  the CircStats package so that it can evaluate multiple
  parameter combinations in the same call.
}

