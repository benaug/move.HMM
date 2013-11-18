\name{rwrpcauchy}
\alias{rwrpcauchy}
\title{Wrapped cauchy random number generation}
\usage{
  rwrpcauchy(n, mu = 0, rho = exp(-1))
}
\arguments{
  \item{n}{The number of random numbers to generate}

  \item{mu}{A value for the mu parameter}

  \item{rho}{A value for the concentration parameter in the
  interval [0,1]}
}
\value{
  A vector of random numbers drawn from a wrapped cauchy
  distribution
}
\description{
  This function generates random number from the wrapped
  cauchy distribution. It is modified from the function in
  the CircStats package so that it can evaluate multiple
  parameter combinations in the same call.
}

