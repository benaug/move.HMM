\name{pwrpcauchy}
\alias{pwrpcauchy}
\title{Wrapped cauchy cdf}
\usage{
  pwrpcauchy(q, mu = 0, rho = 0)
}
\arguments{
  \item{q}{a vector of values where the cdf is to be
  evaluated}

  \item{mu}{A value for the location parameter}

  \item{rho}{A value for the concentration parameter in the
  interval [0,1]}
}
\value{
  A vector of wrapped cauchy cdf values
}
\description{
  This function evaluates the wrapped cauchy cdf.
  Observations must be in the interval [-pi,pi] or [0,2pi].
}

