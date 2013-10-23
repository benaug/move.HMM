\name{dwrpnorm}
\alias{dwrpnorm}
\title{Wrapped normal pdf}
\usage{
  dwrpnorm(theta, mu, rho, acc = 1e-05, tol = acc)
}
\arguments{
  \item{theta}{a vector of values where the pdf is to be
  evaluated}

  \item{mu}{a value for the location parameter in the
  interval [-pi,pi]}

  \item{rho}{a value for the concentration parameter in the
  interval [0,1]}

  \item{acc}{parameter defining the accuracy of the
  estimation of the density Terms are added to the infinite
  summation that defines the density function until
  successive estimates are within acc of each other.}

  \item{tol}{the same as acc}
}
\value{
  A vector of wrapped normal pdf values
}
\description{
  This function evaluates the wrapped normal pdf.  This
  function is modified from the wrapped normal pdf in the
  CircStats package so that it can evaluate more than one
  value at a time
}

