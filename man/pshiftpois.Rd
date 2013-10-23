\name{pshiftpois}
\alias{pshiftpois}
\title{Shifted Poisson cdf}
\usage{
  pshiftpois(x, lambda)
}
\arguments{
  \item{x}{a vector of values where the cdf is to be
  evaluated}

  \item{lambda}{vector of positive means (of an ordinary
  Poisson distribution)}
}
\value{
  A vector of shifted negative binomial cdf values
}
\description{
  This function evaluates the shifted Poisson cdf.  The
  shift is fixed to 1.
}

