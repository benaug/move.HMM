\name{dshiftpois}
\alias{dshiftpois}
\title{Shifted Poisson pdf}
\usage{
  dshiftpois(x, lambda)
}
\arguments{
  \item{x}{a vector of values where the pdf is to be
  evaluated}

  \item{lambda}{vector of positive means (of an ordinary
  Poisson distribution)}
}
\value{
  A vector of shifted negative binomial pdf values
}
\description{
  This function evaluates the shifted Poisson pdf.  The
  shift is fixed to 1.
}

