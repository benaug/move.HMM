\name{dshiftnegbin}
\alias{dshiftnegbin}
\title{Shifted negative binomial pdf}
\usage{
  dshiftnegbin(x, size, prob)
}
\arguments{
  \item{x}{a vector of values where the pdf is to be
  evaluated}

  \item{size}{target for number of successful trials, or
  dispersion parameter (the shape parameter of the gamma
  mixing distribution). Must be strictly positive, need not
  be integer.}

  \item{prob}{probability of success in each trial. 0 <
  prob <= 1.}
}
\value{
  A vector of shifted negative binomial pdf values
}
\description{
  This function evaluates the shifted negative binomial
  pdf.  The shift is fixed to 1.
}

