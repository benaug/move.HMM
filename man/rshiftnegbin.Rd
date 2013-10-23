\name{rshiftnegbin}
\alias{rshiftnegbin}
\title{Generate random deviates shifted negative binomial pdf}
\usage{
  rshiftnegbin(n, size, prob)
}
\arguments{
  \item{n}{the number of random deviates to generate}

  \item{size}{target for number of successful trials, or
  dispersion parameter (the shape parameter of the gamma
  mixing distribution). Must be strictly positive, need not
  be integer.}

  \item{prob}{probability of success in each trial. 0 <
  prob <= 1.}
}
\value{
  A vector of shifted negative binomial random deviates
}
\description{
  This function generates random deviates from the shifted
  negative binomial distribution. The shift is fixed to 1.
}

