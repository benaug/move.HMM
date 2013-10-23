\name{rshiftpois}
\alias{rshiftpois}
\title{Generate random deviates shifted Poisson pdf}
\usage{
  rshiftpois(n, lambda)
}
\arguments{
  \item{n}{the number of random deviates to generate}

  \item{lambda}{vector of positive means (of an ordinary
  Poisson distribution)}
}
\value{
  A vector of shifted Poisson random deviates
}
\description{
  This function generates random deviates from the shifted
  Poisson distribution. The shift is fixed to 1.
}

