\name{rlnorm3}
\alias{rlnorm3}
\title{Shifted lognormal random number generation}
\usage{
  rlnorm3(n, mlog, sdlog, x0 = 0)
}
\arguments{
  \item{n}{the number of random numbers to generate}

  \item{mlog}{mean of the distribution on the log scale}

  \item{sdlog}{standard deviation of the distribution on
  the log scale}

  \item{x0}{shift parameter (must be positive)}
}
\value{
  A vector random numbers from a shifted lognormal
  distribution
}
\description{
  This function generates random numbers from the shifted
  lognormal distribution. The shift must be positive.
}

