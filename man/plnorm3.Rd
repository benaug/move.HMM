\name{plnorm3}
\alias{plnorm3}
\title{Shifted lognormal cdf}
\usage{
  plnorm3(x, mlog, sdlog, x0 = 0)
}
\arguments{
  \item{x}{a vector of values where the cdf is to be
  evaluated}

  \item{mlog}{mean of the distribution on the log scale}

  \item{sdlog}{standard deviation of the distribution on
  the log scale}

  \item{x0}{shift parameter (must be positive)}
}
\value{
  A vector of shifted lognormal cdf values
}
\description{
  This function evaluates the shifted lognormal cdf.  The
  shift must be positive.
}

