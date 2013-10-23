\name{qlnorm3}
\alias{qlnorm3}
\title{Shifted lognormal quantile function}
\usage{
  qlnorm3(x, mlog, sdlog, x0 = 0)
}
\arguments{
  \item{x}{a vector of values where the quantile function
  is to be evaluated}

  \item{mlog}{mean of the distribution on the log scale}

  \item{sdlog}{standard deviation of the distribution on
  the log scale}

  \item{x0}{shift parameter (must be positive)}
}
\value{
  A vector of shifted lognormal quantile function values
}
\description{
  This function evaluates the shifted lognormal quantile
  function.  The shift must be positive.
}

