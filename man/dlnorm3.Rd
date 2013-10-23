\name{dlnorm3}
\alias{dlnorm3}
\title{Shifted lognormal pdf}
\usage{
  dlnorm3(x, mlog, sdlog, x0 = 0)
}
\arguments{
  \item{x}{a vector of values where the pdf is to be
  evaluated}

  \item{mlog}{mean of the distribution on the log scale}

  \item{sdlog}{standard deviation of the distribution on
  the log scale}

  \item{x0}{shift parameter (must be positive)}
}
\value{
  A vector of shifted lognormal pdf values
}
\description{
  This function evaluates the shifted lognormal pdf.  The
  shift must be positive.
}

