\name{move.HSMM.Altman}
\alias{move.HSMM.Altman}
\title{Altman plot}
\usage{
  move.HSMM.Altman(move.HSMM)
}
\arguments{
  \item{move.HSMM}{A move.HSMM object containing a fitted
  HSMM model.}
}
\value{
  A plot
}
\description{
  This function provides a plot from Altman (2004) that
  allows the assessment of goodness-of-fit.  Assuming the
  observed process is stationary, it plots the marginal CDF
  of the data evaluated at the maximum likelihood estimates
  against the empirical CDF values.  The plot should be
  close to a 45 degree straight line.
}

