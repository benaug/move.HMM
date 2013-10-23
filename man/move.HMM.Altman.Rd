\name{move.HMM.Altman}
\alias{move.HMM.Altman}
\title{Altman plot}
\usage{
  move.HMM.Altman(move.HMM)
}
\arguments{
  \item{move.HMM}{A move.HMM object containing a fitted HMM
  model.}
}
\value{
  Altman plot
}
\description{
  This function provides a plot from Altman (2004) that
  allows the assessment of goodness-of-fit.  Assuming the
  observed process is stationary, it plots the marginal CDF
  of the data evaluated at the maximum likelihood estimates
  against the empirical CDF values.  The plotted values
  should be follow a 45 degree straight line.
}

