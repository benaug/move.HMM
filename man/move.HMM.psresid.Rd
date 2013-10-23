\name{move.HMM.psresid}
\alias{move.HMM.psresid}
\title{Compute and plot pseudoresiduals}
\usage{
  move.HMM.psresid(move.HMM, plots = T, returnresids = F)
}
\arguments{
  \item{move.HMM}{A move.HMM object containing a fitted HMM
  model.}

  \item{plots}{A logical indicating whether or not to
  produce plots}

  \item{returnresids}{A logical indicating whether or not
  to return the pseudo-residuals}
}
\value{
  A list of normal pseudoresiduals of length ndist.
}
\description{
  This function, modified from Zuchinni and MacDonald
  (2009), plots the ordinary normal pseudoresiduals.  Three
  plots are produced-a qq-plot, normal pseudoresiduals
  plotted through time, and an ACF plot.  Pseudoresiduals
  for integer valued distributions need to be fixed.
}

