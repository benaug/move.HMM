\name{HSMM.plot}
\alias{HSMM.plot}
\title{Plot a move.HSMM object}
\usage{
  HSMM.plot(move.HSMM, xlim, breaks, by)
}
\arguments{
  \item{move.HSMM}{A move.HSMM object containing a fitted
  HSMM model.}

  \item{xlim}{An optional nstate x 2 matrix containing the
  minimum and maximum x values to be used for plotting each
  distribution.  If no matrix is supplied, the range of the
  plot is the range of the data.  xlim must be specified
  for correct plotting of distributions with support on the
  integers (e.g. poisson; see example).}

  \item{breaks}{An optional vector of length nstate
  containing the number of breaks to use for the histograms
  for each distribution.  It also accepts the breaks
  methods of the hist function.}

  \item{by}{An optional vector of length nstate containing
  the spacing of x points for each distribution.  Default
  is 0.001. by must be specified for correct plotting of
  distributions with support on the integers (e.g. poisson;
  see example).}

  \item{turn}{Parameters determining the transformation for
  circular distributions. turn=1 leads to support on
  (0,2pi) and turn=2 leads to support on (-pi,pi).  For
  animal movement models, the "encamped" state should use
  turn=1 and the "traveling" state should use turn=2}
}
\value{
  plots as described above
}
\description{
  This function plots the fitted distributions against
  histograms of the data for each observation variable.
  Fitted distributions are scaled by the stationary
  distribution of the transition matrix.  Scaled
  distributions are added and this combined distribution is
  plotted in red.  These plots may be helpful in assessing
  lack of fit or convergence problems.  The x scale and
  number of bins may need to be changed for useful plots.
}

