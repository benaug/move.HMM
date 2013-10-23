\name{move.HSMM.ACF}
\alias{move.HSMM.ACF}
\title{ACF plot}
\usage{
  move.HSMM.ACF(move.HSMM, simlength = 10000,
    transforms = NULL, lag.max = 10, ylim = NULL,
    tol = 0.1)
}
\arguments{
  \item{move.HSMM}{A move.HMM object containing a fitted
  HSMM model.}

  \item{simlength}{The number of observations to simulate.
  The ACF from the simulated data will converge to the
  theoretical ACF as simlength goes to infinity}

  \item{transforms}{A list of length ndist that contains
  functions for transforming the data for each
  distribution. Default is NULL, so data are not
  transformed.}

  \item{lag.max}{Maximum lag at which to calculate the acf.
  Default is 10.}

  \item{ylim}{a ndist x 2 matrix with the lower and upper
  bounds for plotting each ACF.  Defaults to (-0.3,0.3).}

  \item{tol}{numeric value indicating the distance between
  the empirical and simulated ACFs plots at each lag
  length. Defaults to 0.1.}
}
\value{
  A vector of shifted negative binomial pdf values
}
\description{
  This function compares the empirical ACFs to those
  simulated from the fitted HSMM.  Experimental.
}

