\name{move.HMM.lalphabeta}
\alias{move.HMM.lalphabeta}
\title{Compute forward and backwards probabilities}
\usage{
  move.HMM.lalphabeta(move.HMM)
}
\arguments{
  \item{move.HMM}{a move.HMM object containing a fitted HMM
  model.}
}
\value{
  A 2 element list containing the forward and backwards
  probabilities.
}
\description{
  This function, modified from Zucchini and MacDonald
  (2009), computes the forward and backwards probabilities
  defined by Equations (4.1) and (4.2) on page 60 in
  Zucchini and MacDonald (2009).  It takes as input a
  move.HMM object.  This function is used to create
  pseudo-residuals and in local decoding.
}

