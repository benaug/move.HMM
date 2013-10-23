\name{move.HSMM.lalphabeta}
\alias{move.HSMM.lalphabeta}
\title{Compute forward and backwards probabilities for a move.HSMM object}
\usage{
  move.HSMM.lalphabeta(move.HSMM)
}
\arguments{
  \item{move.HSMM}{a move.HSMM object containing a fitted
  HSMM model.}
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
  move.HSMM object.
}

