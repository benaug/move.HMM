\name{move.HSMM.CI}
\alias{move.HSMM.CI}
\title{Add confidence intervals to a move.HSMM object}
\usage{
  move.HSMM.CI(move.HSMM, CI = "boot", alpha = 0.05,
    B = 100, cores = 2, ..., useRcpp = FALSE)
}
\arguments{
  \item{move.HSMM}{A fitted move.HSMM object.}

  \item{CI}{A character determining which type of CI is to
  be calculated.  Current options are "FD" for the finite
  differences Hessian and "boot" for parametric
  bootstrapping and percentile CIs.}

  \item{alpha}{Type I error rate for CIs.  alpha=0.05 for
  95 percent CIs}

  \item{B}{Number of bootstrap resamples}

  \item{cores}{Number of cores to be used in parallel
  bootstrapping}

  \item{\dots}{additional arguments to
  \code{\link{move.HSMM.mle}}}

  \item{useRcpp}{Logical indicating whether or not to use
  Rcpp.}
}
\value{
  A list containing the lower and upper confidence bounds
}
\description{
  move.HSMM.CI is used to add confidence intervals to a
  move.HSMM object.  Current options are parametric
  bootstrap percentile CIs or CIs calculated from the
  finite differences Hessian. Bootstrapping utilizes
  parallel processing on a local host.  For each bootstrap
  sample, NAs in the original data are inserted in the
  simulated data in the same positions.  CI's based on the
  finite differences Hessian should generally be treated
  skeptically unless you've done a coverage analysis for
  your model with the amount of data you have.  Bootstrap
  samples are stored in move.HSMM$store boot so they can be
  inspected or combined with more bootstrap samples.
}

