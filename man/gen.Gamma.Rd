\name{gen.Gamma}
\alias{gen.Gamma}
\title{Generate t.p.m. to approximate HSMM with HMM}
\usage{
  gen.Gamma(m, params, PDFs, CDFs)
}
\arguments{
  \item{m}{vector of length nstates indicating the number
  of states to be in each state aggregate (see Langrock and
  Zuchinni 2011).}

  \item{params}{A list of length ndist+1 containing
  matrices of starting parameter values.  If nstates=2, the
  first element of the list must be the starting values for
  the dwell time distribution.  If nstates>2, the first
  element must be the starting values for t.p.m., which
  must have 0's on the diagonal and rowSums=1.}

  \item{PDFs}{A list of PDFs for the ndist distributions.}

  \item{CDFs}{A list of CDFs for the ndist distributions.}
}
\value{
  A vector of shifted negative binomial pdf values
}
\description{
  This function generates a t.p.m. of a HMM that
  approximates a HSMM (see Langrock and Zucchini, 2011)
}

