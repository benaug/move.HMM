\name{move.HMM-package}
\alias{move.HMM-package}
\alias{move.HMM}
\docType{package}
\title{Fit hidden Markov models to animal movement and sensor data via maximum likelihood}
\description{
move.HMM is used to fit hidden Markov models , allowing for multiple observation variables with different distributions. Models can be compared via AICc and fit can be assessed by plotting the fitted models, ordinary normal pseudoresiduals, and goodness-of-fit plots introduced by Altman (2004). Conditional state probabilities can be calculated and global decoding can be performed via the Viterbi algorithm. Support for Rcpp is inlcluded which can greatly speed up fitting HMMs to longer time series and HSMMs, generally.  The package was developed for inferring behavioral states from animal movement and sensor data, but is more widely applicable.
}
\details{
\tabular{ll}{
Package: \tab move.HMM\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2013-10-08\cr
License: \tab GPL(>=2)\cr

See move.HMM.mle and move.HSMM.mle for details about how to use the package.
}
}
\author{
Ben Augustine and Roland Langrock

Maintainer: Ben Augustine <ben.augustine@uky.edu>
}
\references{
Langrock, R., King, R., Matthiopoulos, J., Thomas, L., Fortin, D., & Morales, J. M. (2012). Flexible and practical modeling of animal telemetry data: hidden Markov models and extensions. Ecology, 93(11), 2336-2342.

MacKay Altman, R. (2004). Assessing the Goodness-of-Fit of Hidden Markov Models. Biometrics, 60(2), 444-450.

Zucchini, W., & MacDonald, I. L. (2009). Hidden Markov models for time series: an introduction using R. CRC Press.
}
