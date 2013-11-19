#'Zero-altered Negative Binomial random number generation
#'
#'This function generates random numbers from the zero-altered negative binomial pdf.  This function is modified
#'from the rzanegbin VGAM package so that munb is not a parameter input option. 
#'
#'@param x a vector of values where the pdf is to be evaluated
#'@param size a value for the ordinary Negative Binomial size parameter
#'@param prob  a value for the ordinary Negative Binomial prob parameter
#'@param pobs0 Probability of zero.
#'@return A vector of zero-altered negative binomial pdf values
#'@export
rzanegbin=function (n, size, prob = NULL, pobs0 = 0) 
{
  use.n <- if ((length.n <- length(n)) > 1) 
    length.n
  else if (!is.Numeric(n, integer.valued = TRUE, allowable.length = 1, 
                       positive = TRUE)) 
    stop("bad input for argument 'n'")
  else n
  ans <- rposnegbin(n = use.n, prob = prob, size = size)
  if (length(pobs0) != use.n) 
    pobs0 <- rep(pobs0, len = use.n)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1)) 
    stop("argument 'pobs0' must be between 0 and 1 inclusive")
  ifelse(runif(use.n) < pobs0, 0, ans)
}