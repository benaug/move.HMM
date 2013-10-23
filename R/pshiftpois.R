#'Shifted Poisson cdf
#'
#'This function evaluates the shifted Poisson cdf.  The shift is fixed to 1.
#'@param x a vector of values where the cdf is to be evaluated
#'@param lambda vector of positive means (of an ordinary Poisson distribution)
#'@return A vector of shifted negative binomial cdf values
#'@export
pshiftpois=function(x,lambda){
  if(any(x<1))stop("This distribution only accepts values >=1")
  out=ppois(x-1,lambda)
  return(out)
}