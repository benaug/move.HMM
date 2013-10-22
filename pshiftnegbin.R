#'Shifted negative binomial cdf
#'
#'This function evaluates the shifted negative binomial cdf.  The shift is fixed to 1.
#'@param x a vector of values where the cdf is to be evaluated
#'@param size target for number of successful trials, or dispersion parameter (the shape parameter of the gamma mixing distribution). Must be strictly positive, need not be integer.
#'@param prob probability of success in each trial. 0 < prob <= 1.
#'@return A vector of shifted negative binomial cdf values
#'@export
pshiftnegbin=function(x,size,prob){
  if(any(x<1))stop("This distribution only accepts values >=1")
  out=pnbinom(x-1,size,prob)
  return(out)
}