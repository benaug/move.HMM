#'Shifted negative binomial pdf
#'
#'This function evaluates the shifted negative binomial pdf.  The shift is fixed to 1.
#'
#'@param x a vector of values where the pdf is to be evaluated
#'@param size target for number of successful trials, or dispersion parameter (the shape parameter of the gamma mixing distribution). Must be strictly positive, need not be integer.
#'@param prob probability of success in each trial. 0 < prob <= 1.
#'@return A vector of shifted negative binomial pdf values
#'@export
#'
dshiftnegbin=function(x,size,prob){
  if(any(x<1))stop("This distribution only accepts values >=1")
  out=dnbinom(x-1,size,prob)
  return(out)
}