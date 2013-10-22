#'Generate random deviates shifted negative binomial pdf
#' 
#'This function generates random deviates from the shifted negative binomial distribution.
#'The shift is fixed to 1.
#'@param n the number of random deviates to generate
#'@param size target for number of successful trials, or dispersion parameter (the shape parameter of the gamma mixing distribution). Must be strictly positive, need not be integer.
#'@param prob probability of success in each trial. 0 < prob <= 1.
#'@return A vector of shifted negative binomial random deviates
#'@export
rshiftnegbin=function(n,size,prob){
  out=rnbinom(n,size,prob)+1
  return(out)
}