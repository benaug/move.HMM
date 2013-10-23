#'Generate random deviates shifted Poisson pdf
#' 
#'This function generates random deviates from the shifted Poisson distribution.
#'The shift is fixed to 1.
#'@param n the number of random deviates to generate
#'@param lambda vector of positive means (of an ordinary Poisson distribution)
#'@return A vector of shifted Poisson random deviates
#'@export
rshiftpois=function(n,lambda){
  out=rpois(n,lambda)+1
  return(out)
}