#'Shifted Poisson pdf
#'
#'This function evaluates the shifted Poisson pdf.  The shift is fixed to 1.
#'@param x a vector of values where the pdf is to be evaluated
#'@param lambda vector of positive means (of an ordinary Poisson distribution)
#'@return A vector of shifted negative binomial pdf values
#'@export
dshiftpois=function(x,lambda){
  if(any(x<1))stop("This distribution only accepts values >=1")
  out=dpois(x-1,lambda)
  return(out)
}