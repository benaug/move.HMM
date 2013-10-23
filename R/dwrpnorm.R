#'Wrapped normal pdf
#'
#'This function evaluates the wrapped normal pdf.  This function is modified
#'from the wrapped normal pdf in the CircStats package so that it can 
#'evaluate more than one value at a time
#'
#'@param theta a vector of values where the pdf is to be evaluated
#'@param mu a value for the location parameter in the interval [-pi,pi]
#'@param rho a value for the concentration parameter in the interval [0,1]
#'@param acc parameter defining the accuracy of the estimation of the density Terms are added to the infinite summation that defines the density function until successive estimates are within acc of each other.
#'@param tol the same as acc
#'@return A vector of wrapped normal pdf values
#'@export
dwrpnorm=function(theta, mu, rho, acc = 1e-05, tol = acc) 
{
  if (rho < 0 | rho > 1) 
    stop("rho must be between 0 and 1")
  var <- -2 * log(rho)
  term <- function(theta, mu, var, k) {
    1/sqrt(var * 2 * pi) * exp(-((theta - mu + 2 * pi * k)^2)/(2 * var))
  }
  n=length(theta)
  out=rep(NA,n)
  for(i in 1:n){
    k <- 0
    Next <- term(theta[i], mu, var, k)
    delta <- 1
    while (delta > tol) {
      k <- k + 1
      Last <- Next
      Next <- Last + term(theta[i], mu, var, k) + term(theta[i],mu, var, -k)
      delta <- abs(Next - Last)
    }
  out[i]=Next
  }
  out
}