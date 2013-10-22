#'Shifted lognormal quantile function
#'
#'This function evaluates the shifted lognormal quantile function.  The shift must be 
#'positive.
#'
#'@param x a vector of values where the quantile function is to be evaluated
#'@param mlog mean of the distribution on the log scale
#'@param sdlog standard deviation of the distribution on the log scale
#'@param x0 shift parameter (must be positive)
#'@return A vector of shifted lognormal quantile function values
#'@export

qlnorm3 <- function(x, mlog, sdlog, x0=0) {
  d=qlnorm(x,mlog,sdlog)+x0
  return(d)
}