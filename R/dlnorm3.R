#'Shifted lognormal pdf
#'
#'This function evaluates the shifted lognormal pdf.  The shift must be 
#'positive.
#'
#'@param x a vector of values where the pdf is to be evaluated
#'@param mlog mean of the distribution on the log scale
#'@param sdlog standard deviation of the distribution on the log scale
#'@param x0 shift parameter (must be positive)
#'@return A vector of shifted lognormal pdf values
#'@export
#'
dlnorm3 <- function(x, mlog, sdlog, x0=0) {
  d=rep(0,length(x))
  eval <- which(x>=x0)
  x <- x-x0
  d[eval] <- dlnorm(x[eval],mlog,sdlog)
  return(d)
}