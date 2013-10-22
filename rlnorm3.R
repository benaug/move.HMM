#'Shifted lognormal random number generation
#'
#'This function generates random numbers from the shifted lognormal distribution.
#'The shift must be positive.
#'
#'@param n the number of random numbers to generate
#'@param mlog mean of the distribution on the log scale
#'@param sdlog standard deviation of the distribution on the log scale
#'@param x0 shift parameter (must be positive)
#'@return A vector random numbers from a shifted lognormal distribution
#'@export
rlnorm3 <- function(n, mlog, sdlog, x0=0) {
out=qlnorm3(runif(n),mlog,sdlog,x0)
out
}