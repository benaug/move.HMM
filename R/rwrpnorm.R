#'Wrapped normal random number generation
#'
#'This function generates random numbers from the wrapped normal distribution.
#'It is modified from the functoin in the CircStats package so that it can
#'evaluate multiple parameter combinations in the same call.
#'
#'@param n The number of random numbers to generate
#'@param mu A value for the mu parameter
#'@param rho A value for the concentration parameter in the interval [0,1]
#'@return A vector of random numbers drawn from a wrapped normal distribution
#'@export
rwrpnorm=function(n, mu, rho) 
{
  if(length(mu)!=length(rho))stop("Dimension of mu and rho must be equal")
  if(length(mu)==1){
    mu=rep(mu,n)
  }
  if(length(rho)==1){
    rho=rep(rho,n)
  }

  result=rep(NA,n)
  case1=which(rho==0)
  case2=which(rho==1)
  case3=1:n
  if(length(c(case1,case2))>0){
    case3=case3[-c(case1,case2)]
  }
  if(length(case1)>0){
    result[case1]=runif(length(case1), 0, 2 * pi)
  }
  if(length(case2)>0){
    result[case2]= rep(mu, length(case2))
  }
  if(length(case3)>0){
    sd <- sqrt(-2 * log(rho[case3]))
    result[case3] <- rnorm(n, mu, sd)%%(2 * pi)%%(2 * pi)
  }
  #shift support to (-pi,pi)
  result[result>pi]=result[result>pi]-2*pi
  
  result 
}