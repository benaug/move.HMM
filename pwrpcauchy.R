#'Wrapped cauchy cdf
#'
#'This function evaluates the wrapped cauchy cdf.  Observations must be in
#'the interval [-pi,pi] or [0,2pi].
#'
#'@param q a vector of values where the cdf is to be evaluated
#'@param mu A value for the location parameter
#'@param rho A value for the concentration parameter in the interval [0,1]
#'@return A vector of wrapped cauchy cdf values
#'@export

pwrpcauchy=function(q,mu=0,rho=0){
  Nin=length(q)
  out=rep(NA,Nin)
#fix any obs >pi
  if(any(q>pi,na.rm=T)){
    q[(q>pi)&(!is.na(q))]=q[(q>pi)&(!is.na(q))]-2*pi
  }
  
  zero=mu-pi
  adjustment=integrate(dwrpcauchy,lower=-pi,upper=zero,mu=mu,rho=rho)$value
  #integrate from -pi to 0 then 0 to pi
  below=which(q<zero)
  for(i in 1:Nin){
    if(!is.na(q[i])){
      if(i %in% below){
        out[i]=abs(integrate(dwrpcauchy,lower=-pi,upper=q[i],mu=mu,rho=rho)$value)
      }else{
        out[i]=adjustment+integrate(dwrpcauchy,lower=zero,upper=q[i],mu=mu,rho=rho)$value 
      }
    }
  }
  return(out)
}
