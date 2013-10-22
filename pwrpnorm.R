#'Wrapped normal cdf
#'
#'This function evaluates the wrapped normal cdf.  Observations must be in
#'the interval [-pi,pi] or [0,2pi].
#'
#'@param q a vector of values where the cdf is to be evaluated
#'@param mu A value for the location parameter
#'@param rho A value for the concentration parameter in the interval [0,1]
#'@return A vector of wrapped normal cdf values
#'@export
pwrpnorm=function(q,mu=0,rho=0){
  Nin=length(q)
  out=rep(0,Nin)
  #fix any obs >pi
  if(any(q>pi,na.rm=T)){
    q[(q>pi)&(!is.na(q))]=q[(q>pi)&(!is.na(q))]-2*pi
  }
  zero=mu-pi
  adjustment=integrate(dwrpnorm,lower=-pi,upper=zero,mu=mu,rho=rho)$value
  #integrate from -pi to 0 then 0 to pi
  below=which(q<zero)
    for(i in 1:Nin){
      if(!is.na(q[i])){
        if(i %in% below){
          out[i]=abs(integrate(dwrpnorm,lower=-pi,upper=q[i],mu=mu,rho=rho)$value)
        }else{
          out[i]=adjustment+integrate(dwrpnorm,lower=zero,upper=q[i],mu=mu,rho=rho)$value 
        }
      }
    }
  
  
  return(out)
}