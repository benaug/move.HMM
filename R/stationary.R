#'Calculate the stationary distribution from the t.p.m.
#'
#'This function, is used to get the numerical derivative required
#'to get the covariance matrix of the stationary distribution.
#'
#'@param gamvect A vector containing the MLEs of the t.p.m.
#'@return a vector containing the stationary distribution
stationary=function(gamvect){
  gamma=matrix(gamvect,nrow=sqrt(length(gamvect)))
  station=solve(t(diag(nrow(gamma))-gamma+1),rep(1,nrow(gamma)))
  return(station)
}