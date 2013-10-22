#'ACF plot
#'
#'This function compares the empirical ACFs to those simulated from the fitted HSMM.  Experimental.
#'@param move.HSMM A move.HMM object containing a fitted HSMM model.
#'@param simlength The number of observations to simulate.  The ACF from the simulated data will 
#'converge to the theoretical ACF as simlength goes to infinity
#'@param transforms A list of length ndist that contains functions for transforming the data for each distribution.
#'Default is NULL, so data are not transformed. 
#'@param lag.max Maximum lag at which to calculate the acf.  Default is 10.
#'@param ylim a ndist x 2 matrix with the lower and upper bounds for plotting each ACF.  Defaults to (-0.3,0.3).
#'@param tol numeric value indicating the distance between the empirical and simulated ACFs plots at each lag length. 
#'Defaults to 0.1.
#'@return A vector of shifted negative binomial pdf values
#'@include move.HMM.simulate
#'@export
move.HSMM.ACF=function(move.HSMM,simlength=10000,transforms=NULL,lag.max=10,ylim=NULL,tol=0.1){
  dists <- move.HSMM$dists
  params <- move.HSMM$params
  obs <- move.HSMM$obs
  sim <- move.HSMM.simulate(dists,params,simlength,move.HSMM$nstates)$obs
  ndist <- length(dists)-1
  plotat <- 1:lag.max
  fixy <- is.null(ylim)
  if(is.null(ylim))ylim=matrix(NA,nrow=ndist,ncol=2)
  for(i in 1:ndist){
    if(!is.null(transforms)){
      acf1 <- acf(transforms[[i]](obs[,i]),plot=F,lag.max=lag.max,na.action=na.pass)
      acf2 <- acf(transforms[[i]](sim[,i]),plot=F,na.action=na.pass)
    }else{
      acf1 <- acf(obs[,i],plot=F,lag.max=lag.max,na.action=na.pass)
      acf2 <- acf(sim[,i],plot=F,na.action=na.pass)
    }
    if(fixy){
      ylim[i,] <- c(-0.3,0.3)
      maximum <- max(c(acf1$acf[2:(lag.max+1)]),acf2$acf[2:(lag.max+1)])
      minimum <- min(c(acf1$acf[2:(lag.max+1)]),acf2$acf[2:(lag.max+1)])      
      if(maximum>0.3)ylim[i,2] <- min(maximum*1.2,1)
      if(minimum<=-0.3)ylim[i,1] <- max(minimum*1.2,-1)
    }
    plot(NA,xlim=c(1,lag.max),ylim=ylim[i,],xlab="Lag",ylab="ACF",main=paste(dists[i+1],"Observed vs. Simulated ACFs"))
    abline(c(0,0))
    for(i in 1:length(plotat)){
      lines(x=c(plotat[i]-tol,plotat[[i]]-tol),y=c(0,acf1$acf[i+1]))
    }
    for(i in 1:length(plotat)){
      lines(x=c(plotat[i]+tol,plotat[[i]]+tol),y=c(0,acf2$acf[i+1]),col="red")
    }
    legend(x="bottomright",lwd=c(1,1),col=c("black","red"),cex=0.7,legend=c("Observed","Simulated"))
    par(ask=T)
  }
  par(ask=F)
}