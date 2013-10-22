#'Dwell Plot
#'
#'This function plots the (geometric) pdf of the fitted model against the dwell times from the state sequence
#'predicted by the Viterbi algorithm.
#'
#'@param move.HMM A move.HMM object containing a fitted HMM model.
#'@include move.HMM.viterbi
#'@return A plot of dwell times
#'@export
#'
move.HMM.dwell.plot=function(move.HMM){
  if(move.HMM$nstates==1)stop("Function does not work if nstates=1")
  pred=move.HMM.viterbi(move.HMM)
  dwell=rle(pred)
  states=1:move.HMM$nstates
  stateDwell=vector("list",move.HMM$nstates)
  p=1-diag(move.HMM$params$tmat)
  par(mfrow=c(move.HMM$nstates,1))
  for(i in 1:move.HMM$nstates){
    stateDwell[[i]]=dwell$lengths[dwell$values==i]
    neach=table(stateDwell[[i]])
    total=length(stateDwell[[i]])
    density=neach/total
    plotat1=as.numeric(names(density))
    plotat2=0:(max(as.numeric(names(density)))+3)
    density2=dgeom(plotat2,p[i])
    plotat2=plotat2+1
    plot(density2~plotat2,ylim=c(0,max(density,density2)),xlab="Dwell Time",ylab="Density",main=paste("Predicted vs. Observed Dwell Times State",i))
    for(i in 1:length(plotat1)){
      lines(x=c(plotat1[i],plotat1[[i]]),y=c(0,density[i]))
    }
    legend(x="topright",pch=c(1,NA),legend=c("Predicted","Observed"))
    legend(x="topright",pch=c(NA,"l"),legend=c("Predicted","Observed"),bty="n")
  }
  par(mfrow=c(1,1))
}