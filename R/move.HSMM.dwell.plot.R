#'Dwell Plot
#'
#'This function plots the pdf of the fitted model against the dwell times from the state sequence
#'predicted by the Viterbi algorithm.
#'
#'@param move.HSMM A move.HSMM object containing a fitted HSMM model.
#'@include move.HSMM.viterbi.R
#'@return A plot of dwell times
#'@export
#'
move.HSMM.dwell.plot=function(move.HSMM){
  pred=move.HSMM.viterbi(move.HSMM)
  dwell=rle(pred)
  nstates=move.HSMM$nstates
  params=move.HSMM$params
  states=1:nstates
  stateDwell=vector("list",nstates)
  dists=move.HSMM$dists
  out=Distributions(dists,nstates)
  PDFs=out[[3]]
  par(mfrow=c(nstates,1))
  if(nstates>2){
    params[[1]]=NULL
  }
  for(i in 1:move.HSMM$nstates){
    stateDwell[[i]]=dwell$lengths[dwell$values==i]
    neach=table(stateDwell[[i]])
    total=length(stateDwell[[i]])
    density=neach/total
    plotat1=as.numeric(names(density))
    plotat2=1:(max(as.numeric(names(density)))+3)
    nparams=ncol(params[[1]])
    if(nparams==2){
      density2=PDFs[[1]](plotat2,params[[1]][i,1],params[[1]][i,2])
    }else if(nparams==1){
      density2=PDFs[[1]](plotat2,params[[1]][i,1])
    }else if(nparams==3){
      density2=PDFs[[1]](plotat2,params[[1]][i,1],params[[1]][i,2],params[[1]][i,3])
    }
    plot(density2~plotat2,ylim=c(0,max(density,density2)),xlab="Dwell Time",ylab="Density",main=paste("Predicted vs. Observed Dwell Times State",i))
    for(i in 1:length(plotat1)){
      lines(x=c(plotat1[i],plotat1[[i]]),y=c(0,density[i]))
    }
    legend(x="topright",pch=c(1,NA),legend=c("Predicted","Observed"),cex=0.7)
    legend(x="topright",pch=c(NA,"l"),legend=c("Predicted","Observed"),bty="n",cex=0.7)
  }
  par(mfrow=c(1,1))
}
