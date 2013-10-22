#'Altman plot
#'
#'This function provides a plot from Altman (2004) that allows the assessment
#'of goodness-of-fit.  Assuming the observed process is stationary, it plots
#'the marginal CDF of the data evaluated at the maximum likelihood estimates
#'against the empirical CDF values.  The plotted values should be follow a 45 degree straight
#'line.
#'
#'@param move.HMM A move.HMM object containing a fitted HMM model.
#'@include Distributions.R
#'@include move.HMM.state_probs.R
#'@return Altman plot
#'@export
#'
move.HMM.Altman=function(move.HMM){
  obs=move.HMM$obs
  params=move.HMM$params
  nstates=move.HMM$nstates
  dists=move.HMM$dists
  delta=move.HMM$delta
  ndist=length(dists)
  out=Distributions(dists,nstates)
  CDFs=out[[4]]  
  n=nrow(obs)
  Ffit=vector("list",2)
  Femp=vector("list",2)
  for(i in 1:ndist){
    Ffit[[i]]=numeric(n)
  }
  circ=c("wrpcauchy","wrpnorm")
  for(k in 1:length(dists)){
    if(match(dists[k],circ,nomatch=0)>0){
      if(any(obs[,k]>pi,na.rm=T)){
        obs[,k][(obs[,k]>pi)&(!is.na(obs[,k]))]=obs[,k][(obs[,k]>pi)&(!is.na(obs[,k]))]-2*pi
      }
    }
  }
  #Calculate fitted marginal distribution
  for(k in 1:n){
    #for each distribution 
    for(j in 1:nstates){
      for(i in 1:ndist){
        nparam=max(1,ncol(params[[i+1]]))
        if(nparam==2){
          Ffit[[i]][k]=Ffit[[i]][k]+CDFs[[i]](obs[k,i],params[[i+1]][j,1],params[[i+1]][j,2])*delta[j]
        }else if(nparam==1){
          Ffit[[i]][k]=Ffit[[i]][k]+CDFs[[i]](obs[k,i],params[[i+1]][j])*delta[j]
        }else if(nparam==3){
          Ffit[[i]][k]=Ffit[[i]][k]+CDFs[[i]](obs[k,i],params[[i+1]][j],params[[i+1]][j,2],params[[i+1]][j,3])*delta[j]
        }
      }
    }
  }
  label=names(CDFs)
  par(mfrow=c(ndist,1))
  for(i in 1:ndist){
    Femp[[i]]=numeric(nrow(obs))
    if(out[[6]][i]==1){
      a=ecdf(obs[,i]) 
      Femp[[i]]=a(obs[,i])
    }else if(out[[6]][i]==2){
      a=ecdf(obs[,i]) 
      Femp[[i]]=a(obs[,i])
    }
    plot(Ffit[[i]],Femp[[i]],main=paste(label[i],"Altman Plot"),xlab="Fitted",ylab="Empirical",xlim=c(0,1),ylim=c(0,1))
    abline(a=0,b=1)
  }
}