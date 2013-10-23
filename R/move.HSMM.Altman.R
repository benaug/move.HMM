#'Altman plot
#'
#'This function provides a plot from Altman (2004) that allows the assessment
#'of goodness-of-fit.  Assuming the observed process is stationary, it plots
#'the marginal CDF of the data evaluated at the maximum likelihood estimates
#'against the empirical CDF values.  The plot should be close to a 45 degree straight
#'line.
#'
#'@param move.HSMM A move.HSMM object containing a fitted HSMM model.
#'@include Distributions.R
#'@include move.HSMM.state_probs.R
#'@include gen.Gamma.R
#'@return A plot
#'@export
#'
move.HSMM.Altman=function(move.HSMM){
  obs=move.HSMM$obs
  params=move.HSMM$params
  nstates=move.HSMM$nstates
  dists=move.HSMM$dists
  m=move.HSMM$m1
  ndist=length(dists)-1
  out=Distributions(dists,nstates)
  PDFs=out[[3]]
  CDFs=out[[4]]
  Gamma <- gen.Gamma(m,params,PDFs,CDFs)
  if(nstates>2){
    params[[1]]=NULL
  }
  delta <- solve(t(diag(sum(m))-Gamma+1),rep(1,sum(m)))
  n=nrow(obs)
  Ffit=vector("list",2)
  Femp=vector("list",2)
  for(i in 1:ndist){
    Ffit[[i]]=numeric(n)
  }
  circ=c("wrpcauchy","wrpnorm")
  for(k in 2:length(dists)){
    if(match(dists[k],circ,nomatch=0)>0){
      if(any(obs[,k-1]>pi,na.rm=T)){
        obs[,k-1][(obs[,k-1]>pi)&(!is.na(obs[,k-1]))]=obs[,k-1][(obs[,k-1]>pi)&(!is.na(obs[,k-1]))]-2*pi
      }
    }
  }
  #Calculate fitted marginal distribution
  mstart=c(1,cumsum(m)+1)
  mstart=mstart[-length(mstart)]
  mstop=cumsum(m)
  for(k in 1:n){
    #for each distribution
    for(j in 1:nstates){
      for(i in 1:ndist){
        nparam=max(1,ncol(params[[i+1]]))
        if(nparam==2){
          Ffit[[i]][k]=Ffit[[i]][k]+sum(CDFs[[i+1]](obs[k,i],params[[i+1]][j,1],params[[i+1]][j,2])*delta[mstart[j]:mstop[j]])
        }else if(nparam==1){
          Ffit[[i]][k]=Ffit[[i]][k]+sum(CDFs[[i+1]](obs[k,i],params[[i+1]][j])*delta[mstart[j]:mstop[j]])
        }else if(nparam==3){
          Ffit[[i]][k]=Ffit[[i]][k]+sum(CDFs[[i+1]](obs[k,i],params[[i+1]][j],params[[i+1]][j,2],params[[i+1]][j,3])*delta[mstart[j]:mstop[j]])
        }
      }
    }
  }
  label=names(CDFs)[2:length(CDFs)]
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