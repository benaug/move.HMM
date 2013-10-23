#'Compute and plot pseudoresiduals
#'
#'This function, modified from Zuchinni and MacDonald (2009), plots the ordinary
#'normal pseudoresiduals.  Three plots are produced-a qq-plot, normal pseudoresiduals
#'plotted through time, and an ACF plot.  Pseudoresiduals for integer valued
#'distributions need to be fixed.
#'
#'@param move.HSMM A move.HSMM object containing a fitted HSMM model.
#'@param plots A logical indicating whether or not to produce plots
#'@param returnresids A logical indicating whether or not to return the pseudo-residuals
#'@include Distributions.R
#'@include move.HSMM.state_probs.R
#'@return A list of normal pseudoresiduals of length ndist.
#'@export
#'
move.HSMM.psresid=function(move.HSMM,plots=T,returnresids=F){ 
  params=move.HSMM$params
  nstates=move.HSMM$nstates
  dists=move.HSMM$dists
  obs=move.HSMM$obs
  m=move.HSMM$m1
  ndist=length(dists)-1
  out=Distributions(dists,nstates)
  PDFs=out[[3]]
  CDFs=out[[4]]
  Gamma <- gen.Gamma(m,params,PDFs,CDFs)
  delta <- solve(t(diag(sum(m))-Gamma+1),rep(1,sum(m)))
  delta[delta<0]=1e19
  n <- nrow(obs)
  fb <- move.HSMM.lalphabeta(move.HSMM)
  if(nstates>2){
    params[[1]]=NULL
  }
  la <- t(fb$la)
  lb <- t(fb$lb)
  la <- cbind(log(delta),la)
  lafact <- apply(la ,2,max)
  lbfact <- apply(lb ,2,max)
  w <- matrix(NA ,ncol=n,nrow=sum(m))
  for (i in 1:n){
    foo <- (exp(la[,i]-lafact[i])%*% Gamma)*exp(lb[,i]-lbfact[i])
    w[,i] <- foo/sum(foo)
  }
  resids=vector("list",ndist-1)
  for(i in 1:ndist){
    resids[[i]]=numeric(n)
  }
  mstart=c(1,cumsum(m)+1)
  mstart=mstart[-length(mstart)]
  mstop=cumsum(m)
  for(k in 1:n){
    #for each distribution
    for(j in 1:nstates){
      for(i in 1:(ndist)){
        nparam=max(1,ncol(params[[i+1]]))
        if(nparam==2){
          resids[[i]][k]=resids[[i]][k]+sum(CDFs[[i+1]](obs[k,i],params[[i+1]][j,1],params[[i+1]][j,2])*w[mstart[j]:mstop[j],k])
        }else if(nparam==1){
          resids[[i]][k]=resids[[i]][k]+sum(CDFs[[i+1]](obs[k,i],params[[i+1]][j,1])*w[mstart[j]:mstop[j],k])
        }else if(nparam==3){
          resids[[i]][k]=resids[[i]][k]+sum(CDFs[[i+1]](obs[k,i],params[[i+1]][j,1],params[[i+1]][j,2],params[[i+1]][j,3])*w[mstart[j]:mstop[j],k])
        }
      }
    }
  }
  #remove NA resids
  for(i in 1:ndist){
    rem=which(is.na(obs[,i]))
    if(length(rem)>0){
      resids[[i]]=resids[[i]][-rem]
    }
  }
  #Plots
  if(plots==T){
    label=names(CDFs)[2:length(CDFs)]
    par(mfrow=c(ndist,1))
    #QQ plots
    for(i in 1:ndist){
      if(i==2)par(ask = TRUE)
      qqnorm(qnorm(resids[[i]]),main=paste(label[i],"Q-Q plot"),xlab="",ylab="",xlim=c(-3,3),ylim=c(-3,3))
      abline(a=0,b=1)
    }
    par(ask = TRUE)
    #resids through time
    for(i in 1:ndist){
      plot(qnorm(resids[[i]]),pch=1,main=paste(label[i],"Residuals through time"), ylab="residual")
      abline(h=0)
      lines(loess.smooth(seq(1:length(resids[[i]])),qnorm(resids[[i]])),col="red",lwd=2)
    }
    #acf
    for(i in 1:ndist){
      a=acf(qnorm(resids[[i]]),na.action=na.pass,plot=F)
      plot(a,pch=1,main=paste(label[i],"ACF"))
    }
    par(mfrow=c(1,1))
    par(ask = F)
  }
  if(returnresids==T){
    names(resids)=dists
    resids  
  }
}