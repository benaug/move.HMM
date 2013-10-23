#'Compute and plot pseudoresiduals
#'
#'This function, modified from Zuchinni and MacDonald (2009), plots the ordinary
#'normal pseudoresiduals.  Three plots are produced-a qq-plot, normal pseudoresiduals
#'plotted through time, and an ACF plot.  Pseudoresiduals for integer valued
#'distributions need to be fixed.
#'
#'@param move.HMM A move.HMM object containing a fitted HMM model.
#'@param plots A logical indicating whether or not to produce plots
#'@param returnresids A logical indicating whether or not to return the pseudo-residuals
#'@include Distributions.R
#'@include move.HMM.state_probs.R
#'@return A list of normal pseudoresiduals of length ndist.
#'@export
#'
move.HMM.psresid=function(move.HMM,plots=T,returnresids=F){
  params <- move.HMM$params
  nstates <- move.HMM$nstates
  dists <- move.HMM$dists
  obs <- move.HMM$obs
  delta <- move.HMM$delta
  ndist <- length(dists)
  out <- Distributions(dists,nstates)
  CDFs <- out[[4]]
  n <- nrow(obs)
  fb <- move.HMM.lalphabeta(move.HMM)
  la <- t(fb$la)
  lb <- t(fb$lb)
  la <- cbind(log(delta),la)
  lafact <- apply(la ,2,max)
  lbfact <- apply(lb ,2,max)
  w <- matrix(NA ,ncol=n,nrow=nstates)
  for (i in 1:n){
    foo <- (exp(la[,i]-lafact[i])%*% params[[1]])*exp(lb[,i]-lbfact[i])
    w[,i] <- foo/sum(foo)
  }
  resids <- vector("list",ndist)
  for(i in 1:ndist){
    resids[[i]] <- numeric(n)
  }
  for(k in 1:n){
    #for each distribution 
    for(j in 1:nstates){
      for(i in 1:ndist){
        nparam=max(1,ncol(params[[i+1]]))
        if(nparam==2){
          resids[[i]][k] <- resids[[i]][k]+CDFs[[i]](obs[k,i],params[[i+1]][j,1],params[[i+1]][j,2])*w[j,k]
        }else if(nparam==1){
          resids[[i]][k] <- resids[[i]][k]+CDFs[[i]](obs[k,i],params[[i+1]][j])*w[j,k]
        }else if(nparam==3){
          resids[[i]][k] <- resids[[i]][k]+CDFs[[i]](obs[k,i],params[[i+1]][j],params[[i+1]][j,2],params[[i+1]][j,3])*w[j,k]
        }
      }
    }
  }
  #remove NA resids (from NA observations)
  for(i in 1:ndist){
    rem <- which(is.na(obs[,i]))
    if(length(rem)>0){
      resids[[i]] <- resids[[i]][-rem]
    }
  }
  #Plots
  if(plots==T){
    label <- names(CDFs)
    par(mfrow=c(ndist,1))
    #QQ plots
    par(ask = TRUE)
    for(i in 1:ndist){
      if(i==2)par(ask = TRUE)
      qqnorm(qnorm(resids[[i]]),main=paste(label[i],"Q-Q plot"),xlab="",ylab="",xlim=c(-3,3),ylim=c(-3,3))
      abline(a=0,b=1)
    }
    #resids through time
    for(i in 1:ndist){
      plot(qnorm(resids[[i]]),pch=1,main=paste(label[i],"Residuals through time"), ylab="residual")
      abline(h=0)
      lines(loess.smooth(seq(1:length(resids[[i]])),qnorm(resids[[i]])),col="red",lwd=2)
    }
    #acf
    for(i in 1:ndist){
      a <- acf(qnorm(resids[[i]]),na.action=na.pass,plot=F)
      plot(a,pch=1,main=paste(label[i],"ACF"), ylab="residual")
    }
    par(mfrow=c(1,1))
    par(ask = F)
  }
  if(returnresids==T){
    names(resids) <- dists
    resids  
  }
}