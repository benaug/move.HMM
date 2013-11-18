#'Compute forward and backwards probabilities for a move.HSMM object
#'
#'This function, modified from Zucchini and MacDonald (2009), computes the forward
#'and backwards probabilities defined by Equations (4.1) and (4.2) on page 60
#'in Zucchini and MacDonald (2009).  It takes as input a move.HSMM object.
#'
#'@param move.HSMM a move.HSMM object containing a fitted HSMM model.
#'@include Distributions.R
#'@return A 2 element list containing the forward and backwards probabilities.
#'@export
move.HSMM.lalphabeta <-function(move.HSMM){
  obs <- move.HSMM$obs
  params <- move.HSMM$params
  nstates <- move.HSMM$nstates
  dists <- move.HSMM$dists
  m1 <- move.HSMM$m1
  out <- Distributions(dists,nstates)
  PDFs <- out[[3]]
  CDFs <- out[[4]]
  n <- dim(obs)[1]
  sm1=sum(m1)
  lalpha=lbeta=pre=allprobs<- matrix(rep(1,sm1*n),nrow=n)
  Gamma <- gen.Gamma(m1,params,PDFs,CDFs)
  delta <- solve(t(diag(sm1)-Gamma+1),rep(1,sm1))
  delta[delta<0]=1e-19
  mstart <- c(1,cumsum(m1)+1)
  mstart <- mstart[-length(mstart)]
  mstop <- cumsum(m1)
  if(nstates>2){
    params[[1]]=NULL
  }
  ndists=length(PDFs)
  nparam=unlist(lapply(params,ncol))
  use=!is.na(obs)*1
  for(i in 2:ndists){
    if(nparam[i]==2){
      #for 2 parameter distributions
      for (j in 1:nstates){
        allprobs[use[,i-1],mstart[j]:mstop[j]] <- allprobs[use[,i-1],mstart[j]:mstop[j]]*matrix(rep(PDFs[[i]](obs[use[,i-1],i-1],params[[i]][j,1],params[[i]][j,2]),m1[j]),ncol=m1[j])
      }
    }else if(nparam[i]==1){
      #for 1 parameter distributions. 
      for (j in 1:nstates){
        allprobs[use[,i-1],mstart[j]:mstop[j]] <- allprobs[use[,i-1],mstart[j]:mstop[j]]*matrix(rep(PDFs[[i]](obs[use[,i-1],i-1],params[[i]][j]),m1[j]),ncol=m1[j])
      }
    }else if(nparam[i]==3){
      #for 3 parameter distributions
      for (j in 1:nstates){
        allprobs[use[,i-1],mstart[j]:mstop[j]] <- allprobs[use[,i-1],mstart[j]:mstop[j]]*matrix(rep(PDFs[[i]](obs[use[,i-1],i-1],params[[i]][j,1],params[[i]][j,2],params[[i]][j,3]),m1[j]),ncol=m1[j])
      }
    }
  }
  foo <- delta*allprobs [1,]
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo/sumfoo
  pre[1,] <- foo
  lalpha [1,] <- log(foo)+lscale
  for (i in 2:n){
    foo <- foo %*% Gamma *allprobs[i,]
    sumfoo <- sum(foo)
    lscale <- lscale+log(sumfoo)
    foo <- foo/sumfoo
    pre[i,] <- foo
    lalpha[i,] <- log(foo) +lscale
  }
  lbeta[n,] <- rep(0,nstates)
  foo <- rep (1/nstates,nstates)
  lscale <- log(nstates)
  for (i in (n-1) :1){
    foo <- Gamma  %*%( allprobs[i+1,]*foo)
    lbeta[i,] <- log(foo) +lscale
    sumfoo <- sum(foo)
    foo <- foo/sumfoo
    lscale <- lscale+log(sumfoo)
  }
  list(la=lalpha ,lb=lbeta,pre=pre)
}