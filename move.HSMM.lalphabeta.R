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
  m <- move.HSMM$m1
  out <- Distributions(dists,nstates)
  PDFs <- out[[3]]
  CDFs <- out[[4]]
  n <- dim(obs)[1]
  lalpha=lbeta=pre <- matrix(rep(1,(sum(m))*n),nrow=n)
  Gamma <- gen.Gamma(m,params,PDFs,CDFs)
  delta <- solve(t(diag(sum(m))-Gamma+1),rep(1,sum(m)))
  delta[delta<0]=1e19
  if(nstates>2){
    params[[1]]==NULL
  }
  allprobs <- matrix(rep(1,(sum(m))*n),nrow=n)
  mstart <- c(1,cumsum(m)+1)
  mstart <- mstart[-length(mstart)]
  mstop <- cumsum(m)
  if(nstates>2){
    params[[1]]=NULL
  }
  for (k in 1:n){
    for (j in 1:nstates){
      for(i in 2:length(PDFs)){
        nparam=max(1,ncol(params[[i]]))
        if(nparam==2){
          #for 2 parameter distributions
          allprobs[k,mstart[j]:mstop[j]] <- allprobs[k,mstart[j]:mstop[j]]*rep(ifelse(is.na(obs[k,i-1]),1,PDFs[[i]](obs[k,i-1],params[[i]][j,1],params[[i]][j,2])),m[j])
        }else if(nparam==1){
          #for 1 parameter distributions. 
          allprobs[k,mstart[j]:mstop[j]] <- allprobs[k,mstart[j]:mstop[j]]*rep(ifelse(is.na(obs[k,i-1]),1,PDFs[[i]](obs[k,i-1],params[[i]][j])),m[j])
        }else if(nparam==3){
          #for 3 parameter distributions
          allprobs[k,mstart[j]:mstop[j]] <- allprobs[k,mstart[j]:mstop[j]]*rep(ifelse(is.na(obs[k,i-1]),1,PDFs[[i]](obs[k,i-1],params[[i]][j,1],params[[i]][j,2],params[[i]][j,3])),m[j])
        }
      } #i index
    } # j index
  } # k index
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