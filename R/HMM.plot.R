
#'Plot a move.HMM object
#'
#'This function plots the fitted distributions against histograms of the data
#'for each observation variable.  Fitted distributions are scaled by the 
#'stationary distribution of the transition matrix.  Scaled distributions are
#'added and this combined distribution is plotted in red.  These plots may be
#'helpful in assessing lack of fit or convergence problems.  The x scale and
#'number of bins may need to be changed for useful plots.
#'
#'@param move.HMM A move.HMM object containing a fitted HMM model.
#'@param xlim An optional nstate x 2 matrix containing the minimum and maximum x values to be used
#'for plotting each distribution.  If no matrix is supplied, the range of the 
#'plot is the range of the data.  xlim must be specified for correct
#'plotting of distributions with support on the integers (e.g. poisson; see example).
#'@param breaks An optional  vector of length nstate containing the number
#'of breaks to use for the histograms for each distribution.  It also accepts
#'the breaks methods of the hist function.
#'@param by An optional vector of length nstate containing the spacing of x
#'points for each distribution.  Default is 0.001. by must be specified for correct
#'plotting of distributions with support on the integers (e.g. poisson; see example).
#'@return plots as described above
#'@include Distributions.R
#'@export

#Plot fitted distributions agains histograms of data.  Distributions weighted by
#stationary distrubution
HMM.plot=function(move.HMM,xlim,breaks,by=NULL){
  obs=move.HMM$obs
  dists=move.HMM$dists
  params=move.HMM$params
  nstates=move.HMM$nstates
  out=Distributions(dists,nstates)
  PDFs=out[[3]]
  ndist=length(dists)
  n=nrow(obs)
  discrete=which(dists%in%c("poisson","nbinom","geometric","logarithmic","binom","pospoisson","posnegbin","zapois","zanegbin"))
  if(is.null(by)){
    by=rep(0.001,ndist)
    if(length(discrete>0)){
      by[discrete]=rep(1,length(discrete))
    }
  }
  if(length(discrete>0)){
    for(i in length(discrete)){
      if(!all(floor(xlim[discrete[i],])==xlim[discrete[i],]))stop("Make sure discrete RVs have discrete support in 'xlim'")
      if(!all(floor(by[discrete[i]])==by[discrete[i]]))stop("Make sure discrete RVs are incremented by integers in 'by'")
    }
  }
  
  x=vector("list",ndist)  
  y=vector("list",ndist)
  if(missing(xlim)){
    mins=apply(obs,2,min,na.rm=T)
    maxs=apply(obs,2,max,na.rm=T)
    xlim=cbind(mins,maxs)
  }
  if(missing(breaks)){
    breaks=rep("sturges",ndist)
  }
  #get PDF values for each distribution
  par(mfrow=c(ndist,1))
  par(mfrow=c(1,1))
  for(i in 1:ndist){
    x[[i]]=seq(xlim[i,1],xlim[i,2],by[i])
    y[[i]]=matrix(NA,nrow=length(x[[i]]),ncol=nstates+1)
    for(j in 1:nstates){
      nparam=max(1,ncol(params[[i+1]]))
      if(nparam==2){
        y[[i]][,j]=PDFs[[i]](x[[i]],params[[i+1]][j,1],params[[i+1]][j,2])
      }else if(nparam==1){
        y[[i]][,j]=PDFs[[i]](x[[i]],params[[i+1]][j,1])          
      }else if(nparam==3){
        y[[i]][,j]=PDFs[[i]](x[[i]],params[[i+1]][j,1],params[[i+1]][j,2],params[[i+1]][j,3])          
      }
    }
    #plots
    deltarep=matrix(rep(move.HMM$delta,length(x[[i]])),byrow=T,ncol=nstates)
    y[[i]][,1:nstates]=deltarep*y[[i]][,1:nstates]
    if(nstates>1){
      y[[i]][,nstates+1]=rowSums(y[[i]][,1:nstates])
    }else{
      y[[i]][,nstates+1]=y[[i]][,1]
    }
    y2=unlist(y[[i]])
    if(!all(is.finite(y2)))stop("x value supplied not in support of distribution")
    a=hist(obs[,i],plot=F,breaks=breaks[i])
    if(max(c(y2,max(a$density)))<2*max(a$density)){
      ymax=max(c(y2,max(a$density)))*1.1      
    }else{
      ymax=max(a$density)*1.2
    }
    par(ask = TRUE)
    hist(obs[,i],freq=F,breaks=breaks[i],xlim=xlim[i,],ylim=c(0,ymax),main=dists[i],xlab="x")
    col=c(rep("black",nstates),"red")
    for(j in 1:(nstates+1)){
      lines(x[[i]],y[[i]][,j],col=col[j])  
    }
  }
  par(mfrow=c(1,1))
  par(ask = FALSE)
}