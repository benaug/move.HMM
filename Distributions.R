#'Parameter transformations and distribution functions
#'
#'This function contains supported PDFs, CDFs, and random deviate generating
#'functions along with appropriate parameter transformations and inverse
#'transformations needed for maximum likelihood estimation.  It takes a vector of
#'distributions and the number of hidden states as input and outputs the
#'appropriate functions for those distributions.  This function is for internal use.
#'Additional distributions can be added by modifying this function.
#'
#'@param dists A length ndist vector of distributions from the following list:
#'weibull, gamma, exponential, normal, lognormal, lnorm3, posnorm,
#'invgamma, rayleigh, f, ncf, dagum, frechet, beta, geometric, logarithmic, binom, poisson, nbinom,
#'zapois, wrpcauchy, wrpnorm,shiftpois,shiftnegbin
#'@param nstates Number of hidden states
#'@param turn Parameters determining the transformation for circular distributions.
#'turn=1 yields support on (0,2pi) and turn=2 yields support on (-pi,pi).  For
#'animal movement models, the "encamped" state should use turn=1 and the "traveling"
#'state should use turn=2.
#'@return A list containing transformations, inverse transformations, PDFs, CDFs,
#'and functions to generate random deviates for the input distributions.
#'@export
#'@include dshiftnegbin.R
#'
Distributions=function(dists,nstates,turn=NULL){

#List supported distributions
validDists=c("weibull","gamma","exponential","normal","lognormal","lnorm3","posnorm","invgamma","rayleigh","f","ncf","dagum","frechet","beta","geometric","logarithmic","binom","poisson","nbinom","zapois","pospois","posnegbin","shiftpois","shiftnegbin","posgeom","wrpcauchy","wrpnorm")
if(!all(dists %in% validDists))stop("At least one distribution is not currently supported or not spelled correctly.  Check spelling and capitalization against the list in Distributions help file.")

#Preallocate
transforms=vector("list",length(validDists))
inv.transforms=vector("list",length(validDists))
names(transforms)=validDists
names(inv.transforms)=validDists

#Assign appropriate transformations for all available dists
transforms$weibull=c(log,log)
transforms$gamma=c(log,function(x){log(1/x)})
transforms$exponential=c(log)
transforms$normal=c(identity,log)
transforms$lognormal=c(identity,log)
transforms$lnorm3=c(identity,log,log)
transforms$posnorm=c(identity,log)
transforms$invgamma=c(log,log)
transforms$rayleigh=c(log)
transforms$f=c(log,log)
transforms$ncf=c(log,log,log)
transforms$dagum=c(log,log,log)
transforms$frechet=c(log,log,log)
transforms$beta=c(log,log)
transforms$geometric=c(qlogis)
transforms$logarithmic=c(qlogis)
transforms$binom=c(log,qlogis)
transforms$poisson=c(log)
transforms$nbinom=c(log,qlogis)
transforms$zapois=c(log,log)
transforms$pospois=c(log)
transforms$posnegbin=c(log,qlogis)
transforms$shiftpois=c(log)
transforms$shiftnegbin=c(log,qlogis)
transforms$posgeom=c(qlogis)
transforms$wrpcauchy=vector("list",nstates+1)
transforms$wrpnorm=vector("list",nstates+1)
if(!is.null(turn)){
  for(i in 1:(nstates)){
    if(turn[i]==2){ #go straight
      transforms$wrpcauchy[[i]]=function(x){log((pi+x)/(pi-x))}
      transforms$wrpnorm[[i]]=function(x){log((pi+x)/(pi-x))}
    }else if(turn[i]==1){  #turn around
      transforms$wrpcauchy[[i]]=function(x){log((x)/(2*pi-x))}
      transforms$wrpnorm[[i]]=function(x){log((x)/(2*pi-x))}
    }
  }
}
transforms$wrpcauchy[[nstates+1]]=qlogis
transforms$wrpnorm[[nstates+1]]=qlogis

#Assign appropriate inverse transformations for all available dists
inv.transforms$weibull=c(exp,exp)
inv.transforms$gamma=c(exp,function(x){1/exp(x)})
inv.transforms$exponential=c(exp)
inv.transforms$normal=c(identity,exp)
inv.transforms$lognormal=c(identity,exp)
inv.transforms$lnorm3=c(identity,exp,exp)
inv.transforms$posnorm=c(identity,exp)
inv.transforms$invgamma=c(exp,exp)
inv.transforms$rayleigh=c(exp)
inv.transforms$f=c(exp,exp)
inv.transforms$ncf=c(exp,exp,exp)
inv.transforms$dagum=c(exp,exp,exp)
inv.transforms$frechet=c(exp,exp,exp)
inv.transforms$geometric=c(plogis)
inv.transforms$logarithmic=c(plogis)
inv.transforms$beta=c(exp,exp)
inv.transforms$binom=c(exp,plogis)
inv.transforms$poisson=c(exp)
inv.transforms$nbinom=c(exp,plogis)
inv.transforms$zapoisson=c(exp,exp)
inv.transforms$pospois=c(exp)
inv.transforms$posnegbin=c(exp,plogis)
inv.transforms$shiftpois=c(exp)
inv.transforms$shiftnegbin=c(exp,plogis)
inv.transforms$posgeom=c(plogis)
inv.transforms$wrpcauchy=vector("list",nstates+1)
inv.transforms$wrpnorm=vector("list",nstates+1)

if(!is.null(turn)){
  for(i in 1:(nstates)){
    if(turn[i]==1){
      inv.transforms$wrpcauchy[[i]]=function(x){2*pi*plogis(x)}
      inv.transforms$wrpnorm[[i]]=function(x){2*pi*plogis(x)}
    }else{
      inv.transforms$wrpcauchy[[i]]=function(x){pi*(exp(x)-1)/(exp(x)+1)}
      inv.transforms$wrpnorm[[i]]=function(x){pi*(exp(x)-1)/(exp(x)+1)}
    }
  }
}
inv.transforms$wrpcauchy[[nstates+1]]=plogis
inv.transforms$wrpnorm[[nstates+1]]=plogis

#Assign appropriate PDFs for all dists
PDFs=vector("list",length(validDists))
names(PDFs)=validDists
PDFs$weibull=dweibull
PDFs$gamma=dgamma
PDFs$exponential=dexp
PDFs$normal=dnorm
PDFs$lognormal=dlnorm
PDFs$lnorm3=dlnorm3
PDFs$posnorm=dposnorm
PDFs$invgamma=densigamma
PDFs$rayleigh=drayleigh
PDFs$f=df
PDFs$ncf=df
PDFs$dagum=ddagum
PDFs$frechet=dfrechet
PDFs$beta=dbeta
PDFs$geometric=dgeom
PDFs$logarithmic=dlog
PDFs$binom=dbinom
PDFs$poisson=dpois
PDFs$nbinom=dnbinom
PDFs$zapois=dzapois
PDFs$pospois=dpospois
PDFs$posnegbin=dposnegbin
PDFs$shiftpois=dshiftpois
PDFs$shiftnegbin=dshiftnegbin
PDFs$posgeom=dposgeom
PDFs$wrpcauchy=dwrpcauchy
PDFs$wrpnorm=dwrpnorm

#Assign appropriate CDFs for all dists
CDFs=vector("list",length(validDists))
names(CDFs)=validDists
CDFs$weibull=pweibull
CDFs$gamma=pgamma
CDFs$exponential=pexp
CDFs$normal=pnorm
CDFs$lognormal=plnorm
CDFs$lnorm3=plnorm3
CDFs$posnorm=pposnorm
CDFs$invgamma=pigamma
CDFs$rayleigh=prayleigh
CDFs$f=pf
CDFs$ncf=pf
CDFs$dagum=pdagum
CDFs$frechet=pfrechet
CDFs$beta=pbeta
CDFs$geometric=pgeom
CDFs$logarithmic=plog
CDFs$binom=pbinom
CDFs$poisson=ppois
CDFs$nbinom=pnbinom
CDFs$zapois=pzapois
CDFs$pospois=ppospois
CDFs$posnegbin=pposnegbin
CDFs$shiftpois=pshiftpois
CDFs$shiftnegbin=pshiftnegbin
CDFs$posgeom=pposgeom
CDFs$wrpcauchy=pwrpcauchy
CDFs$wrpnorm=pwrpnorm

#Assign appropriate CDFs for all dists
generate=vector("list",length(validDists))
names(generate)=validDists
generate$weibull=rweibull
generate$gamma=rgamma
generate$exponential=rexp
generate$normal=rnorm
generate$lognormal=rlnorm
generate$lnorm3=rlnorm3
generate$posnorm=rposnorm
generate$invgamma=rigamma
generate$rayleigh=rrayleigh
generate$f=rf
generate$ncf=rf
generate$dagum=rdagum
generate$frechet=rfrechet
generate$beta=rbeta
generate$geometric=rgeom
generate$logarithmic=rlog
generate$binom=rbinom
generate$poisson=rpois
generate$nbinom=rnbinom
generate$zapois=rzapois
generate$pospois=rpospois
generate$posnegbin=rposnegbin
generate$shiftpois=rshiftpois
generate$shiftnegbin=rshiftnegbin
generate$posgeom=rposgeom
generate$wrpcauchy=rwrpcauchy
generate$wrpnorm=rwrpnorm

#Test to see if dists are supported
for(i in 1:length(dists)){
  if(sum(is.element(validDists,dists[i]))==0)stop(paste("Dist",i,"not currently supported.  Try checking the spelling and capitalization."))
}

#Which distributions did the user request?
pick=which(is.element(validDists,dists))
ord=match(dists,validDists[pick])
pick=pick[ord]

#Extract appropriate transformations, PDFs, CDFs, etc.
transforms.used=transforms[pick]
inv.transforms.used=inv.transforms[pick]
PDFs.used=PDFs[pick]
CDFs.used=CDFs[pick]
generate.used=generate[pick]

#Assign appropriate transformations types (need to flag circular dists)
type=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2)
transforms.used$type=type[pick]
inv.transforms.used$type=type[pick]
type.used=type[pick]

#How many parameters do the distributions have? (Used to make sure user inputs the correct #)
npar=c(2,2,1,2,2,3,2,2,1,2,3,3,3,2,1,1,2,1,2,2,1,2,1,2,1,2,2)
npar.used=npar[pick]

return(list(transforms.used,inv.transforms.used,PDFs.used,CDFs.used,generate.used,type.used,npar.used))
}

