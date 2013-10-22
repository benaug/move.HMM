moveHMM=function(x,...)UseMethod("moveHMM")

moveHMM.default <- function(obs,dists,params,nstates,stepm,iterlim=200,turn=NULL,hessian=FALSE){
  #Get appropriate linearizing transformations and PDFs 
  out=get.Transforms.PDFs(dists,nstates,turn)
  transforms=out[[1]]
  inv.transforms=out[[2]]
  PDFs=out[[3]]
  skeleton=params
  parvect <- move.HMM.pn2pw(transforms,params,nstates)
  
  
  mod <- nlm(move.HMM.mllk,parvect,obs,print.level=2,stepmax=stepm,PDFs=PDFs,skeleton=skeleton,inv.transforms=inv.transforms,nstates=nstates,iterlim=iterlim,hessian=F)
  mllk <- mod$minimum
  pn <- move.HMM.pw2pn(inv.transforms,mod$estimate,skeleton,nstates)
  params=pn$params
  k=length(parvect)
  AIC=2*k+2*mllk
  AICc=AIC+(2*k*(k+1))/(nrow(obs)-k-1)
  list(params=params,Gamma=pn$Gamma,delta=pn$delta,mllk=mllk,AICc=AICc)
}

print.moveHMM=function(x,...){
  print(params)
  print(AICc)
}

iterlim=200
out=moveHMM(obs,dists,params,nstates,stepm,iterlim,turn)


files=c("move.HMM.pw2pn.R","move.HMM.pn2pw.R","move.HMM.mllk.R","move.HMM.mle.R","move.HMM.lalpha.R","move.HMM.psres.R","pwrpcauchy.R","pwrpnorm.R","dwrpnorm.R","get.Transforms.PDFs.R","psresid.R","move.HMM.lalphabeta.R","move.HMM.state_probs.R")
  package.skeleton(name="moveHMM",code_files=files)