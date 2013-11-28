
##' @S3method coef move.HSMM
##' @S3method coef move.HMM
coef.move.HMM <- coef.move.HSMM <-
    function(object, fmt=c("list","vector","matrix"), ...) {
        fmt <- match.arg(fmt)
        switch(fmt,list=object[["params"]],
               vector=object[["params"]][,1],
               matrix=object[["params"]])
    }

simulateMove <- function(object, nsim=1, seed=NULL,
                         nsteps=nrow(object[["obs"]]),
                         delta=object[["delta"]],
                         newparams=NULL,
                         ...)  ## for generic compatibility
{
        if (!is.null(seed)) set.seed(seed)
        if (nsim>1) return(replicate(nsim,simulate(object,delta=delta,
                                                   newparams=newparams),
                                     simplify=FALSE))
        ## TO DO: check that format of 'newparams' is correct ...
        if (!is.null(newparams))
            object$params <- newparams
        mfun <- switch(class(object),
                       move.HMM=move.HMM.simulate,
                       move.HSMM=move.HSMM.simulate)
        res <- with(object,mfun(dists,params,n=nsteps,
                                nstates,delta))
        ## now convert the results to a data frame
        with(res,setNames(data.frame(obs,states),
                      c("steps","turns","state")))
    }

##' @S3method simulate move.HSMM
##' @S3method simulate move.HMM

simulate.move.HSMM <- simulate.move.HMM <- simulateMove

