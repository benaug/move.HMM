#'Zero-altered Negative Binomial PDF
#'
#'This function evaluates zero-altered negative binomial pdf.  This function is modified
#'from the dzanegbin VGAM package so that munb is not a parameter input option. 
#'
#'@param x a vector of values where the pdf is to be evaluated
#'@param size a value for the ordinary Negative Binomial size parameter
#'@param prob  a value for the ordinary Negative Binomial prob parameter
#'@param pobs0 Probability of zero.
#'@return A vector of zero-altered negative binomial pdf values
#'@export
dzanegbin=function (x, size, prob = NULL, pobs0 = 0, log = FALSE){
  if (!is.logical(log.arg <- log) || length(log) != 1) 
    stop("bad input for argument 'log'")
  rm(log)
  LLL <- max(length(x), length(pobs0), length(prob), length(size))
  if (length(x) != LLL) 
    x <- rep(x, len = LLL)
  if (length(pobs0) != LLL) 
    pobs0 <- rep(pobs0, len = LLL)
  if (length(prob) != LLL) 
    prob <- rep(prob, len = LLL)
  if (length(size) != LLL) 
    size <- rep(size, len = LLL)
  ans <- rep(0, len = LLL)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1)) 
    stop("argument 'pobs0' must be in [0,1]")
  if (!is.Numeric(prob, positive = TRUE)) 
    stop("argument 'prob' must be in (0,Inf)")
  if (!is.Numeric(size, positive = TRUE)) 
    stop("argument 'size' must be in (0,Inf)")
  index0 <- x == 0
  if (log.arg) {
    ans[index0] <- log(pobs0[index0])
    ans[!index0] <- log1p(-pobs0[!index0]) + dposnegbin(x[!index0], 
                                                        prob = prob[!index0], size = size[!index0], log = TRUE)
  }
  else {
    ans[index0] <- pobs0[index0]
    ans[!index0] <- (1 - pobs0[!index0]) * dposnegbin(x[!index0], 
                                                      prob = prob[!index0], size = size[!index0])
  }
  ans
}