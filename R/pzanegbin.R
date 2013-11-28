#'Zero-altered Negative Binomial CDF
#'
#'This function evaluates zero-altered negative binomial cdf.  This function is modified
#'from the pzanegbin VGAM package so that munb is not a parameter input option. 
#'
#'@param q a vector of values where the cdf is to be evaluated
#'@param size a value for the ordinary Negative Binomial size parameter
#'@param prob  a value for the ordinary Negative Binomial prob parameter
#'@param pobs0 Probability of zero.
#'@return A vector of zero-altered negative binomial pdf values
#'@export
pzanegbin=function (q, size, prob = NULL, pobs0 = 0) {
  LLL <- max(length(q), length(pobs0), length(prob), length(size))
  if (length(q) != LLL) 
    q <- rep(q, len = LLL)
  if (length(pobs0) != LLL) 
    pobs0 <- rep(pobs0, len = LLL)
  if (length(prob) != LLL) 
    prob <- rep(prob, len = LLL)
  if (length(size) != LLL) 
    size <- rep(size, len = LLL)
  ans <- rep(0, len = LLL)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1)) 
    stop("argument 'pobs0' must be in [0,1]")
  qindex <- (q > 0)
  ans[qindex] <- pobs0[qindex] + (1 - pobs0[qindex]) * pposnegbin(q[qindex], 
                                                                  size = size[qindex], prob = prob[qindex])
  ans[q < 0] <- 0
  ans[q == 0] <- pobs0[q == 0]
  ans
}
