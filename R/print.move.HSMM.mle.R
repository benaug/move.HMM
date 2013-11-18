#'Print a move.HSMM object
#'
#'This function prints the parameter values and AICc from a move.HSMM
#'object.
#'
#'@param x A move.HSMM object containing a fitted HSMM model.
#'
#'@param digits	a non-null value for digits specifies the minimum number of significant digits to be printed in values. The default, NULL, uses getOption(digits). (For the interpretation for complex numbers see signif.) Non-integer values will be rounded down, and only values greater than or equal to 1 and no greater than 22 are accepted.
#'@param quote	logical, indicating whether or not strings (characters) should be printed with surrounding quotes.
#'@param na.print	a character string which is used to indicate NA values in printed output, or NULL (see 'Details').
#'@param print.gap	a non-negative integer ??? 1024, or NULL (meaning 1), giving the spacing between adjacent columns in printed vectors, matrices and arrays.
#'@param right	logical, indicating whether or not strings should be right aligned. The default is left alignment.
#'@param max	a non-null value for max specifies the approximate maximum number of entries to be printed. The default, NULL, uses getOption(max.print); see that help page for more details.
#'@param useSource	logical, indicating whether to use source references or copies rather than deparsing language objects. The default is to use the original source if it is available.
#'@param ...	further arguments to be passed to or from other methods. They are ignored in this function.
#'@return Parameter values and AICc from a move.HSMM object.
#'@method print move.HSMM
#'@export
print.move.HSMM=function(x, digits = NULL, quote = TRUE,
                            na.print = NULL, print.gap = NULL, right = FALSE,
                            max = NULL, useSource = TRUE, ...){
  cat("HSMM with",x$nstates,"Hidden States\n")
  cat("\nParameter estimates\n")
  print(round(x$parout,5))
  cat("\n * derived parameter \n ** 0 by def'n\n\n")
  cat("\n Approximate Stationary Distribution\n")
  print(x$delta)
  cat("\nApproximation\n\n")
  a=cbind(1:x$nstates,x$m1)
  colnames(a)=c("Aggregate","States")
  rownames(a)=rep("",x$nstates)
  print(a)
  cat("\nlogL:",x$mllk )
  cat("\n# estimated parameters:",x$npar )
  cat("\nAICc:",x$AICc,"\n" )
}