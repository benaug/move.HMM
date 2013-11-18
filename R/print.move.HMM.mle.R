#'Print a move.HMM object
#'
#'This function prints the parameter values and AICc from a move.HMM
#'object.
#'
#'@param x A move.HMM object containing a fitted HMM model.
#'
#'@param digits	a non-null value for digits specifies the minimum number of significant digits to be printed in values. The default, NULL, uses getOption(digits). (For the interpretation for complex numbers see signif.) Non-integer values will be rounded down, and only values greater than or equal to 1 and no greater than 22 are accepted.
#'@param quote	logical, indicating whether or not strings (characters) should be printed with surrounding quotes.
#'@param na.print	a character string which is used to indicate NA values in printed output, or NULL (see 'Details').
#'@param print.gap	a non-negative integer ??? 1024, or NULL (meaning 1), giving the spacing between adjacent columns in printed vectors, matrices and arrays.
#'@param right	logical, indicating whether or not strings should be right aligned. The default is left alignment.
#'@param max	a non-null value for max specifies the approximate maximum number of entries to be printed. The default, NULL, uses getOption(max.print); see that help page for more details.
#'@param useSource	logical, indicating whether to use source references or copies rather than deparsing language objects. The default is to use the original source if it is available.
#'@param ...	further arguments to be passed to or from other methods. They are ignored in this function.
#'@return Parameter values and AICc from a move.HMM object.
#'@method print move.HMM
#'@export
print.move.HMM=function(x, digits = NULL, quote = TRUE,
                            na.print = NULL, print.gap = NULL, right = FALSE,
                            max = NULL, useSource = TRUE, ...){
  cat("HMM with",x$nstates,"Hidden States\n")
  cat("\nParameter estimates (* are derived from other parameters)\n")
  print(round(x$parout,5))
  cat("\nlogL:",x$mllk )
  cat("\n# estimated parameters:",x$npar )
  cat("\nAICc:",x$AICc ,"\n")
}