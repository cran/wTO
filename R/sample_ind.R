#' @title sample_ind
#' @description Computes the resufling of the expression values for the IDs.
#' @param x Column to be resampled
#' @param dfExpression data.frame object containing the genes expression on the rows and the individuals (Individuals in the Columns)
#' @keywords internal


sample_ind = function(x, dfExpression){
  z = base::sample(dfExpression[,x], replace = F)
  return(z)
}
