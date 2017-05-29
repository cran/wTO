#' @title wTO
#' @description Calculates the weighted topologycal overlap (wTO)
#' between a set of Nodes and the Overlapping nodes. This function implements the method from \cite{nowick2009differences}.
#' @param A Is the weighted adjency matrix (correlation matrix).
#' @param sign ("abs", "sign") if the user wants to use the absolute correlation or the signed correlation.
#' @return A matrix containing the wTO values.
#' @export
#' @references \url{http://www.pnas.org/content/106/52/22358.full.pdf}


wTO = function(A,  sign = c("abs", "sign")){

  if(sign %in% c("abs", "absolute")){
    A = abs(A)
  }
  A_TF = as.data.frame(subset(A, select = row.names(A)))
  C = as.matrix(A) %*% t(A)

  W = C + A_TF ###
  K  = matrix(NA, nrow(A_TF), ncol(A_TF))
  KI = rowSums(abs(A), na.rm = T)
  for( ii in 1: nrow(A_TF)){
    for( jj in 1: ncol(A_TF)){
      K[ii,jj] = min(KI[ii], KI[jj])
    }
  }
  WTO = W / (K + 1 - abs(A_TF))
  return(WTO)
}
