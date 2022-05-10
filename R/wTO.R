#' @title wTO
#' @description Calculates the weighted topologycal overlap (wTO)
#' between a set of Nodes and the Overlapping nodes. This function implements the method from Nowick (2009).
#' @param A Is the weighted adjency matrix (correlation matrix).
#' @param sign ("abs", "sign") if the user wants to use the absolute correlation or the signed correlation.
#' @return A matrix containing the wTO values.
#' @export
#' @references Katja Nowick, Tim Gernat, Eivind Almaas and Lisa Stubbs (2009) <doi:10.1073/pnas.0911376106>
#' @author Deisy Morselli Gysi <deisy at bioinf.uni-leipzig.de>
#' @importFrom Rfast Crossprod transpose


wTO = function(A,  sign = c("abs", "sign")){

  if(sign %in% c("abs", "absolute")){
    A = abs(A)
  }
  A_TF = as.data.frame(subset(A, select = row.names(A)))
  C = Rfast::Crossprod(A, Rfast::transpose(A))
  
  W = C + A_TF ###
  K  = matrix(NA, nrow(A_TF), ncol(A_TF))
  KI = rowSums(abs(A), na.rm = T)
  for( ii in 1: nrow(A_TF)){
    for( jj in 1: ncol(A_TF)){
      K[ii,jj] = min(KI[ii], KI[jj])
    }
  }
  WTO = round(W / (K + 1 - abs(A_TF)),3)
  return(WTO)
}
