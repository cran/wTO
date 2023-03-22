#' @title CorrelationOverlap
#' @description This function computes the correlation between Nodes and the Overlapping Nodes of interest.
#' @param Data data.frame containing the expression data. Nodes on the Rows, Individuals on the Columns. Don't forget to give the names to the Nodes and to the Individuals. Nodes must have the row.names() with the Node Name.
#' @param Overlap  A vector containg the names of the Nodes of interest.
#' @param method  Spearman ("s", "spearman") or Pearson ("p", "pearson") correlation
#' @rdname CorrelationOverlap
#' @author Deisy Morselli Gysi <deisy at bioinf.uni-leipzig.de>

#' @export
#' @importFrom stats cor
#' @importFrom HiClimR fastCor
#' @importFrom Rfast transpose


CorrelationOverlap = function(Data, Overlap, method ){
  
  Overlap = as.character(Overlap)
  if(method %in% c("pearson", "p")){
    COR = suppressMessages(suppressWarnings(HiClimR::fastCor(t((as.matrix(Data))), 
                                                             nSplit = 10, 
                                                             upperTri = FALSE,  verbose = F)
                                            ))
  } else{
    COR = suppressMessages(suppressWarnings(stats::cor(t(Data), 
                                                       method = method,  
                                                       use = "pairwise.complete.obs")))  
  }
  
  diag(COR) <- 0
  COR[is.na(COR)] = 0
  # Final_Correlation = subset(COR, row.names(COR) %in% Overlap)
  # 
  Final_Correlation = COR[row.names(COR) %in% Overlap, colnames(COR) %in% Overlap]
  return(Final_Correlation)
}

