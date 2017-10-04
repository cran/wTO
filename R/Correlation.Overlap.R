#' @title Correlation.Overlap
#' @description This function computes the correlation between Nodes and the Overlapping Nodes of interest.
#' @param Data data.frame containing the expression data. Nodes on the Rows, Individuals on the Columns. Don't forget to give the names to the Nodes and to the Individuals. Nodes must have the row.names() with the Node Name.
#' @param Overlap  A vector containg the names of the Nodes of interest.
#' @param method  Spearman ("s", "spearman") or Pearson ("p", "pearson") correlation
#' @rdname Correlation.Overlap
#' @export
#' @importFrom stats cor


Correlation.Overlap = function(Data, Overlap, method ){
  COR = stats::cor(t(Data), method = method,  use = "all.obs")
  diag(COR) <- 0
  COR[is.na(COR)] = 0
  Final_Correlation = subset(COR, row.names(COR) %in% Overlap)
  return(Final_Correlation)
}

