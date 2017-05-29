#' @title Correlation.Overlap
#' @description This function computes the correlation between Nodes and the Overlapping Nodes of interest.
#' @param Data data.frame containing the expression data. Genes on the Rows, Individuals on the Columns. Don't forget to give the names to the Genes and to the Individuals. Genes must have the row.names() with the Gene Name.
#' @param Overlap  A vector containg the names of the Gene Regulatory Factors (GRFs).
#' @param method  Spearman ("s", "spearman") or Pearson ("p", "pearson") correlation
#' @return List containing
#' \item{Descriptive}{the quantiles of the correlation and mean.}
#' \item{Final_Correlation}{the correlation matrix where the diagonal is 0.}
#' @rdname Correlation.Overlap
#' @export
#' @importFrom stats cor

Correlation.Overlap = function(Data, Overlap, method ){
  dfExpression = Data
  GRF = Overlap
  if (method %in% c("s", "spearman")){
    # message("\nCorrelation: Sp")


    FINAL_SP = stats::cor(t(dfExpression), method = "s",  use = "all.obs")
    diag(FINAL_SP) <- 0
    FINAL_SP[is.na(FINAL_SP)] = 0

    Final_Correlation = subset(FINAL_SP, row.names(FINAL_SP) %in% GRF)
    # Descriptive = as.data.frame((rbind( as.matrix(stats::quantile(Final_Correlation, na.rm = T)),
    #                                     as.matrix(mean(Final_Correlation, na.rm = T) ))))
    # row.names(Descriptive) = c("Spearman0",
    #                            "Spearman25",
    #                            "Spearman50",
    #                            "Spearman75",
    #                            "Spearman100",
    #                            "Spearmanmean")

  }
  if (method  %in% c("p", "pearson")){
    # message("\nCorrelation: Pearson")
    FINAL_PEARSON = stats::cor(t(dfExpression), method = "p", use = "all.obs")
    diag(FINAL_PEARSON) <- 0
    FINAL_PEARSON[is.na(FINAL_PEARSON)] = 0

    Final_Correlation = subset(FINAL_PEARSON, row.names(FINAL_PEARSON) %in% GRF)
    # Descriptive = as.data.frame((rbind( as.matrix(stats::quantile(Final_Correlation, na.rm = T)),
    #                                     as.matrix(mean(Final_Correlation, na.rm = T) ))))
    # row.names(Descriptive) = c("Pearson0",
    #                            "Pearson25",
    #                            "Pearson50",
    #                            "Pearson75",
    #                            "Pearson100",
    #                            "Pearsonmean")

  }
  output = list(Descriptive= NULL, Final_Correlation)
  return(output)
}

