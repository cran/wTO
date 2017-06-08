
#' @title wTO.aux.each
#' @description wTO.aux.each calculte the wTO for each one of the resamplings.
#' @keywords internal
#' @importFrom stats na.exclude
#' @importFrom data.table setkeyv
#' @param n Number of bootstraps / reshuffles to be run for the estimatives of the "Threshold" or "pval".
#' @param Data data.frame containing the count / expression data for the correlation.
#' @param Overlap Nodes of interested, where the Overlapping weights will be computed.
#' @param method Type of the correlation that should be used. "s" / "spearman" will compute the rank spearman correlation, "p" / "pearson" will compute the linear correlation. If no value is given, the default is to use "s".
#' @param method_resampling method of the resampling. Bootstrap or Reshuffle. Bootstrap null hypothesis is that the wTO is random, and Reshuffle tests if the wTO is equal to zero.


wTO.aux.each = function (n, Data, Overlap, method, method_resampling, lag){
  dfExpression = Data
  GRF = unique(Overlap)
  if(method_resampling == "Bootstrap"){
    real_Genes = sample(dfExpression, replace = T)
  }
  if(method_resampling == "BlockBootstrap"){
    nsampl = ifelse (ncol(dfExpression) %% lag == 0, ncol(dfExpression) %/% lag, ncol(dfExpression) %/% lag +1)
    Y = sample(1:nsampl, size = nsampl, replace =  T)
    Vect = Y*lag
    i = lag - 1
    while( i > 0){
      Vect = cbind(Vect , Y*lag - i)
      i = i - 1
    }

    SAMPLES = c(Vect)
    SAMPLES[SAMPLES > ncol(dfExpression)] <- NA
    SAMPLE = stats::na.exclude(SAMPLES)
    real_Genes =  dfExpression[,SAMPLE]
    row.names(real_Genes)=row.names(dfExpression)

  }
  else if(method_resampling == "Reshuffle" ){
    real_Genes = as.data.frame(lapply(1:ncol(dfExpression), sample_ind, dfExpression = dfExpression))
    names(real_Genes)=names(dfExpression)
    row.names(real_Genes)=row.names(dfExpression)
  }


  Saving = Correlation.Overlap(Data = real_Genes, Overlap = GRF, method = method)
  WTO_abs = wTO(A = Saving[[2]],  sign = "abs")
  WTO_sig = wTO(A = Saving[[2]],  sign = "sign")

  message(".", appendLF = F)
  Cor_star = wTO.in.line(WTO_sig)
  Cor_star_abs = wTO.in.line(WTO_abs)
  names(Cor_star) = c ("Node.1", "Node.2", "wTo_sign")
  names(Cor_star_abs) = c ("Node.1", "Node.2", "wTo_abs")
  data.table::setkeyv(Cor_star, c("Node.1", "Node.2"))
  data.table::setkeyv(Cor_star_abs, c("Node.1", "Node.2"))
  
  return(Cor_star[Cor_star_abs])
}

