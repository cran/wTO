
#' @title wTO.aux.each
#' @description wTO.aux.each calculate the wTO for each one of the resamplings.
#' @keywords internal
#' @importFrom stats na.exclude
#' @importFrom data.table setkeyv
#' @param n Number of bootstraps / reshuffles to be run for the estimates of the "Threshold" or "pval".
#' @param Data data.frame containing the count / expression data for the correlation.
#' @param Overlap Nodes of interested, where the Overlapping weights will be computed.
#' @param method Type of the correlation that should be used. "s" / "spearman" will compute the rank spearman correlation, "p" / "pearson" will compute the linear correlation. If no value is given, the default is to use "s".
#' @param method_resampling method of the resampling. Bootstrap or Reshuffle. Bootstrap null hypothesis is that the wTO is random, and Reshuffle tests if the wTO is equal to zero.


wTO.aux.each = function (n, Data, Overlap, method, method_resampling, lag, ID){
  if(method_resampling == "Bootstrap"){
    real_Genes = sample(Data, replace = T)
  }
  if(method_resampling == "BlockBootstrap"){
    if(!is.null(lag)){
      nsampl = ifelse (ncol(Data) %% lag == 0, ncol(Data) %/% lag, ncol(Data) %/% lag +1)
      Y = sample(1:nsampl, size = nsampl, replace =  T)
      Vect = Y*lag
      i = lag - 1
      while( i > 0){
        Vect = cbind(Vect , Y*lag - i)
        i = i - 1
      }
      
      SAMPLES = c(Vect)
      SAMPLES[SAMPLES > ncol(Data)] <- NA
      SAMPLE = stats::na.exclude(SAMPLES)
      real_Genes =  Data[,SAMPLE]
      row.names(real_Genes)=row.names(Data)
      
    }
    
    if(!is.null(ID)){
      ID %<>% as.factor
      bootID = sample(levels(ID), replace = TRUE)
      
      
      Data_boot = subset(Data, select = ID %in% bootID[1])
      
      for (k in 2:length(bootID)){
        real_Genes = cbind(Data_boot,
                          subset(Data, select = ID %in% bootID[k]))
      }
    }
  }
  else if(method_resampling == "Reshuffle" ){
    real_Genes = as.data.frame(lapply(1:ncol(Data), FUN = sample_ind, dfExpression = Data))
    names(real_Genes)=names(Data)
    row.names(real_Genes)=row.names(Data)
  }
  
  
  Saving = CorrelationOverlap(Data = real_Genes, Overlap = Overlap, method = method)
  WTO_abs = wTO(A = Saving,  sign = "abs")
  WTO_sig = wTO(A = Saving,  sign = "sign")
  
  # message(".", appendLF = F)
  Cor_star = wTO.in.line(WTO_sig)
  Cor_star_abs = wTO.in.line(WTO_abs)
  names(Cor_star) = c ("Node.1", "Node.2", "wTo_sign")
  names(Cor_star_abs) = c ("Node.1", "Node.2", "wTo_abs")
  data.table::setkeyv(Cor_star, c("Node.1", "Node.2"))
  data.table::setkeyv(Cor_star_abs, c("Node.1", "Node.2"))
  
  return(Cor_star[Cor_star_abs])
}

