#' @title wTO.Complete
#' @author Deisy Morselli Gysi <deisy at bioinf.uni-leipzig.de>

#' @param k Number of threads to be used for computing the weight Topological Overlap. Default is set to 1.
#' @param n Number of resamplings, used to compute the empirical distribuitions of the links. Default is set to 100.
#' @param Data data.frame containing the count / expression data for the correlation.
#' @param Overlap Set of nodes of interest, where the Overlapping weights will be computed.
#' @param method Type of the correlation that should be used. "s" / "spearman" will compute the rank spearman correlation, "p" / "pearson" will compute the linear correlation. If no value is given, the default is to use "p".
#' @param method_resampling method of the resampling. Bootstrap, BlockBootstrap or Reshuffle. Bootstrap null hypothesis is that the wTO is random, and Reshuffle tests if the wTO is equal to zero.
#' @param pvalmethod method to compute the multiple test correction for the pvalue. for more information check the function \code{\link[stats]{p.adjust}}.
#' @param savecor T/F if need to save the correlation.
#' @param expected.diff Difference expected between the real wTO and resampled wTO By default, it is set to 0.2.
#' @param lag time dependency, lag, if you are using the BlockedBootstrap.
#' @param ID ID of the samples for the blocked bootstrap (for repeated measures).

#' @param normalize T/F Should the data be normalized?
#' @param plot T/F Should the diagnosis plot be plotted?
#'
#' @description Compute the wTO and also the bootstraps. Proposed at: 	arXiv:1711.04702
#' @return a list with results.
#' \itemize{
#' \item wTO is a data.frame containig the Nodes, the wTO computed using the signed correlations, the pvalue and the adj.pvalue.
#' \item abs.wTO is a data.frame containig the Nodes, the wTO computed using the absolute correlations, the pvalue and the adj.pvalue.
#' \item Correlation is a data.frame containing the correlation between all the nodes.
#' \item Empirical.Quantile quantile values for the empirical distribution.
#' \item Quantile quantile values for the sample distribution.
#' }
#' @importFrom  parallel makeCluster clusterExport clusterApplyLB  stopCluster
#' @importFrom  data.table rbindlist dcast
#' @importFrom  som normalize
#' @importFrom  stats cor p.adjust reshape
#' @importFrom  graphics plot axis par abline legend



#'
#'
#' @examples
#' \dontrun{
#' # Using spearman rank correlation and bonferroni correction for the pvalues.
#' wTO.Complete( k =8, n = 1000, Data = Microarray_Expression1,
#'  Overlap = ExampleGRF$x, method = "s", pvalmethod = "bonferroni")
#'  # Changing the resampling method to Reshuffle.
#' wTO.Complete( k =1, n = 1000, Data = Microarray_Expression1,
#' Overlap = ExampleGRF$x, method_resampling = "Reshuffle")
#'  # Changing the resampling method to BlockBootstrap, with a lag of 2.
#'  row.names(metagenomics_abundance) = metagenomics_abundance$OTU
#'  metagenomics_abundance = metagenomics_abundance[,-1]
#' wTO.Complete( k =1, n = 1000, Data = metagenomics_abundance, method = "s",
#' Overlap = row.names(metagenomics_abundance), method_resampling = "BlockBootstrap", lag = 2)
#' wTO.Complete( k =2, n = 1000, Data = Microarray_Expression1, method = "s",
#' Overlap = ExampleGRF$x, method_resampling = "BlockBootstrap", ID = rep(1:9,each = 2))
#' X = wTO.Complete( k =1, n = 1000, Data = Microarray_Expression1,
#' Overlap = ExampleGRF$x, method = "p", plot = FALSE)
#'  }

#' @export




wTO.Complete = function(k = 1 ,n = 100, Data , Overlap = row.names(Data),
                        method = "p", method_resampling = "Bootstrap",
                        pvalmethod = "BH", savecor = F,
                        expected.diff = 0.20, lag = NULL, ID = NULL, 
                        normalize = F, plot = T){
  N = k
  Overlap = unique(as.character(Overlap))
  `%ni%` <- Negate(`%in%`)
  ##### Messages
  if(is.numeric(k) == F){
    stop("k must be numeric.")
  }
  if(k <= 0){
    stop("k must be greater than 0.")
  }
  if(is.numeric(n) == F){
    stop("n must be numeric.")
  }
  if(n <= 0){
    stop("n must be greater than 0.")
  }
  if(is.data.frame(Data) == F){
    stop("Data must be a data.frame.")
  }
  
  if(method %ni% c("s", "spearman", "p", "pearson")){
    stop('Method must be: "s", "spearman", "p" or "pearson".')
  }
  
  if(method_resampling %ni% c("Bootstrap", "Reshuffle", "BlockBootstrap")){
    stop('Method must be: "Bootstrap", "BlockBootstrap" or "Reshuffle".')
  }
  if(method_resampling %in% "BlockBootstrap"){
    if (is.null(lag)&is.null(ID)){
      stop('If you want to use the "BlockBootstrap" please give a lag or the indivuals ID.')
    }
    if(!is.null(lag)&!is.null(ID)){
      stop('If you want to use the "BlockBootstrap" please give a lag OR the indivuals ID.')
    }
  }
  if(pvalmethod %ni% c ('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none')){
    stop("pvalmethod must be:  'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr' or 'none'")
  }
  
  if(normalize %ni% c (T, F)){
    stop("normalize must be:  TRUE or FALSE.")
  }
  
  if(normalize == T){
    Data.n = as.data.frame(som::normalize(Data))
    row.names(Data.n)= row.names(Data)
    Data = Data.n
  }
  
  DIM_Overlap = nrow(subset(Data, row.names(Data) %in% Overlap))
  if(DIM_Overlap == 0){
    stop('There is no overlapping nodes. Please check your input "Overlap"')
  }
  if(!is.null(DIM_Overlap)){
    message(paste('There are',DIM_Overlap, "overlapping nodes,",dim(Data)[1],
                  "total nodes and" , dim(Data)[2],"individuals." ))
  }
  
  
  message("This function might take a long time to run. Don't turn off the computer.")
  PAR = par()
   
  ## For the original data
  # real_Genes = Data
  Saving = CorrelationOverlap(Data = Data, Overlap = Overlap, method = method)
  WTO_abs = wTO(A = Saving,  sign = "abs")
  WTO_sign = wTO(A = Saving,  sign = "sign")
  Cor_real = wTO.in.line(WTO_sign)
  Cor_real_abs = wTO.in.line(WTO_abs)
  names(Cor_real) = names(Cor_real_abs) <-c("Node.1", "Node.2", "wTO_0")
  idcol = c("Node.1", "Node.2")
  rm("WTO_abs")
  rm("WTO_sign")
  data.table::setkeyv(Cor_real, c("Node.1", "Node.2"))
  data.table::setkeyv(Cor_real_abs, c("Node.1", "Node.2"))
  
  Orig = cbind(Rep = 0, Cor_real[Cor_real_abs])
  names(Orig)= c("Rep","Node.1", "Node.2", "wTO_sign", "wTO_abs")
  reps_rest = n
  ### If only one node
  if ( k == 1){
    a = 0
    while ( reps_rest  > 0){
      # message(a)
      K = 1:min(N, reps_rest)
      
      # K = 1:n
      OUTPUT = lapply(K, wTO.aux.each, Data= Data,
                      Overlap = Overlap, method = method, ID, lag = lag, method_resampling= method_resampling)
      
      ALL  =  data.table::rbindlist(OUTPUT, idcol = idcol)
      names(ALL) = names(Orig) =  c("Rep", "Node.1", "Node.2", "wTO_sign" ,"wTO_abs")
      
      ALL_DT_sig = data.table::dcast(ALL, Node.1 + Node.2  ~ Rep, value.var = "wTO_sign")
      ALL_DT_abs = data.table::dcast(ALL, Node.1 + Node.2  ~ Rep, value.var = "wTO_abs")
      
      if ( a == 0){
        
        Ps1 = rowSums(ALL_DT_sig[,-c(1:2)] < Orig$wTO_sign - expected.diff)
        Ps2 = rowSums(ALL_DT_sig[,-c(1:2)] > Orig$wTO_sign + expected.diff)
        Ps = Ps1 + Ps2
        
        Pa1 = rowSums(ALL_DT_abs[,-c(1:2)] < Orig$wTO_abs - expected.diff)
        Pa2 = rowSums(ALL_DT_abs[,-c(1:2)] > Orig$wTO_abs + expected.diff)
        Pa = Pa1 + Pa2
        
        TAB_SIGN = as.data.frame(table(unlist(round(ALL_DT_sig[,-c(1:2)], 2))))
        TAB_ABS = as.data.frame(table(unlist(round(ALL_DT_abs[,-c(1:2)], 2))))
      }
      if ( a > 0){
        Ps1 = rowSums(ALL_DT_sig[,-c(1:2)] < Orig$wTO_sign - expected.diff)
        Ps2 = rowSums(ALL_DT_sig[,-c(1:2)] > Orig$wTO_sign + expected.diff)
        Ps = Ps + Ps1 + Ps2
        
        Pa1 = rowSums(ALL_DT_abs[,-c(1:2)] < Orig$wTO_abs - expected.diff)
        Pa2 = rowSums(ALL_DT_abs[,-c(1:2)] > Orig$wTO_abs + expected.diff)
        Pa = Pa + Pa1 + Pa2
        
        TAB_SIGN_aux = as.data.frame(table(unlist(round(ALL_DT_sig[,-c(1:2)], 2))))
        TAB_ABS_aux = as.data.frame(table(unlist(round(ALL_DT_abs[,-c(1:2)], 2))))
        
        TAB_SIGN = plyr::join(TAB_SIGN, TAB_SIGN_aux, by = "Var1")
        TAB_SIGN = data.frame(Var1 = TAB_SIGN$Var1,
                              Sum = rowSums(TAB_SIGN[,-1]))
        TAB_ABS = plyr::join(TAB_ABS, TAB_ABS_aux, by = "Var1")
        TAB_ABS = data.frame(Var1 = TAB_ABS$Var1,
                             Sum = rowSums(TAB_ABS[,-1]))
      }
      rm("ALL_DT_sig", "ALL_DT_abs", "ALL", "OUTPUT")
      
      reps_rest = (reps_rest - N)
      a = a +1
      
    }
  }
  else if ( k > 1){
    
    WTO = new.env()
    assign("Data", Data, envir = WTO)
    assign("Overlap", Overlap, envir = WTO)
    assign("method", method, envir = WTO)
    assign("CorrelationOverlap", CorrelationOverlap, envir = WTO)
    assign("wTO", wTO, envir = WTO)
    assign("wTO.in.line", wTO, envir = WTO)
    assign("wTO.aux.each", wTO.aux.each, envir = WTO)
    assign("method_resampling", method_resampling, envir = WTO)
    assign("sample_ind", sample_ind, envir = WTO)
    assign("lag", lag, envir = WTO)
    assign("ID", ID, envir = WTO)
    cl = parallel::makeCluster(k)
    parallel::clusterExport(cl, "Data", envir = WTO)
    parallel::clusterExport(cl, "wTO.in.line", envir = WTO)
    parallel::clusterExport(cl, "lag", envir = WTO)
    parallel::clusterExport(cl, "Overlap", envir = WTO)
    parallel::clusterExport(cl, "method", envir = WTO)
    parallel::clusterExport(cl, "CorrelationOverlap", envir = WTO )
    parallel::clusterExport(cl, "wTO", envir = WTO)
    parallel::clusterExport(cl, 'wTO.aux.each', envir = WTO)
    parallel::clusterExport(cl, 'method_resampling', envir = WTO)
    parallel::clusterExport(cl, 'sample_ind', envir = WTO)
    # message("cluster")
    # K = 1:n
    
    a = 0
    while ( reps_rest  > 0){
      # message(a)
      K = 1:min(N, reps_rest)
      
      OUTPUT = parallel::clusterApply(cl, K, wTO.aux.each , Data= Data,
                                      Overlap = Overlap, ID, lag = lag, method = method, method_resampling= method_resampling)
      ALL  =  data.table::rbindlist(OUTPUT, idcol = idcol)
      names(ALL) = names(Orig) =  c("Rep", "Node.1", "Node.2", "wTO_sign" ,"wTO_abs")
      
      ALL_DT_sig = data.table::dcast(ALL, Node.1 + Node.2  ~ Rep, value.var = "wTO_sign")
      ALL_DT_abs = data.table::dcast(ALL, Node.1 + Node.2  ~ Rep, value.var = "wTO_abs")
      
      
      if ( a == 0){
        
        Ps = rowSums(ALL_DT_sig[,-c(1:2)] < Orig$wTO_sign - expected.diff) +
          rowSums(ALL_DT_sig[,-c(1:2)] > Orig$wTO_sign + expected.diff)
        
        Pa = rowSums(ALL_DT_abs[,-c(1:2)] < Orig$wTO_abs - expected.diff) +
          rowSums(ALL_DT_abs[,-c(1:2)] > Orig$wTO_abs + expected.diff)
        
        
        TAB_SIGN = as.data.frame(table(unlist(round(ALL_DT_sig[,-c(1:2)], 2))))
        TAB_ABS = as.data.frame(table(unlist(round(ALL_DT_abs[,-c(1:2)], 2))))
      }
      if ( a > 0){
        Ps = Ps + rowSums(ALL_DT_sig[,-c(1:2)] < Orig$wTO_sign - expected.diff) +
          rowSums(ALL_DT_sig[,-c(1:2)] > Orig$wTO_sign + expected.diff)
        Pa = Pa + rowSums(ALL_DT_abs[,-c(1:2)] < Orig$wTO_abs - expected.diff) +
          rowSums(ALL_DT_abs[,-c(1:2)] > Orig$wTO_abs + expected.diff)
        
        # message(Pa)
        # message(Ps)
        
        TAB_SIGN_aux = as.data.frame(table(unlist(round(ALL_DT_sig[,-c(1:2)], 2))))
        TAB_ABS_aux = as.data.frame(table(unlist(round(ALL_DT_abs[,-c(1:2)], 2))))
        
        TAB_SIGN = plyr::join(TAB_SIGN, TAB_SIGN_aux, by = "Var1")
        TAB_SIGN = data.frame(Var1 = TAB_SIGN$Var1,
                              Sum = rowSums(TAB_SIGN[,-1]))
        TAB_ABS = plyr::join(TAB_ABS, TAB_ABS_aux, by = "Var1")
        TAB_ABS = data.frame(Var1 = TAB_ABS$Var1,
                             Sum = rowSums(TAB_ABS[,-1]))
      }
      rm("ALL_DT_sig", "ALL_DT_abs", "ALL", "OUTPUT")
      
      reps_rest = (reps_rest - N)
      a = a +1
    }
    
    parallel::stopCluster(cl)
  }
  
  message("Simulations are done.")
  
  message("Computing p-values")
  Orig$pval_sig = Ps / n
  Orig$pval_abs = Pa / n
  
  
  if(method_resampling == "Reshuffle"){
    Orig$pval_sig = 1- Orig$pval_sig
    Orig$pval_abs = 1- Orig$pval_abs
  }
  
  Orig$Padj_sig = (stats::p.adjust(Orig$pval_sig, method = pvalmethod))
  Orig$Padj_abs = (stats::p.adjust(Orig$pval_abs, method = pvalmethod))
  
  
  ## Running the correlation
  if( savecor == T){
    Total_Correlation = as.data.frame(stats::cor(t(Data), method = method))
    Total_Correlation = wTO.in.line(Total_Correlation)
    names(Total_Correlation) = c("Node.1", "Node.2", "Cor")
  }
  if( savecor == F){
    Total_Correlation = NULL
  }
  
  TAB_SIGN_aux = as.data.frame(table(round(Orig$wTO_sign,2)))
  TAB_ABS_aux = as.data.frame(table(round(Orig$wTO_abs,2)))
  TAB_SIGN = plyr::join(TAB_SIGN, TAB_SIGN_aux, by = "Var1")
  TAB_ABS = plyr::join(TAB_ABS, TAB_ABS_aux, by = "Var1")
  
  message("Computing cutoffs")
  if(plot == TRUE){
    graphics::par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(3,2))
    
  }
  
  Cutoffs = Cut.off(TAB_SIGN, "wTO - Resampling", plot = plot)
  Cutoffs_abs = Cut.off(TAB_ABS, "|wTO| - Resampling",  plot = plot)
  
  
  Orig = Orig[, -"Rep"]
  
  
  Orig$wTO_abs = as.numeric(Orig$wTO_abs)
  Orig$wTO_sign = as.numeric(Orig$wTO_sign)
  Orig$pval_abs = as.numeric(Orig$pval_abs)
  Orig$pval_sig = as.numeric(Orig$pval_sig)
  Orig$Padj_abs = as.numeric(Orig$Padj_abs)
  Orig$Padj_sig = as.numeric(Orig$Padj_sig)
  
  Quantiles = rbind(
    Cutoffs$Empirical.Quantile,
    Cutoffs$Quantile ,
    Cutoffs_abs$Empirical.Quantile,
    Cutoffs_abs$Quantile)
  row.names(Quantiles) = c(  'Empirical.Quantile',
                             'Quantile',
                             'Empirical.Quantile.abs',
                             'Quantile.abs')
  
  tQ = as.data.frame(t(Quantiles))
  output = list(wTO = Orig,
                Correlation = Total_Correlation,
                Quantiles = Quantiles
  )
  col = ifelse(Orig$pval_sig < 0.05 & Orig$pval_abs < 0.05, "red",
               ifelse(Orig$pval_sig < 0.05, "orange",
                      ifelse (Orig$pval_abs < 0.05, "yellow", "black")))
  if(plot == T){
    # par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(3,1))
    graphics::plot(Orig$wTO_sign, Orig$wTO_abs, axes = F,
                   xlab = "|wTO|", ylab = "wTO",
                   main = "|wTO| vs wTO", pch = ".", xlim = c(-1,1), ylim = c(0,1),
                   col.main = "steelblue2", col.lab = "steelblue2", col = col)
    graphics::axis(1, las = 1, cex.axis = 0.6, col = "steelblue",
                   
                   col.ticks = "steelblue3", col.axis = "steelblue")
    graphics::axis(2, las = 1, cex.axis = 0.6, col = "steelblue",col.ticks = "steelblue3", col.axis = "steelblue")
    
    graphics::legend(c(0.9,0), c ("p-value < 0.05", 'wTO sign & |wTO|',
                                  'wTO sign','|wTO|'),
                     col = c("transparent","red", "orange", "yellow"), pch = 16, bty = "n",
                     inset=c(-0.8,0), cex = 0.5 )
    
    
    graphics::par(xpd=FALSE)
    graphics::abline( h = 0,  lty = 2, col = "gray50")
    graphics::abline(v = 0,  lty = 2, col = "gray50")
    
    
    graphics::plot(Orig$wTO_sign, Orig$pval_sig, axes = F,
                   xlab = "wTO", ylab = "p-value", ylim = c(0,1), xlim = c(-1,1), col.main = "steelblue2", col.lab = "steelblue2",
                   main = "wTO vs p-value",
                   pch = 16)
    graphics::axis(1, las = 1, cex.axis = 0.6, col = "steelblue",
                   
                   col.ticks = "steelblue3", col.axis = "steelblue")
    graphics::axis(2, las = 1, cex.axis = 0.6, col = "steelblue",col.ticks = "steelblue3", col.axis = "steelblue")
    
    graphics::par(xpd=FALSE)
    graphics::abline( v = tQ$Empirical.Quantile,  lty = 2, col = c("red", "orange", "yellow", "yellow", "orange", "red"))
    graphics::par(xpd=T)
    graphics::legend(c(0.9,0), c ("Empirical Quantiles", '0.1%','2.5%','10%','90%','97.5%','99.9%'),
                     col = c("white", "red", "orange", "yellow", "yellow", "orange", "red"), lwd = 2, bty = "n",
                     inset=c(-0.8,0), cex = 0.5 )
    
    graphics::par(xpd=FALSE)
    graphics::plot(Orig$wTO_abs, Orig$pval_abs, axes = F,
                   xlab = "|wTO|", ylab = "p-value", ylim = c(0,1), xlim = c(0,1),
                   main = "|wTO| vs p-value",
                   pch = 16, col.main = "steelblue2", col.lab = "steelblue2")
    graphics::axis(1, las = 1, cex.axis = 0.6, col = "steelblue",
                   
                   col.ticks = "steelblue3", col.axis = "steelblue")
    graphics::axis(2, las = 1, cex.axis = 0.6, col = "steelblue",col.ticks = "steelblue3", col.axis = "steelblue")
    
    
    graphics::abline( v = tQ$Empirical.Quantile.abs,  lty = 2, col = c("red", "orange", "yellow", "yellow", "orange", "red"))
    graphics::par(xpd=T)
    graphics::legend(c(0.9,0), c ("Empirical Quantiles", '0.1%','2.5%','10%','90%','97.5%','99.9%'),
                     col = c("white", "red", "orange", "yellow", "yellow", "orange", "red"), lwd = 2, bty = "n",
                     inset=c(-0.8,0), cex = 0.5 )
  }
  class(output)<- append('wTO', class(output))
  message("Done!")
  return(output)
}