#' @title wTO.Complete

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
#' @param normalize T/F Should the data be normalized?
#'
#'
#' @description Compute the wTO and also the bootstraps.
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
#' wTO.Complete( k =8, n = 1000, Data = ExampledfExpression,
#'  Overlap = ExampleGRF$x, method = "s", pvalmethod = "bonferroni")
#'  # Changing the resampling method to Reshuffle.
#' wTO.Complete( k =1, n = 20, Data = ExampledfExpression,
#' Overlap = ExampleGRF$x, method_resampling = "Reshuffle")
#'  # Changing the resampling method to BlockBootstrap, with a lag of 2.
#' wTO.Complete( k =1, n = 20, Data = ExampledfExpression, method = "s",
#' Overlap = ExampleGRF$x, method_resampling = "BlockBootstrap", lag = 2)
#'  }
#' wTO.Complete( k =2, n = 8, Data = ExampledfExpression,
#' Overlap = ExampleGRF$x, method = "p")
#' @export




wTO.Complete = function(k = 1 ,n = 100, Data , Overlap ,
                        method = "p", method_resampling = "Bootstrap",
                        pvalmethod = "BH", savecor = F,
                        expected.diff = 0.20, lag = NULL,
                        normalize = F){
  dfExpression = Data
  rm(Data)
  GRF = as.character(Overlap)
  rm(Overlap)
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
  if(is.data.frame(dfExpression) == F){
    stop("dfExpression must be a data.frame.")
  }
  
  if(method %ni% c("s", "spearman", "p", "pearson")){
    stop('Method must be: "s", "spearman", "p" or "pearson".')
  }
  
  if(method_resampling %ni% c("Bootstrap", "Reshuffle", "BlockBootstrap")){
    stop('Method must be: "Bootstrap", "BlockBootstrap" or "Reshuffle".')
  }
  if(method_resampling %in% "BlockBootstrap"){
    if (is.null(lag)){
      stop('If you want to use the "BlockBootstrap" please give a lag.')
    }
    
  }
  if(pvalmethod %ni% c ('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none')){
    stop("pvalmethod must be:  'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr' or 'none'")
  }
  
  if(normalize %ni% c (T, F)){
    stop("normalize must be:  TRUE or FALSE.")
  }
  
  if(normalize == T){
    dfExpression.n = as.data.frame(som::normalize(dfExpression))
    row.names(dfExpression.n)= row.names(dfExpression)
    dfExpression = dfExpression.n
  }
  
  DIM_GRF = nrow(subset(dfExpression, row.names(dfExpression) %in% GRF))
  if(DIM_GRF == 0){
    stop('There is no overlapping nodes. Please check your input "Overlap"')
  }
  if(!is.null(DIM_GRF)){
    message(paste('There are',DIM_GRF, "overlapping nodes,",dim(dfExpression)[1],
                  "total nodes and" , dim(dfExpression)[2],"individuals." ))
  }
  
  message("This function might take a long time to run. Don't turn off the computer.")
  PAR = par()
  graphics::par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(3,2))
  
  ## For the original data
  # real_Genes = dfExpression
  Saving = Correlation.Overlap(Data = dfExpression, Overlap = GRF, method = method)
  WTO_abs = wTO(A = Saving[[2]],  sign = "abs")
  WTO_sig = wTO(A = Saving[[2]],  sign = "sign")
  Cor_real = wTO.in.line(WTO_sig)
  Cor_real_abs = wTO.in.line(WTO_abs)
  names(Cor_real) = names(Cor_real_abs) <-c("Node.1", "Node.2", "wTO_0")
  
  ### If only one node
  if ( k == 1){
    K = 1:n
    OUTPUT = lapply(K, wTO.aux.each, Data= dfExpression,
                    Overlap = GRF, method = method, lag = lag, method_resampling= method_resampling)
  }
  else if ( k > 1){
    WTO = new.env()
    assign("dfExpression", dfExpression, envir = WTO)
    assign("GRF", GRF, envir = WTO)
    assign("method", method, envir = WTO)
    assign("Correlation.Overlap", Correlation.Overlap, envir = WTO)
    assign("wTO", wTO, envir = WTO)
    assign("wTO.in.line", wTO, envir = WTO)
    assign("wTO.aux.each", wTO.aux.each, envir = WTO)
    assign("method_resampling", method_resampling, envir = WTO)
    assign("sample_ind", sample_ind, envir = WTO)
    assign("lag", lag, envir = WTO)
    cl = parallel::makeCluster(k)
    parallel::clusterExport(cl, "dfExpression", envir = WTO)
    parallel::clusterExport(cl, "wTO.in.line", envir = WTO)
    parallel::clusterExport(cl, "lag", envir = WTO)
    parallel::clusterExport(cl, "GRF", envir = WTO)
    parallel::clusterExport(cl, "method", envir = WTO)
    parallel::clusterExport(cl, "Correlation.Overlap", envir = WTO )
    parallel::clusterExport(cl, "wTO", envir = WTO)
    parallel::clusterExport(cl, 'wTO.aux.each', envir = WTO)
    parallel::clusterExport(cl, 'method_resampling', envir = WTO)
    parallel::clusterExport(cl, 'sample_ind', envir = WTO)
    # message("cluster")
    K = 1:n
    OUTPUT = parallel::clusterApplyLB(cl, K, wTO.aux.each , Data= dfExpression,
                                      Overlap = GRF, lag = lag, method = method, method_resampling= method_resampling)
    parallel::stopCluster(cl)
  }
  idcol = c("Node.1", "Node.2")
  message("Simulations are done.")
  data.table::setkeyv(Cor_real, c("Node.1", "Node.2"))
  data.table::setkeyv(Cor_real_abs, c("Node.1", "Node.2"))

  Orig = cbind(Rep = 0, Cor_real[Cor_real_abs])
  names(Orig)= c("Rep","Node.1", "Node.2", "wTO_sign", "wTO_abs")
  ALL  =  data.table::rbindlist(OUTPUT, idcol = idcol)
  
  rm(OUTPUT)
  names(ALL) = names(Orig) =  c("Rep", "Node.1", "Node.2", "wTO_sign" ,"wTO_abs")
  ALL = rbind(Orig, ALL)
  ALL_DT_sig = data.table::dcast(ALL, Node.1 + Node.2  ~ Rep, value.var = "wTO_sign")
  ALL_DT_abs = data.table::dcast(ALL, Node.1 + Node.2  ~ Rep, value.var = "wTO_abs")
  
  real_sig= as.numeric(as.matrix(ALL_DT_sig[,3]))
  real_abs= as.numeric(as.matrix(ALL_DT_abs[,3]))
  
  Z_sig = matrix(as.numeric(as.matrix(ALL_DT_sig[, - c(1:3)])), ncol = n)
  Z_abs = matrix(as.numeric(as.matrix(ALL_DT_abs[, - c(1:3)])), ncol = n)
 
  
  message("Computing p-values")
  P1 = rowSums(Z_sig < real_sig - expected.diff)
  P2 = rowSums(Z_sig > real_sig + expected.diff)
  Orig$pval_sign = (P1 + P2) / n
  
  P1 = rowSums(Z_abs < real_abs - expected.diff)
  P2 = rowSums(Z_abs > real_abs + expected.diff)
  Orig$pval_abs = (P1 + P2) / n
  
  
  if(method_resampling == "Reshuffle"){
    Orig$pval_sign = 1- Orig$pval_sign
    Orig$pval_abs = 1- Orig$pval_abs
  }
  
  Orig$Padj_sign = (stats::p.adjust(Orig$pval_sign, method = pvalmethod))
  Orig$Padj_abs = (stats::p.adjust(Orig$pval_abs, method = pvalmethod))

  
  ## Running the correlation
  if( savecor == T){
    Total_Correlation = as.data.frame(stats::cor(t(dfExpression), method = method))
    Total_Correlation = wTO.in.line(Total_Correlation)
    names(Total_Correlation) = c("Node.1", "Node.2", "Cor")
  }
  if( savecor == F){
    Total_Correlation = NULL
  }
  Cutoffs = Cut.off(ALL_DT_sig, "wTO - Resampling")
  Cutoffs_abs = Cut.off(ALL_DT_abs, "|wTO| - Resampling")
  
  
  Orig = Orig[, -"Rep"]
  
  
  Orig$wTO_abs = as.numeric(Orig$wTO_abs)
  Orig$wTO_sign = as.numeric(Orig$wTO_sig)
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
  
  # par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(3,1))
  graphics::plot(as.matrix(WTO_sig), as.matrix(WTO_abs), axes = F,
                 xlab = "|wTO|", ylab = "wTO",
                 main = "|wTO| vs wTO", pch = 16, xlim = c(-1,1), ylim = c(0,1),
                 col.main = "steelblue2", col.lab = "steelblue2")
  graphics::axis(1, las = 1, cex.axis = 0.6, col = "steelblue",
                 
                 col.ticks = "steelblue3", col.axis = "steelblue")
  graphics::axis(2, las = 1, cex.axis = 0.6, col = "steelblue",col.ticks = "steelblue3", col.axis = "steelblue")
  
  graphics::par(xpd=FALSE)
  graphics::abline( h = 0,  lty = 2, col = "gray50")
  graphics::abline(v = 0,  lty = 2, col = "gray50")
  
  
  graphics::plot(Orig$wTO_sign, Orig$pval_sign, axes = F,
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
  return(output)
}
