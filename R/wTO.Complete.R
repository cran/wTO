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
#' @importFrom  data.table rbindlist
#' @importFrom  som normalize
#' @importFrom  stats cor p.adjust
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
#' wTO.Complete( k =2, n = 10, Data = ExampledfExpression,
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
    message(paste('There are',DIM_GRF, "overlapping nodes." ))
  }


  message("This function might take a long time to run. Don't turn off the computer.")

  graphics::par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(3,2))

  ## For the original data
  # real_Genes = dfExpression
  Saving = wTO::Correlation.Overlap(Data = dfExpression, Overlap = GRF, method = method)
  WTO_abs = wTO::wTO(A = Saving[[2]],  sign = "abs")
  WTO_sig = wTO::wTO(A = Saving[[2]],  sign = "sign")
  Cor_real = wTO.in.line(WTO_sig)
  Cor_real_abs = wTO.in.line(WTO_abs)
  names(Cor_real) = names(Cor_real_abs) <-c("Node.1", "Node.2", "wTO_0")
  SAVE_VALUES = Cor_real
  SAVE_VALUES_ABS = Cor_real_abs



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
    assign("Correlation.Overlap", wTO::Correlation.Overlap, envir = WTO)
    assign("wTO", wTO::wTO, envir = WTO)
    assign("wTO.aux.each", wTO.aux.each, envir = WTO)
    assign("method_resampling", method_resampling, envir = WTO)
    assign("sample_ind", sample_ind, envir = WTO)
    assign("lag", lag, envir = WTO)
    cl = parallel::makeCluster(k)
    parallel::clusterExport(cl, "dfExpression", envir = WTO)
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

  ALL  = data.frame(data.table::rbindlist(OUTPUT, use.names = T))
  SAVE_VALUES = data.frame(SAVE_VALUES, ALL$Cor_star)
  SAVE_VALUES_ABS = data.frame(SAVE_VALUES_ABS, ALL$Cor_star_abs)

  Z = as.matrix(SAVE_VALUES[, - c(1:3)])

  Z_ABS = as.matrix(SAVE_VALUES_ABS[, - c(1:3)])



  P1 = rowSums((Z) <= SAVE_VALUES[,3]- expected.diff)
  P2 = rowSums((Z) >= SAVE_VALUES[,3]+ expected.diff)
  Cor_real$P = (P1 + P2) / ncol(Z)

  if(method_resampling == "Reshuffle"){
    Cor_real$P = 1- Cor_real$P
  }

  Cor_real$Padj = (stats::p.adjust(Cor_real$P, method = pvalmethod))


  P1 = rowSums((Z) <= SAVE_VALUES_ABS[,3]- expected.diff)
  P2 = rowSums((Z) >= SAVE_VALUES_ABS[,3]+ expected.diff)
  Cor_real_abs$P = (P1 + P2) / ncol(Z_ABS)
  if(method_resampling == "Reshuffle"){
    Cor_real_abs$P = 1- Cor_real_abs$P
  }
  Cor_real_abs$Padj = (stats::p.adjust(Cor_real_abs$P, method = pvalmethod))

  ## Running the correlation
  if( savecor == T){
    Total_Correlation = as.data.frame(stats::cor(t(dfExpression), method = method))
    Total_Correlation = wTO.in.line(Total_Correlation)
    names(Total_Correlation) = c("Node.1", "Node.2", "Cor")
  }
  if( savecor == F){
    Total_Correlation = NULL
  }
  Cutoffs = Cut.off(SAVE_VALUES, "wTO - Resampling")
  Cutoffs_abs = Cut.off(SAVE_VALUES_ABS, "|wTO| - Resampling")

  names(Cor_real) = names(Cor_real_abs) = c("Node.1", "Node.2", "wTO", "pval", "pval.adj")
  output = list(wTO = Cor_real,
                abs.wTO = Cor_real_abs,
                # SAVE = SAVE_VALUES,
                # SAVE2 = SAVE_VALUES_ABS,
                Correlation = Total_Correlation,
                Empirical.Quantile= Cutoffs$Empirical.Quantile,
                Quantile = Cutoffs$Quantile ,
                Empirical.Quantile.abs= Cutoffs_abs$Empirical.Quantile,
                Quantile.abs = Cutoffs_abs$Quantile
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


  graphics::plot(output$wTO$wTO, output$wTO$pval, axes = F,
       xlab = "wTO", ylab = "p-value", ylim = c(0,1), xlim = c(-1,1), col.main = "steelblue2", col.lab = "steelblue2",
       main = "wTO vs p-value",
       pch = 16)
  graphics::axis(1, las = 1, cex.axis = 0.6, col = "steelblue",

                 col.ticks = "steelblue3", col.axis = "steelblue")
  graphics::axis(2, las = 1, cex.axis = 0.6, col = "steelblue",col.ticks = "steelblue3", col.axis = "steelblue")

  graphics::par(xpd=FALSE)
  graphics::abline( v = output$Empirical.Quantile,  lty = 2, col = c("red", "orange", "yellow", "yellow", "orange", "red"))
  graphics::par(xpd=T)
  graphics::legend(c(0.9,0), c ("Empirical Quantiles", '0.1%','2.5%','10%','90%','97.5%','99.9%'),
         col = c("white", "red", "orange", "yellow", "yellow", "orange", "red"), lwd = 2, bty = "n",
         inset=c(-0.8,0), cex = 0.5 )

  graphics::par(xpd=FALSE)
  graphics::plot(output$abs.wTO$wTO, output$abs.wTO$pval, axes = F,
       xlab = "|wTO|", ylab = "p-value", ylim = c(0,1), xlim = c(0,1),
       main = "|wTO| vs p-value",
       pch = 16, col.main = "steelblue2", col.lab = "steelblue2")
  graphics::axis(1, las = 1, cex.axis = 0.6, col = "steelblue",

                 col.ticks = "steelblue3", col.axis = "steelblue")
  graphics::axis(2, las = 1, cex.axis = 0.6, col = "steelblue",col.ticks = "steelblue3", col.axis = "steelblue")


  graphics::abline( v = output$Empirical.Quantile.abs,  lty = 2, col = c("red", "orange", "yellow", "yellow", "orange", "red"))
  graphics::par(xpd=T)
  graphics::legend(c(0.9,0), c ("Empirical Quantiles", '0.1%','2.5%','10%','90%','97.5%','99.9%'),
         col = c("white", "red", "orange", "yellow", "yellow", "orange", "red"), lwd = 2, bty = "n",
         inset=c(-0.8,0), cex = 0.5 )
  return(output)
}
