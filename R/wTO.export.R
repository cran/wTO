#' @title wTO.export
#' @aliases wTO.export
#' @description Exports the significative interactions, the wTO weight and pvalues into a .txt file, tab separeted. This file can be imported in other visualization tools (Cytoscape for example).

#' @param DATA Output from the function wTO.Complete or wTO.Consensus.
#' @param path Path and file name where the .txt file should be saved.
#' @param sign Should the network contain the results for the signed network or unsigned? Only for data coming from wTO.Complete.
#' @param pvalue cutoff p-value for the network. Only for data coming from wTO.Complete.
#' @param padj cutoff adjusted p-value for the network. Only for data coming from wTO.Complete.
#' @param prop.NA cutoff proportion of NAs for the network. Only for data coming from wTO.Consensus.
#' @importFrom utils write.table
#'
#' @export
#'
#' @examples
#' \dontrun{
#' EXAMPLE =  wTO.Complete( k =1, n = 200, Data = Microarray_Expression2,
#'                                       Overlap = ExampleGRF$x, method = "p")
#' wTO.export(EXAMPLE , './EXAMPLE.txt')
#'
#' #Selection of only the significative ones for the Consensus
#' Ex_k1_cor_p_boot_p005_sig = subset(EXAMPLE$wTO,
#' EXAMPLE$wTO$pval_sig < 0.05,
#' select = c("Node.1", "Node.2", "wTO_sign"))
#' Ex_k1_cor_p_boot_p005_abs = subset(EXAMPLE$wTO,
#' EXAMPLE$wTO$pval_abs < 0.05,
#' select = c("Node.1", "Node.2", "wTO_abs"))
#' # Constructing the consensus network
#' CN = wTO.Consensus(data = list(Ex_k1_cor_p_boot_p005_sig,
#' Ex_k1_cor_p_boot_p005_abs))
#' wTO.export(CN, './CN.txt')
#' ### You can store the result on the workspace.
#' y = wTO.export(CN, './CN.txt')
#' head(y)
#' }

#'
#'
#'
wTO.export = function(DATA, path, sign = TRUE, pvalue = 0.05, padj = 0.05, prop.NA = 0.5){
  if(any(class(DATA) %in% 'wTO')){
    if(sign == TRUE){
      save = subset(DATA$wTO, DATA$wTO$pval_sig < pvalue & DATA$wTO$Padj_sig < padj,
                    select = c("Node.1" ,  "Node.2" ,  "wTO_sign" , "pval_sig" , "Padj_sig"))
    }
    if(sign == FALSE){
      save = subset(DATA$wTO, DATA$wTO$pval_abs < pvalue & DATA$wTO$Padj_abs,
                    select = c("Node.1" ,  "Node.2" ,  "wTO_abs" , "pval_abs" , "Padj_abs"))
    }
    names(save) = c(c("Node.1" ,  "Node.2" ,  "wTO" , "pval" , "Padj"))
  }
  if(any(class(DATA) %in% 'wTOCN')){
    allowedNA = prop.NA
      save = subset(DATA, DATA$prop.NA < allowedNA,
                    select = c("Node.1" ,  "Node.2" ,  "wTO_Cons" , "prop.NA" ))


    names(save) = c("Node.1" ,  "Node.2" ,  "wTO" , "prop.NA")
  }

  write.table(save, path, quote = F, sep = '\t', row.names = FALSE)
  return(invisible(save))
}
