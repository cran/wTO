#' @title wTO.Consensus
#' @aliases wTO.Consensus
#' @description Consensus requires a list of data.frame containing the pair of nodes, and the wTO values for all networks that need to be joined. Reference: 	arXiv:1711.04702
#' @param data list of data.frame containing the "Node.1", "Node.2" and "wTO".
#' @author Deisy Morselli Gysi <deisy at bioinf.uni-leipzig.de>

#' @export
#' @importFrom plyr join_all
#' @importFrom stats pchisq

#' @examples
#' \dontrun{
#'EXAMPLE =  wTO.Complete( k =1, n = 200, Data = Microarray_Expression1,
#'                                       Overlap = ExampleGRF$x, method = "p")
#'
#' # Constructing the consensus network
#' data = list(data.frame(Node.1 = EXAMPLE$wTO$Node.1,
#' Node.2 = EXAMPLE$wTO$Node.2,
#' wto_sig = EXAMPLE$wTO$wTO_sign,
#' pvalsig = EXAMPLE$wTO$pval_sig),
#' data.frame(Node.1 = EXAMPLE$wTO$Node.1,
#' Node.2 = EXAMPLE$wTO$Node.2,
#' wtoabs = EXAMPLE$wTO$wTO_abs,
#' pvalabs = EXAMPLE$wTO$pval_abs) )
#' CONS = wTO.Consensus(data)
#' 
#' }


wTO.Consensus = function(data){
  if (!is.list(data)){
    stop("data must be a list of data.frames.")
  }
  ### Weight
  weight = pval = nodes = list()
  
  for ( i in 1:length(data)){
    weight[[i]] = data[[i]][,1:3]
    pval[[i]] = data[[i]][,c(1:2,4)]
    ID = unique(c(levels(data[[i]]$Node.1), levels(data[[i]]$Node.2)))
    nodes[[i]] = data.frame(ID =  ID)
    names(weight[[i]])[3] = paste0(names(weight[[i]][3]), i)
    names(pval[[i]])[3] = paste0(names(pval[[i]][3]), i)
  }
  
  weight = plyr::join_all(weight, type = 'full')
  pval = plyr::join_all(pval, type = 'full')
  nodes = plyr::join_all(nodes, type = 'inner')

  message(paste('Total common nodes:', nrow(nodes)))
  weight = subset(weight, weight$Node.1 %in% nodes$ID & weight$Node.2 %in% nodes$ID)
  pval = subset(pval, pval$Node.1 %in% nodes$ID & pval$Node.2 %in% nodes$ID)
  pval[is.na(pval)] <- 1
  weight[is.na(weight)] <- 0.01
  
  wTOCN = CN_aux(weight[, -c(1:2)])
  pvalue_fisher = fishermethod(pval[, -c(1:2)])
  
  Out = data.frame(Node.1 = pval[,1], Node.2 = pval[,2],
                   CN = wTOCN, pval.fisher = pvalue_fisher)
  return(Out)
}


fishermethod = function(data_x){
  chi = rowSums(log(data_x))*-2
  pval = sapply(chi, function(x) stats::pchisq(x, 2*ncol(data_x), lower.tail = FALSE))
  return(pval)
}

CN_aux = function(data_x){
  abs_x = apply(data_x, 2, abs)
  sum_abs_x = apply(abs_x, 1, sum)
  div = (abs_x/sum_abs_x) * data_x
  wTO_cons = apply(div, 1, sum)
  return(wTO_cons)
}