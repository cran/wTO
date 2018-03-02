#' @title wTO.Consensus
#' @aliases wTO.Consensus
#' @description Consensus requires a list of data.frame containing the pair of nodes, and the wTO values for all networks that need to be joined. Reference: 	arXiv:1711.04702
#' @param data list of data.frame containing the "Node.1", "Node.2" and "wTO".
#' @param full Missing links should be considered zero?
#' @author Deisy Morselli Gysi <deisy at bioinf.uni-leipzig.de>

#' @export
#' @importFrom plyr join_all

#' @examples
#' \dontrun{
#'EXAMPLE =  wTO.Complete( k =1, n = 200, Data = ExampledfExpression,
#'                                       Overlap = ExampleGRF$x, method = "p")
#'
#' #Selection of only the significative ones for the Consensus
#' Ex_k1_cor_p_boot_p005_sig = subset(EXAMPLE$wTO,
#' EXAMPLE$wTO$pval_sig < 0.05,
#' select = c("Node.1", "Node.2", "wTO_sign"))
#' Ex_k1_cor_p_boot_p005_abs = subset(EXAMPLE$wTO,
#' EXAMPLE$wTO$pval_abs < 0.05,
#' select = c("Node.1", "Node.2", "wTO_abs"))
#' # Constructing the consensus network
#' CONS = wTO.Consensus(data = list(Ex_k1_cor_p_boot_p005_sig,
#' Ex_k1_cor_p_boot_p005_abs))
#' }


wTO.Consensus = function(data, full = FALSE){
  if (!is.list(data)){
    stop("data must be a list of data.frames.")
  }

  Nodes = list()
  for ( i in 1:length(data)){
    names(data[[i]])[1:2] = c('Node.1', 'Node.2')
    names(data[[i]])[3]= paste('wTO', i, sep = '_')
    Nodes[[i]] = data.frame(ID = unique(c(as.character(data[[i]]$Node.1),
                                          as.character(data[[i]]$Node.2))))
  }
  Nodes = plyr::join_all(Nodes, type = 'inner')
  if (full == FALSE){
    data_xx = plyr::join_all(data, by = c("Node.1", "Node.2"),
                             type = "inner")
  }
  else if(full == TRUE){
    data_xx = plyr::join_all(data, by = c("Node.1", "Node.2"),
                             type = "full")
  }
  data_xx = subset(data_xx, data_xx$Node.1 %in% Nodes$ID & data_xx$Node.2 %in% Nodes$ID)
  wTO_NA = apply(data_xx[, -c(1, 2)], 1, function(x) sum(is.na(x)))
  wTO_NA = wTO_NA/ ncol(data_xx[, -c(1, 2)])
  data_xx[is.na(data_xx)] <- 0
  data_x = data_xx[, -c(1, 2)]
  abs_x = apply(data_x, 2, abs)
  sum_abs_x = apply(abs_x, 1, sum)
  div = (abs_x/sum_abs_x) * data_x
  wTO_cons = apply(div, 1, sum)

  cons_wto = data.table::data.table(data.frame(Node.1 = data_xx$Node.1, Node.2 = data_xx$Node.2,
                        wTO_Cons = wTO_cons), prop.NA = wTO_NA)
  class(cons_wto)<- append('wTOCN', class(cons_wto))
  return(cons_wto)
}
