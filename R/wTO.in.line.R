#' @title wTO.in.line
#'
#' @param d correlation matrix to be converted into the line format.
#' @description Transforms a correlation matrix into the line format.
#' @return the wTO matrix into a data.frame: Node1, Node2 and wTO.
#' @export
#'

wTO.in.line <-function(d){
  names2= matrix(names(d), nrow = nrow(d), ncol = ncol(d), byrow = T)
  names3= matrix(row.names(d), nrow = nrow(d), ncol = ncol(d))

  Genes.1 <- names2[upper.tri(names2)]
  Genes.2 <- names3[upper.tri(names3)]
  M.Genes.1 <- apply(cbind(Genes.1, Genes.2), 1, min)
  M.Genes.2<- apply(cbind(Genes.1, Genes.2), 1, max)
  # M.nomes = paste(M.Genes.1, M.Genes.2, sep = "~")


  M.sup <- d[upper.tri(d)]
  corre=as.data.frame(cbind(M.Genes.1 ,M.Genes.2, M.sup))
  names(corre)<-c("Node.1", "Node.2", "wTO")
  row.names(corre)= paste(M.Genes.1, M.Genes.2, sep = "<->")
  corre$wTO = as.numeric(as.matrix(corre$wTO))

  return(corre)
}
