#' @title wTO.in.line
#'
#' @param d correlation matrix to be converted into the line format.
#' @description Transforms a correlation matrix into the line format.
#' @return the wTO matrix into a data.frame: Node1, Node2 and wTO.
#' @importFrom data.table as.data.table
#' @importFrom stats na.omit
#' @importFrom reshape2 melt
#' @author Deisy Morselli Gysi <deisy at bioinf.uni-leipzig.de>

#' @export
#'

wTO.in.line <-function(d){
  # names2= matrix(names(d), nrow = nrow(d), ncol = ncol(d), byrow = T)
  # names3= matrix(row.names(d), nrow = nrow(d), ncol = ncol(d))
  #
  # Genes.1 <- names2[upper.tri(names2)]
  # Genes.2 <- names3[upper.tri(names3)]
  # M.Genes.1 <- apply(cbind(Genes.1, Genes.2), 1, min)
  # M.Genes.2<- apply(cbind(Genes.1, Genes.2), 1, max)
  # # M.nomes = paste(M.Genes.1, M.Genes.2, sep = "~")
  #
  #
  # M.sup <- d[upper.tri(d)]
  # corre=data.table::as.data.table(cbind(M.Genes.1 ,M.Genes.2, M.sup))
  # names(corre)<-c("Node.1", "Node.2", "wTO")
  # # row.names(corre)= paste(M.Genes.1, M.Genes.2, sep = "<->")
  # # corre$wTO = as.numeric(as.matrix(corre$wTO))
  #

  # system.time(correlations<-cor(mydata,use="pairwise.complete.obs"))#get correlation matrix
  upperTriangle<-upper.tri(d, diag=F) #turn into a upper triangle
  d.upperTriangle<-d #take a copy of the original cor-mat
  d.upperTriangle[!upperTriangle]<-NA#set everything not in upper triangle o NA
  d_melted<-data.table::as.data.table(stats::na.omit(reshape2::melt(as.matrix(d.upperTriangle), value.name ="correlationCoef"))) #use melt to reshape the matrix into triplets, na.omit to get rid of the NA rows
  names(d_melted)<-c("Node.1", "Node.2", "wTO")


  return(d_melted)
}
