#' @title wTO.rep_measure
#' @author Deisy Morselli Gysi <deisy at bioinf.uni-leipzig.de>

#' @param n Number of resamplings, used to compute the empirical distribuitions of the links. Default is set to 100.
#' @param Data data.frame containing the count / expression data for the correlation.
#' @param Overlap Set of nodes of interest, where the Overlapping weights will be computed.
#' @param sign Should the wTO be signed?
#' @param delta expected difference between the real wTO and the bootstraped.
#' @param ID a vector with the individuals identification

#' @description Compute the wTO for a repeated measures expermiment and also the bootstraps. Proposed at arXiv:1711.04702. This is a quicker version of the wTO.Complete. It doesn'T contain diagnose plots nor a parallel version.
#' @importFrom  parallel makeCluster clusterExport clusterApplyLB  stopCluster
#' @importFrom  data.table rbindlist dcast
#' @importFrom  som normalize
#' @importFrom  stats cor p.adjust reshape pchisq
#' @importFrom  graphics plot axis par abline legend
#' @import  magrittr 
#' @export
#' @examples 
#' 
#' #wTO.rep_measure(Data = Microarray_Expression1, ID = rep(c(1:9),2), 
#' #Overlap = ExampleGRF$x)


wTO.rep_measure = function(Data, Overlap = row.names(Data), ID,
                           sign = 'sign', 
                           delta = 0.2, n = 10){
  Overlap = unique(as.character(Overlap))
  `%ni%` <- Negate(`%in%`)
  ID = as.factor(ID)
  ##### Messages
  
  if(is.numeric(n) == F){
    stop("n must be numeric.")
  }
  if(n <= 0){
    stop("n must be greater than 0.")
  }
  if(is.data.frame(Data) == F){
    stop("Data must be a data.frame.")
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
  
  Datat = t(Data)
  Cor = matrix(0, nrow = ncol(Datat), ncol = ncol(Datat)) %>% data.frame()
  names(Cor) = row.names(Cor)= names(Datat)
  message('Starting correlations.')
  for( i in 1:ncol(Datat)){
    for(j in i:(ncol(Datat) )){
      if( i == j){
        Cor[i,j] = 0
      }
      else{
        Cor[i,j] = Cor[j,i] = suppressWarnings(rmcor(ID,Datat[,i],Datat[,j]))
      }
    }
  }
  Cor[is.na(Cor)] = 0
  names(Cor) = row.names(Cor)= colnames(Datat)
  wtomelt0 = subset(Cor, row.names(Cor) %in% Overlap) %>% wTO::wTO(., sign)
  
  `%>%` <- magrittr::`%>%`
  . <- NULL
  for ( B in 1:n){
    message('.', appendLF = FALSE)
    
    bootID = sample(levels(ID), replace = TRUE)
    
    Data_boot = subset(Datat, ID == bootID[1])
    for (k in 2:length(bootID)){
      Data_boot = rbind(Data_boot,
                        subset(Datat, ID == bootID[k]))
    }
    
    Cor = matrix(0, nrow = ncol(Data_boot), ncol = ncol(Data_boot)) %>% data.frame()
    names(Cor) = row.names(Cor)= colnames(Data_boot)
    for( i in 1:ncol(Data_boot)){
      for(j in 1:(ncol(Data_boot) )){
        Cor[i,j] = Cor[j,i] = suppressWarnings( rmcor(ID,Data_boot[,i],Data_boot[,j]))
      }
    }
    Cor[is.na(Cor)] = 0
    res = subset(Cor, row.names(Cor) %in% Overlap) %>% wTO::wTO(., sign)
    
    U  = (res < wtomelt0 - delta) + (res > wtomelt0 + delta)
    
    if ( B == 1){
      out = U}
    if (B != 1){
      out = out + U
    }
  }
  
  wtomelt0 = wTO.in.line(wtomelt0)
  cor      = wTO.in.line(out)
  pval = data.table::data.table(wtomelt0, pval = cor$wTO/n)
  message('Done!')
  return(pval)
}
