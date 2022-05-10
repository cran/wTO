#' @title wTO.fast
#' @author Deisy Morselli Gysi <deisy at bioinf.uni-leipzig.de>

#' @param n Number of resamplings, used to compute the empirical distribuitions of the links. Default is set to 100.
#' @param Data data.frame containing the count / expression data for the correlation.
#' @param Overlap Set of nodes of interest, where the Overlapping weights will be computed.
#' @param method Type of the correlation that should be used. "s" / "spearman" will compute the rank spearman correlation, "p" / "pearson" will compute the linear correlation. If no value is given, the default is to use "p".
#' @param sign Should the wTO be signed?
#' @param delta expected difference between the real wTO and the bootstraped.
#' @param method_resampling method of the resampling. Bootstrap or BlockBootstrap.If the second is used, please give the lag (time dependency among the data).
#' @param lag Time dependency for the blocked bootstrap (for time series).
#' @param ID ID of the samples for the blocked bootstrap (for repeated measures).
#' 
#' @description Compute the wTO and also the bootstraps. Proposed at arXiv:1711.04702. This is a quicker version of the wTO.Complete. It doesn't contain diagnose plots nor a parallel version.
#' @importFrom  parallel makeCluster clusterExport clusterApplyLB  stopCluster
#' @importFrom  data.table rbindlist dcast
#' @importFrom  som normalize
#' @importFrom  stats cor p.adjust reshape pchisq
#' @importFrom  graphics plot axis par abline legend
#' @import  magrittr 
#' @export
#' @examples 
#' # wTO.fast(Data = Microarray_Expression1,
#' # Overlap = ExampleGRF$x, 
#' # method = "p")
#' 
#' # For a time series with lag = 4
#' # wTO.fast(Data = Microarray_Expression1,
#' # Overlap = ExampleGRF$x, 
#' # method = "p", 
#' # method_resampling = 'BlockBootstrap', 
#' # lag = 4)
#' 
#' # For a study where the individuals were measured multiple times.
#' # wTO.fast(Data = Microarray_Expression1,
#' # Overlap = ExampleGRF$x, 
#' # method = "p", 
#' # method_resampling = 'BlockBootstrap', 
#' # ID = rep(1:9, each= 2))

wTO.fast = function(Data, 
                    Overlap = row.names(Data), 
                    method = 'p', sign = 'sign', 
                    delta = 0.2, n = 10,
                    method_resampling = 'Bootstrap', lag = NULL, ID = NULL){
  Overlap = unique(as.character(Overlap))
  `%ni%` <- Negate(`%in%`)
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
  
  if(method %ni% c("s", "spearman", "p", "pearson")){
    stop('Method must be: "s", "spearman", "p" or "pearson".')
  }
  
  if(method_resampling %ni% c("Bootstrap", "BlockBootstrap")){
    stop('Method must be: "Bootstrap" or "BlockBootstrap".')
  }
  if(method_resampling %in% "BlockBootstrap"){
    if (is.null(lag)&is.null(ID)){
      stop('If you want to use the "BlockBootstrap" please give a lag or the indivuals ID.')
    }
    if(!is.null(lag)&!is.null(ID)){
      stop('If you want to use the "BlockBootstrap" please give a lag OR the indivuals ID.')
    }
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
  
  wtomelt0 =  wTO::CorrelationOverlap(Data = Data,
                                      Overlap = Overlap, 
                                      method = method) %>%
    wTO::wTO(., sign)
  `%>%` <- magrittr::`%>%`
  . <- NULL
  for ( i in 1:n){
    message(' ',i,' ', appendLF = FALSE)
    
    if(method_resampling == 'BlockBootstrap'){
      if (!is.null(lag)){
        nsampl = ifelse (ncol(Data) %% lag == 0, ncol(Data) %/% lag, ncol(Data) %/% lag +1)
        Y = sample(1:nsampl, size = nsampl, replace =  T)
        Vect = Y*lag
        j = lag - 1
        while( j > 0){
          Vect = cbind(Vect , Y*lag - j)
          j = j - 1
        }
        
        SAMPLES = c(Vect)
        SAMPLES[SAMPLES > ncol(Data)] <- NA
        SAMPLE = stats::na.exclude(SAMPLES)
        Data_boot = Data[,SAMPLE]
      }
      
      if(!is.null(ID)){
        ID %<>% as.factor
        bootID = sample(levels(ID), replace = TRUE)
        
        
        Data_boot = subset(Data, select = ID %in% bootID[1])
        
        for (k in 2:length(bootID)){
          Data_boot = cbind(Data_boot,
                            subset(Data, select = ID %in% bootID[k]))
        }
      }
      
      res = wTO::CorrelationOverlap(Data = Data_boot, Overlap = Overlap, method = method) %>% 
        wTO::wTO(., sign)
      
    }
    else if (method_resampling != 'BlockBootstrap'){
      res = wTO::CorrelationOverlap(Data = Data[,sample(1:ncol(Data), replace  = TRUE)], Overlap = Overlap, method = method) %>% 
        wTO::wTO(., sign)
    }
    
    
    U  = (res < wtomelt0 - delta) + (res > wtomelt0 + delta)
    if ( i == 1){
      out = U}
    if (i != 1){
      out = out + U
    }
    rm(res)
    rm (U)
  }
  
  wtomelt0 = wTO.in.line(wtomelt0)
  cor      = wTO.in.line(out)
  adj.pval = p.adjust(cor$wTO/n, method = 'BH')
  
  pval = data.table::data.table(wtomelt0, pval = cor$wTO/n, pval.adj = adj.pval)
  
  message('Done!')
  return(pval)
}
