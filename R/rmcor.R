
#' @title wTO.aux.each
#' @description Computes the repeated measures correlation. based on Blendman Altman (1995). Implemented by <https://doi.org/10.3389/fpsyg.2017.00456>
#' @keywords internal


rmcor <- function (ID, 
                   Measure1, 
                   Measure2)
                                                                
{
  newdat <- stats::na.omit(data.frame(ID, Measure1, 
                                      Measure2))
  Participant <- newdat$ID
  Measure1 <- newdat$Measure1
  Measure2 <- newdat$Measure2
  
  lmmodel <- stats::lm(Measure2 ~ Participant + Measure1)
  lmslope <- stats::coef(lmmodel)["Measure1"]
  errordf <- lmmodel$df.residual
  corrsign <- sign(lmslope)
  type3rmcorr <- stats::drop1(lmmodel, ~., test = "F")
  SSFactor <- type3rmcorr$"Sum of Sq"[3]
  SSresidual <- type3rmcorr$RSS[1]
  rmcorrvalue <- as.numeric(corrsign * sqrt(SSFactor/(SSFactor + 
                                                        SSresidual)))
  return(rmcorrvalue)
}
