#---------------------------------------------------------------------#
#------CPO - Conditional Predictive Ordinate--------------------------#
#---------------------------------------------------------------------#
#'CPO of the GLW Finite Mixture Model
#'
#' \code{cpo_glw} estimates the CPO based on a MCMC sample of the posterior distribuiton.
#' 
#' @param parm: sample of posterior parameters.
#' @param Data: the same data as in \code{\link{glwfm}}.
#' @param mixture: indicates the finite mixture to be fitted to the data. Options are: 'glw', 'gl', 'gw', 'lw', 'g', 'l', 'w'.
#' 
#' @return Returns the Conditional Predictive Ordinate of the GLW finite mixture model for every observation.
#' 
#' @export
cpo_glw <- function(parm, Data, mixture) {
  if(nchar(mixture)>1){
   cpo <- matrix(nrow=nrow(parm), ncol=nrow(Data$X)) 
   for(i in 1:nrow(parm)) cpo[i,] <- 1/exp(loglik(parm=parm[i,], Data=Data, mixture))
   out <- 1/colMeans(cpo)
  } else {
   if(mixture=='g') {out <- gamma.cpo(parm=parm, Data=Data)}
   if(mixture=='l') {out <- lnorm.cpo(parm=parm, Data=Data)}
   if(mixture=='w') {out <- weibull.cpo(parm=parm, Data=Data)}
   }
  return(out)
  }

#---------------------------------------------------------------------#
#------Gamma----------------------------------------------------------#
#---------------------------------------------------------------------#
gamma.cpo <- function(parm, Data) {
  cpo <- matrix(nrow=nrow(parm), ncol=nrow(Data$X)) 
  for(i in 1:nrow(parm)) cpo[i,] <- 1/exp(gamma.loglik(parm=parm[i,], Data=Data))
  out <- 1/colMeans(cpo)
  return(out)
}

#---------------------------------------------------------------------#
#------Weibull--------------------------------------------------------#
#---------------------------------------------------------------------#
weibull.cpo <- function(parm, Data) {
  cpo <- matrix(nrow=nrow(parm), ncol=nrow(Data$X)) 
  for(i in 1:nrow(parm)) cpo[i,] <- 1/exp(weibull.loglik(parm=parm[i,], Data=Data))
  out <- 1/colMeans(cpo)
  return(out)
}

#---------------------------------------------------------------------#
#------LogNormal------------------------------------------------------#
#---------------------------------------------------------------------#
lnorm.cpo <- function(parm, Data) {
  cpo <- matrix(nrow=nrow(parm), ncol=nrow(Data$X)) 
  for(i in 1:nrow(parm)) cpo[i,] <- 1/exp(lnorm.loglik(parm=parm[i,], Data=Data))
  out <- 1/colMeans(cpo)
  return(out)
}

