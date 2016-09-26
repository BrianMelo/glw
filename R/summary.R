#------------------------------------------------------------------------------------------------------------------#
#------------HPD Intervals-----------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#
#' HPD Intervals
#' 
#' Computes the High Posterior Density Intervals for each column of a given matrix - requires TeachingDemos package
#' 
#' @param parm: matrix where each column is a sample from a parameter posterior distribution
#' @param conf: credibility level
#' 
#' @return A matrix where each row is the HPD of each column of the matrix \code{parm}
#' 
#' @seealso \code{\link[TeachingDemos]{hpd}} and \code{\link[TeachingDemos]{emp.hpd}} in \code{TeachingDemos} package.
#' 
#' @export
colHpds <- function(parm, conf=0.95) {
  ic <- matrix(0, nrow=ncol(parm), ncol=2)
  for(i in 1:ncol(parm)) ic[i,] <- TeachingDemos::emp.hpd(parm[,i], conf=conf)
  return(ic)
}