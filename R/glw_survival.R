#------------------------------------------------------------------------------------------------------------------------#
#---------------Estimated survival functions og the GLW finite mixture model---------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------#
#'Estimated Survival Function of GLW Mixture
#'
#'\code{urv_glw} estimates the survival function of the GLW finite mixture regression model for a 
#' observation with vector of covariates \code{X} using a sample from the posterior distribution 
#' of the GLW parameters.
#'
#' @param parm: resulting sample from MCMC simulation.
#' @param y: response variable (survival time).
#' @param X: vector of covariates for one observation.
#' @param mixture: indicates the finite mixture to be fitted to the data. Options are: 'glw', 'gl', 'gw', 'lw'.
#' @param by: increment of time.
#' 
#' @return A list with two obects: \code{surv} and \code{time}. \code{$surv} contains the estimated survival 
#' probabilities and \code{$time} the respective time. If \code{by}=NULL then time is equal to y.
#' 
#' @export
surv_glw <- function(parm, y, X, mixture, by=NULL) {
  out <- eval(parse(text=paste('S.', mixture, '(parm=parm,y=y, X=X, by=by)', sep='')))
  return(out)
}     

S.lw <- function(parm, y, X, by) {
  
  if(!is.null(by)) {tempo <- seq(min(y), max(y), by=by)} else{tempo <- y}
  
  s <- matrix(nrow=nrow(parm), ncol=length(tempo)) 
  J <- ncol(parm)-3
  r <- length(X)
  
  Betas <- as.matrix(parm[,1:J])
  Sigma <- parm[,(J+1)]
  Pesos <- parm[,(J+2):(J+3)]
  
  xbeta <- 0
  if( J>r ) {
    if(r>1) xbeta <- t(X[2:r]%*%t(Betas[,3:J]))
    mu.lnorm <- exp(Betas[,1] + xbeta)
    mu.weibu <- exp(Betas[,2] + xbeta)
    mu.gamma <- mu.lnorm
  }
  if(J==r) {
    if(r>1) xbeta <- X[2:r]%*%t(Betas[,2:J])
    mu.gamma <- exp(Betas[,1] + xbeta)
    mu.lnorm <- mu.gamma
    mu.weibu <- mu.gamma
  }
  
  for(i in 1:nrow(parm)){
    parms <- Parm3(mu.gamma[i], mu.lnorm[i], mu.weibu[i], Sigma[i])
    s[i,] <- Pesos[i,1]*plnorm(tempo, parms$mu_l, parms$sig_l, lower.tail=F) + 
      Pesos[i,2]*pweibull(tempo, shape=parms$gama, scale=parms$delta, lower.tail=F)
    
  }
  out <- list(surv=colMeans(s), time=tempo)
  return(out)
}

S.gw <- function(parm, y, X, by) {
  
  if(!is.null(by)) {tempo <- seq(min(y), max(y), by=by)} else{tempo <- y}
  
  s <- matrix(nrow=nrow(parm), ncol=length(tempo)) 
  J <- ncol(parm)-3
  r <- length(X)
  
  Betas <- as.matrix(parm[,1:J])
  Sigma <- parm[,(J+1)]
  Pesos <- parm[,(J+2):(J+3)]
  
  xbeta <- 0
  if( J>r ) {
    if(r>1) xbeta <- t(X[2:r]%*%t(Betas[,3:J]))
    mu.gamma <- exp(Betas[,1] + xbeta)
    mu.lnorm <- mu.gamma
    mu.weibu <- exp(Betas[,2] + xbeta)
  }
  if(J==r) {
    if(r>1) xbeta <- X[2:r]%*%t(Betas[,2:J])
    mu.gamma <- exp(Betas[,1] + xbeta)
    mu.lnorm <- mu.gamma
    mu.weibu <- mu.gamma
  }
  
  for(i in 1:nrow(parm)){
    parms <- Parm3(mu.gamma[i], mu.lnorm[i], mu.weibu[i], Sigma[i])
    s[i,] <- Pesos[i,1]*pgamma(tempo, shape=parms$alpha, scale=parms$beta, lower.tail=F) + 
      Pesos[i,2]*pweibull(tempo, shape=parms$gama, scale=parms$delta, lower.tail=F)
    
  }
  out <- list(surv=colMeans(s), time=tempo)
  return(out)
}

S.gl <- function(parm, y, X, by) {
  
  if(!is.null(by)) {tempo <- seq(min(y), max(y), by=by)} else{tempo <- y}
  
  s <- matrix(nrow=nrow(parm), ncol=length(tempo)) 
  J <- ncol(parm)-3
  r <- length(X)
  
  Betas <- as.matrix(parm[,1:J])
  Sigma <- parm[,(J+1)]
  Pesos <- parm[,(J+2):(J+3)]
  
  xbeta <- 0
  if( J>r ) {
    if(r>1) xbeta <- t(X[2:r]%*%t(Betas[,3:J]))
    mu.gamma <- exp(Betas[,1] + xbeta)
    mu.lnorm <- exp(Betas[,2] + xbeta)
    mu.weibu <- mu.gamma
  }
  if(J==r) {
    if(r>1) xbeta <- X[2:r]%*%t(Betas[,2:J])
    mu.gamma <- exp(Betas[,1] + xbeta)
    mu.lnorm <- mu.gamma
    mu.weibu <- mu.gamma
  }
  
  for(i in 1:nrow(parm)){
    parms <- Parm3(mu.gamma[i], mu.lnorm[i], mu.weibu[i], Sigma[i])
    s[i,] <- Pesos[i,1]*pgamma(tempo, shape=parms$alpha, scale=parms$beta, lower.tail=F) + 
      Pesos[i,2]*plnorm(tempo, parms$mu_l, parms$sig_l, lower.tail=F)
    
  }
  out <- list(surv=colMeans(s), time=tempo)
  return(out)
}