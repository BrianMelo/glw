#------------------------------------------------------------------------------------------------------------------------#
#---------------Estimated survival functions og the GLW finite mixture modelo--------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------#
#'Estimated Survival Function of GLW Promotion Time Cure Rate Model with covariates
#'
#'\code{S} estimates the survival function of the Standard Cure Rate Model, with covariates, using the GLW FM distribution for the susceptible individuals 
#' and with vector of covariates \code{X} using a sample from the posterior distribution 
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
Spop_ptcr <- function(parm, y, X, mixture, by=NULL) {
  out <- eval(parse(text=paste('Spop_ptcr.', mixture, '(parm=parm, y=y, X=X, by=by)', sep='')))
  return(out)
}

Spop_ptcr.glw <- function(parm, y, X, by) {
  
  if(!is.null(by)) {tempo <- seq(min(y), max(y), by=by)} else{tempo <- y}
  
  s <- matrix(nrow=nrow(parm), ncol=length(tempo)) 
  J <- ncol(parm)
  r <- length(X)
  
  Mu <- parm[,1]
  Sigma <- parm[,2]
  Pesos <- parm[,3:5]
  Betas <- as.matrix(parm[,6:J])
  
  xbeta <- X%*%t(Betas)
  eta <- exp(xbeta)
  
  for(i in 1:nrow(parm)){
    parms <- Par(Mu[i], Sigma[i])
    s[i,] <- exp(-eta[i]*(Pesos[i,1]*pgamma(tempo, shape=parms$alpha, scale=parms$beta) + 
                            Pesos[i,2]*plnorm(tempo, parms$mu_l, parms$sig_l) + 
                            Pesos[i,3]*pweibull(tempo, shape=parms$gama, scale=parms$delta)))
  }
  out <- list(surv=colMeans(s), time=tempo)
  return(out)
}

Spop_ptcr.gl <- function(parm, y, X, by) {
  
  if(!is.null(by)) {tempo <- seq(min(y), max(y), by=by)} else{tempo <- y}
  
  s <- matrix(nrow=nrow(parm), ncol=length(tempo)) 
  J <- ncol(parm)
  r <- length(X)
  
  Mu <- parm[,1]
  Sigma <- parm[,2]
  Pesos <- parm[,3:4]
  Betas <- as.matrix(parm[,5:J])
  
  xbeta <- X%*%t(Betas)
  eta <- exp(xbeta)
  
  for(i in 1:nrow(parm)){
    parms <- Par(Mu[i], Sigma[i])
    s[i,] <- exp(-eta[i]*(Pesos[i,1]*pgamma(tempo, shape=parms$alpha, scale=parms$beta) + 
                            Pesos[i,2]*plnorm(tempo, parms$mu_l, parms$sig_l)))
  }
  out <- list(surv=colMeans(s), time=tempo)
  return(out)
}

Spop_ptcr.gw <- function(parm, y, X, by) {
  
  if(!is.null(by)) {tempo <- seq(min(y), max(y), by=by)} else{tempo <- y}
  
  s <- matrix(nrow=nrow(parm), ncol=length(tempo)) 
  J <- ncol(parm)
  r <- length(X)
  
  Mu <- parm[,1]
  Sigma <- parm[,2]
  Pesos <- parm[,3:4]
  Betas <- as.matrix(parm[,5:J])
  
  xbeta <- X%*%t(Betas)
  eta <- exp(xbeta)
  
  for(i in 1:nrow(parm)){
    parms <- Par(Mu[i], Sigma[i])
    s[i,] <- exp(-eta[i]*(Pesos[i,1]*pgamma(tempo, shape=parms$alpha, scale=parms$beta) + 
                            Pesos[i,2]*pweibull(tempo, shape=parms$gama, scale=parms$delta)))
  }
  out <- list(surv=colMeans(s), time=tempo)
  return(out)
}

Spop_ptcr.lw <- function(parm, y, X, by) {
  
  if(!is.null(by)) {tempo <- seq(min(y), max(y), by=by)} else{tempo <- y}
  
  s <- matrix(nrow=nrow(parm), ncol=length(tempo)) 
  J <- ncol(parm)
  r <- length(X)
  
  Mu <- parm[,1]
  Sigma <- parm[,2]
  Pesos <- parm[,3:4]
  Betas <- as.matrix(parm[,5:J])
  
  xbeta <- X%*%t(Betas)
  eta <- exp(xbeta)
  
  for(i in 1:nrow(parm)){
    parms <- Par(Mu[i], Sigma[i])
    s[i,] <- exp(-eta[i]*(Pesos[i,1]*plnorm(tempo, parms$mu_l, parms$sig_l) + 
                            Pesos[i,2]*plnorm(tempo, parms$mu_l, parms$sig_l)))
  }
  out <- list(surv=colMeans(s), time=tempo)
  return(out)
}

Spop_ptcr.g <- function(parm, y, X, by) {
  
  if(!is.null(by)) {tempo <- seq(min(y), max(y), by=by)} else{tempo <- y}
  
  s <- matrix(nrow=nrow(parm), ncol=length(tempo)) 
  J <- ncol(parm)
  r <- length(X)
  
  Mu <- parm[,1]
  Sigma <- parm[,2]
  Betas <- as.matrix(parm[,3:J])
  
  xbeta <- X%*%t(Betas)
  eta <- exp(xbeta)
  
  for(i in 1:nrow(parm)){
    parms <- Par(Mu[i], Sigma[i])
    s[i,] <- exp(-eta[i]*pgamma(tempo, shape=parms$alpha, scale=parms$beta))
  }
  out <- list(surv=colMeans(s), time=tempo)
  return(out)
}

Spop_ptcr.l <- function(parm, y, X, by) {
  
  if(!is.null(by)) {tempo <- seq(min(y), max(y), by=by)} else{tempo <- y}
  
  s <- matrix(nrow=nrow(parm), ncol=length(tempo)) 
  J <- ncol(parm)
  r <- length(X)
  
  Mu <- parm[,1]
  Sigma <- parm[,2]
  Betas <- as.matrix(parm[,3:J])
  
  xbeta <- X%*%t(Betas)
  eta <- exp(xbeta)
  
  for(i in 1:nrow(parm)){
    parms <- Par(Mu[i], Sigma[i])
    s[i,] <- exp(-eta[i]*plnorm(tempo, parms$mu_l, parms$sig_l))
  }
  out <- list(surv=colMeans(s), time=tempo)
  return(out)
}

Spop_ptcr.w <- function(parm, y, X, by) {
  
  if(!is.null(by)) {tempo <- seq(min(y), max(y), by=by)} else{tempo <- y}
  
  s <- matrix(nrow=nrow(parm), ncol=length(tempo)) 
  J <- ncol(parm)
  r <- length(X)
  
  Mu <- parm[,1]
  Sigma <- parm[,2]
  Betas <- as.matrix(parm[,3:J])
  
  xbeta <- X%*%t(Betas)
  eta <- exp(xbeta)
  
  for(i in 1:nrow(parm)){
    parms <- Par(Mu[i], Sigma[i])
    s[i,] <- exp(-eta[i]*pweibull(tempo, shape=parms$gama, scale=parms$delta))
  }
  out <- list(surv=colMeans(s), time=tempo)
  return(out)
}
