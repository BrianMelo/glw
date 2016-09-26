#------------------------------------------------------------------------------------------------------------------------#
#---------------Estimated survival functions og the GLW finite mixture modelo--------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------#
#'Estimated Survival Function of GLW Standard Cure Rate Model with covariates
#'
#'\code{Spop_scr} estimates the survival function of the Standard Cure Rate Model, with covariates, using the GLW FM distribution for the susceptible individuals 
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
Spop_scr <- function(parm, y, X, mixture, by=NULL) {
  out <- eval(parse(text=paste('Spop_scr.', mixture, '(parm=parm,y=y, X=X, by=by)', sep='')))
  return(out)
}

Spop_scr.glw <- function(parm, y, X, by) {
  
  if(!is.null(by)) {tempo <- seq(min(y), max(y), by=by)} else{tempo <- y}
  
  s <- matrix(nrow=nrow(parm), ncol=length(tempo)) 
  J <- ncol(parm)
  r <- length(X)
  
  Mu <- parm[,1]
  Sigma <- parm[,2]
  Pesos <- parm[,3:5]
  Betas <- as.matrix(parm[,6:J])
  
  xbeta <- X%*%t(Betas)
  Pi <- exp(xbeta)/(1+exp(xbeta))
  
  for(i in 1:nrow(parm)){
    parms <- Par(Mu[i], Sigma[i])
    s[i,] <- Pi[i] + (1-Pi[i])*(Pesos[i,1]*pgamma(tempo, shape=parms$alpha, scale=parms$beta, lower.tail=F) + 
                                  Pesos[i,2]*plnorm(tempo, parms$mu_l, parms$sig_l, lower.tail=F) + 
                                  Pesos[i,3]*pweibull(tempo, shape=parms$gama, scale=parms$delta, lower.tail=F))
  }
  out <- list(surv=colMeans(s), time=tempo)
  return(out)
}

Spop_scr.gl <- function(parm, y, X, by) {
  
  if(!is.null(by)) {tempo <- seq(min(y), max(y), by=by)} else{tempo <- y}
  
  s <- matrix(nrow=nrow(parm), ncol=length(tempo)) 
  J <- ncol(parm)
  r <- length(X)
  
  Mu <- parm[,1]
  Sigma <- parm[,2]
  Pesos <- parm[,3:4]
  Betas <- as.matrix(parm[,5:J])
  
  xbeta <- X%*%t(Betas)
  Pi <- exp(xbeta)/(1+exp(xbeta))
  
  for(i in 1:nrow(parm)){
    parms <- Par(Mu[i], Sigma[i])
    s[i,] <- Pi[i] + (1-Pi[i])*(Pesos[i,1]*pgamma(tempo, shape=parms$alpha, scale=parms$beta, lower.tail=F) + 
                                  Pesos[i,2]*plnorm(tempo, parms$mu_l, parms$sig_l, lower.tail=F))
  }
  out <- list(surv=colMeans(s), time=tempo)
  return(out)
}

Spop_scr.gw <- function(parm, y, X, by) {
  
  if(!is.null(by)) {tempo <- seq(min(y), max(y), by=by)} else{tempo <- y}
  
  s <- matrix(nrow=nrow(parm), ncol=length(tempo)) 
  J <- ncol(parm)
  r <- length(X)
  
  Mu <- parm[,1]
  Sigma <- parm[,2]
  Pesos <- parm[,3:4]
  Betas <- as.matrix(parm[,5:J])
  
  xbeta <- X%*%t(Betas)
  Pi <- exp(xbeta)/(1+exp(xbeta))
  
  for(i in 1:nrow(parm)){
    parms <- Par(Mu[i], Sigma[i])
    s[i,] <- Pi[i] + (1-Pi[i])*(Pesos[i,1]*pgamma(tempo, shape=parms$alpha, scale=parms$beta, lower.tail=F) + 
                                  Pesos[i,2]*pweibull(tempo, shape=parms$gama, scale=parms$delta, lower.tail=F))
  }
  out <- list(surv=colMeans(s), time=tempo)
  return(out)
}

Spop_scr.lw <- function(parm, y, X, by) {
  
  if(!is.null(by)) {tempo <- seq(min(y), max(y), by=by)} else{tempo <- y}
  
  s <- matrix(nrow=nrow(parm), ncol=length(tempo)) 
  J <- ncol(parm)
  r <- length(X)
  
  Mu <- parm[,1]
  Sigma <- parm[,2]
  Pesos <- parm[,3:4]
  Betas <- as.matrix(parm[,5:J])
  
  xbeta <- X%*%t(Betas)
  Pi <- exp(xbeta)/(1+exp(xbeta))
  
  for(i in 1:nrow(parm)){
    parms <- Par(Mu[i], Sigma[i])
    s[i,] <- Pi[i] + (1-Pi[i])*(Pesos[i,1]*plnorm(tempo, parms$mu_l, parms$sig_l, lower.tail=F) + 
                                  Pesos[i,2]*plnorm(tempo, parms$mu_l, parms$sig_l, lower.tail=F))
  }
  out <- list(surv=colMeans(s), time=tempo)
  return(out)
}


Spop_scr.g <- function(parm, y, X, by) {
  
  if(!is.null(by)) {tempo <- seq(min(y), max(y), by=by)} else{tempo <- y}
  
  s <- matrix(nrow=nrow(parm), ncol=length(tempo)) 
  J <- ncol(parm)
  r <- length(X)
  
  Mu <- parm[,1]
  Sigma <- parm[,2]
  Betas <- as.matrix(parm[,3:J])
  
  xbeta <- X%*%t(Betas)
  Pi <- exp(xbeta)/(1+exp(xbeta))
  
  al <- Mu^2/Sigma
  be <- Sigma/Mu
  
  for(i in 1:nrow(parm)){
    
    s[i,] <- Pi[i] + (1-Pi[i])*pgamma(tempo, shape=al[i], scale=be[i], lower.tail=F)
  }
  out <- list(surv=colMeans(s), time=tempo)
  return(out)
}

Spop_scr.l <- function(parm, y, X, by) {
  
  if(!is.null(by)) {tempo <- seq(min(y), max(y), by=by)} else{tempo <- y}
  
  s <- matrix(nrow=nrow(parm), ncol=length(tempo)) 
  J <- ncol(parm)
  r <- length(X)
  
  Mu <- parm[,1]
  Sigma <- parm[,2]
  Betas <- as.matrix(parm[,3:J])
  
  xbeta <- X%*%t(Betas)
  Pi <- exp(xbeta)/(1+exp(xbeta))
  
  mu_l  <- log(Mu^2/sqrt(Mu^2+Sigma))
  sig_l <- sqrt(log((Mu^2+Sigma)/Mu^2))
  
  for(i in 1:nrow(parm)){
    
    s[i,] <- Pi[i] + (1-Pi[i])*plnorm(tempo, mu_l[i], sig_l[i], lower.tail=F)
  }
  out <- list(surv=colMeans(s), time=tempo)
  return(out)
}

Spop_scr.w <- function(parm, y, X, by) {
  
  if(!is.null(by)) {tempo <- seq(min(y), max(y), by=by)} else{tempo <- y}
  
  s <- matrix(nrow=nrow(parm), ncol=length(tempo)) 
  J <- ncol(parm)
  r <- length(X)
  
  Mu <- parm[,1]
  Sigma <- parm[,2]
  Betas <- as.matrix(parm[,3:J])
  
  xbeta <- X%*%t(Betas)
  Pi <- exp(xbeta)/(1+exp(xbeta))
  
  for(i in 1:nrow(parm)){
    parms <- Par(Mu[i], Sigma[i])
    s[i,] <- Pi[i] + (1-Pi[i])*pweibull(tempo, shape=parms$gama, scale=parms$delta, lower.tail=F)
  }
  out <- list(surv=colMeans(s), time=tempo)
  return(out)
}