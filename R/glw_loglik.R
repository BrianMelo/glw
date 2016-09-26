#---------------------------------------------------------------------#
#------Individual log likelihood--------------------------------------#
#---------------------------------------------------------------------#
#' Individual Log-Likelihood
#' 
#' \code{loglik} computes the log likelihood of every observations for a given vector of parameters
#' 
#' @param parm: vector of parameters
#' @param Data: the same data as in \code{\link{glwfm}}.
#' @param mixture: indicates the finite mixture to be fitted to the data. Options are: 'glw', 'gl', 'gw', 'lw'
#' 
#' @return Returns a vector with the individual contribution of each observations to the log likelihood
#' 
#' @export
loglik <- function(parm, Data, mixture) {
  
  if(is.null(Data$X)) Data$X <- matrix(0, nrow=length(Data$y), ncol=1)
  
  y <- Data$y
  D <- Data$d
  X <- Data$X
  yf <- Data$yf ; if(is.null(yf)) yf <- rep(0, nrow(X))
  rint <- Data$rint
  K <- ifelse(mixture=='glw', 3, 2)
  J <- ifelse(rint==T, K-1, 0) + ncol(X)
  
  parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0, p=rep(0,3)))
  
  pos.beta  <- grep("beta", parm.names)
  pos.sigma <- grep("sigma", parm.names)      
  pos.p     <- grep("p", parm.names)
  
  ### Parameters
  beta <- parm[pos.beta]
  p <- parm[pos.p]
  sigma <- parm[pos.sigma]
  
  if(K==2){
    if(length(grep('gl', mixture))==1) { p[3] <- 0 }
    if(length(grep('gw', mixture))==1) { p[3] <- p[2] ; p[2] <- 0 }
    if(length(grep('lw', mixture))==1) { p[2:3] <- p[1:2] ; p[1] <- 0 }
  }
  
  if(Data$rint==T & K==2){
    if(length(grep('gl', mixture))==1) { beta[4:(J+1)] <- beta[3:J] ; beta[3] <- 0 }
    if(length(grep('gw', mixture))==1) { beta[3:(J+1)] <- beta[2:J] ; beta[2] <- 0 }
    if(length(grep('lw', mixture))==1) { beta[2:(J+1)] <- beta[1:J] ; beta[1] <- 0 }
    J <- J+1
  }
  
  ### Log-Likelihood
  xbeta <- 0 
  if(rint==T){
    if(ncol(X)>1) xbeta <- X[,2:ncol(X)]%*%as.matrix(beta[4:J])
    mu.gamma <- exp(beta[1] + xbeta)
    mu.lnorm <- exp(beta[2] + xbeta)
    mu.weibu <- exp(beta[3] + xbeta)
  }
  if(rint==F) {
    if(ncol(X)>1) xbeta <- X[,2:ncol(X)]%*%as.matrix(beta[2:J])
    mu.gamma <- exp(beta[1] + xbeta)
    mu.lnorm <- mu.gamma
    mu.weibu <- mu.gamma
  }
  parms <- Parm3(mu.gamma, mu.lnorm, mu.weibu, sigma)
  
  a1 <- ifelse(D==1, p[1]*dgamma(y, shape=parms$alpha, scale=parms$beta) + 
                 p[2]*dlnorm(y, parms$mu_l, parms$sig_l) + 
                 p[3]*dweibull(y, shape=parms$gama, scale=parms$delta),1)
  
  a2 <- ifelse(D==0, p[1]*pgamma(y, shape=parms$alpha, scale=parms$beta, lower.tail=F) + 
                 p[2]*plnorm(y, parms$mu_l, parms$sig_l, lower.tail=F) + 
                 p[3]*pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=F),1)
  
  a3 <- ifelse(D==2, p[1]*pgamma(y, shape=parms$alpha, scale=parms$beta, lower.tail=T) + 
                 p[2]*plnorm(y, parms$mu_l, parms$sig_l, lower.tail=T) + 
                 p[3]*pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T),1)
  
  a4 <- ifelse(D==3, p[1]*(pgamma(y, shape=parms$alpha, scale=parms$beta, lower.tail=T) - pgamma(yf, shape=parms$alpha, scale=parms$beta, lower.tail=T)) + 
                 p[2]*(plnorm(y, parms$mu_l, parms$sig_l, lower.tail=T) - plnorm(y, parms$mu_l, parms$sig_l, lower.tail=T)) + 
                 p[3]*(pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T) - pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T)),1)
  
  LL <- log(a1) + log(a2) + log(a3) + log(a4)
  
  return(matrix(LL, nrow=1))
}

#---------------------------------------------------------------------#
#------LogNormal------------------------------------------------------#
#---------------------------------------------------------------------#
lnorm.loglik <- function(parm, Data) {
  
  y <- Data$y
  D <- Data$d
  X <- Data$X
  yf <- Data$yf ; if(is.null(yf)) yf <- rep(0, nrow(X))
  J <- ncol(X)
  
  parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
  
  pos.beta  <- grep("beta", parm.names)
  pos.sigma <- grep("sigma", parm.names)      
  
  ### Parameters
  beta <- parm[pos.beta]
  sigma <- parm[pos.sigma]
  
  ### Log-Likelihood
  mu <- exp(X%*%as.matrix(beta))
  mu_l  <- log(mu^2/sqrt(mu^2+sigma))
  sig_l <- sqrt(log((mu^2+sigma)/mu^2))
  
  a1 <- ifelse(D==1, dlnorm(y, mu_l, sig_l), 1)
  a2 <- ifelse(D==0, plnorm(y, mu_l, sig_l, lower.tail=F), 1)
  a3 <- ifelse(D==2, plnorm(y, mu_l, sig_l, lower.tail=T), 1)
  a4 <- ifelse(D==3, plnorm(y, mu_l, sig_l, lower.tail=T) - plnorm(y, mu_l, sig_l, lower.tail=T), 1)
  
  LL <- log(a1) + log(a2) + log(a3) + log(a4)
  
  return(matrix(LL, nrow=1))
}

#---------------------------------------------------------------------#
#------Weibull--------------------------------------------------------#
#---------------------------------------------------------------------#
weibull.loglik <- function(parm, Data) {
  
  y <- Data$y
  D <- Data$d
  X <- Data$X
  yf <- Data$yf ; if(is.null(yf)) yf <- rep(0, nrow(X))
  J <- ncol(X)
  
  parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
  
  pos.beta  <- grep("beta", parm.names)
  pos.sigma <- grep("sigma", parm.names)      
  
  ### Parameters
  beta <- parm[pos.beta]
  sigma <- parm[pos.sigma]
  
  ### Log-Likelihood
  mu <- exp(X%*%as.matrix(beta))
  parms <- Parm3(mu, mu, mu, sigma)
  
  a1 <- ifelse(D==1, dweibull(y, shape=parms$gama, scale=parms$delta), 1)
  a2 <- ifelse(D==0, pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=F), 1)
  a3 <- ifelse(D==2, pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T), 1)
  a4 <- ifelse(D==3, pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T) - pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T), 1)
  
  LL <- log(a1) + log(a2) + log(a3) + log(a4)
  
  return(matrix(LL, nrow=1))
}

#---------------------------------------------------------------------#
#------Gamma----------------------------------------------------------#
#---------------------------------------------------------------------#
gamma.loglik <- function(parm, Data) {
  
  y <- Data$y
  D <- Data$d
  X <- Data$X
  yf <- Data$yf ; if(is.null(yf)) yf <- rep(0, nrow(X))
  J <- ncol(X)
  
  parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
  
  pos.beta  <- grep("beta", parm.names)
  pos.sigma <- grep("sigma", parm.names)      
  
  ### Parameters
  beta <- parm[pos.beta]
  sigma <- parm[pos.sigma]
  
  ### Log-Likelihood
  mu <- exp(X%*%as.matrix(beta))
  al <- mu^2/sigma
  be <- sigma/mu
  
  a1 <- ifelse(D==1, dgamma(y, shape=al, scale=be), 1)
  a2 <- ifelse(D==0, pgamma(y, shape=al, scale=be, lower.tail=F), 1)
  a3 <- ifelse(D==2, pgamma(y, shape=al, scale=be, lower.tail=T), 1)
  a4 <- ifelse(D==3, pgamma(y, shape=al, scale=be, lower.tail=T) - pgamma(yf, shape=al, scale=be, lower.tail=T), 1)
  
  LL <- log(a1) + log(a2) + log(a3) + log(a4)
  
  return(matrix(LL, nrow=1))
}


