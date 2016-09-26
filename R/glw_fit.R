#------------------------------------------------------------------------------------------------------------------------#
#---------------Functions to generate a MCMC sample from the GLW finite mixture model------------------------------------#
#------------------------------------------------------------------------------------------------------------------------#
#' Fits the GLW Finite Mixture
#' 
#' \code{glwfm} generates a MCMC sample from the GLW finite mixture regression model
#' 
#' @param parm0: vector of initial values, in the order betas, variance, weights.
#' @param Data: list of the data to be analyzed. Must contain: \cr
#'          y: survival times. \cr
#'          d: censorship indicator. The code is 0=right censored, 1=event at time, 2=left censored, 3=interval censored. \cr
#'          X: matrix of covariates (model matrix) or vector of 1's when tere are no covariates.\cr
#'          rint: True or False indicating when to consider different intercept for each component of the finite mixture.
#' @param Specs: list of simulation parameters. Must contain: (Same as Laplaces Demon). \cr
#'        Iterations: number of iterations to be simulated. \cr
#'        Status: number of simulations before showing the status. \cr
#'        Thinning: thinning of the generated sample. \cr
#'        Algorithm: algorith used to generate the posterior samples. \cr
#'        LogFile: path to print the simulation status.
#' @param Specs2: list of simulation parameters specific to the Algorithm chosen (Same as Laplaces Demon).
#' @param prior: list of log of prior distribution for beta an sigma and vector of Dirichlet parameters for the weight. 
#'        Default: list(beta='dnormv(beta, 0, 100, log=T)', sigma='dexp(sigma, 1, log=T)', alpha==c(1,1,1)).
#' @param mixture: indicates the finite mixture to be fitted to the data. Options are: 'glw', 'gl', 'gw', 'lw'.
#' 
#' @return A LaplacesDemon object which contains also a MCMC sample of the posterior distribution of the GLW parameters.
#' 
#' @export
glwfm <- function(parm0=NULL, Data, Specs, Specs2=list(Adaptive=1000, Periodicity=100), prior=NULL, mixture='glw') {
  
  y <- Data$y
  N <- length(y)
  X <- as.matrix(Data$X)
  d <- Data$d
  yf <- Data$yf ; if(is.null(Data$yf)) yf <- rep(0, N)
  rint <- Data$rint
  K <- length(grep('g', mixture)) + length(grep('l', mixture)) + length(grep('w', mixture))
  J <- ifelse(rint==T, K-1, 0) + ncol(X)
  
  mon.names  <- "LP"
  parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0, p=rep(0,K)))
  
  pos.beta  <- grep("beta", parm.names)
  pos.p     <- grep("p", parm.names)
  pos.sigma <- grep("sigma", parm.names)
  
  MyData <- list(prior, J=J, K=K, X=X, D=d, y=y, yf=yf, rint=rint, mon.names=mon.names, parm.names=parm.names, pos.beta=pos.beta, pos.p=pos.p, pos.sigma=pos.sigma)
  
  if(is.null(parm0)) parm0 <- c(rep(0,J), 1, rep(1/K,K))
  
  Model <- eval(parse(text=paste('Model.', mixture, sep='')))
  
  Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values=parm0, Covar=Specs$Covar, Iterations=Specs$Iterations, Status=Specs$Status, Thinning=Specs$Thinning,  Algorithm=Specs$Algorithm, Specs=Specs2, LogFile=Specs$LogFile)
  return(Fit)
}


Model.gl <- function(parm, Data) {
  
  ### Parameters
  beta <- parm[Data$pos.beta]
  
  p <- interval(parm[Data$pos.p], 0, 1)
  p <- p/sum(p)
  parm[Data$pos.p] <- p
  
  sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.sigma] <- sigma
  
  ### Log-Prior
  beta.prior  <- ifelse(is.null(Data[[1]]$beta), sum(dnormv(beta, 0, 1000, log=TRUE)), sum(eval(parse(1, text=Data[[1]]$beta))) )
  sigma.prior <- ifelse(is.null(Data[[1]]$sigma), sum(dexp(sigma, 1, log=TRUE)), sum(eval(parse(1, text=Data[[1]]$sigma))) )
  p.prior     <- ifelse(is.null(Data[[1]]$alpha), ddirichlet(p, c(1,1), log=TRUE), ddirichlet(p, Data[[1]]$alpha, log=TRUE) )
  
  ### Log-Likelihood
  X <- Data$X
  y <- Data$y
  D <- Data$D
  yf <- Data$yf
  J  <- Data$J
  
  xbeta <- 0
  if(Data$rint==T){
    if(ncol(X)>1) xbeta <- X[,2:ncol(X)]%*%as.matrix(beta[3:J])
    mu.gamma <- exp(beta[1] + xbeta)
    mu.lnorm <- exp(beta[2] + xbeta)
    mu.weibu <- mu.gamma
  }
  if(Data$rint==F) {
    if(ncol(X)>1) xbeta <- X[,2:ncol(X)]%*%as.matrix(beta[2:J])
    mu.gamma <- exp(beta[1] + xbeta)
    mu.lnorm <- mu.gamma
    mu.weibu <- mu.gamma
  }
  parms <- Parm3(mu.gamma, mu.lnorm, mu.weibu, sigma)
  
  a1 <- ifelse(D==1, p[1]*dgamma(y, shape=parms$alpha, scale=parms$beta) + p[2]*dlnorm(y, parms$mu_l, parms$sig_l),1)
  a2 <- ifelse(D==0, p[1]*pgamma(y, shape=parms$alpha, scale=parms$beta, lower.tail=F) + p[2]*plnorm(y, parms$mu_l, parms$sig_l, lower.tail=F),1)
  a3 <- ifelse(D==2, p[1]*pgamma(y, shape=parms$alpha, scale=parms$beta, lower.tail=T) + p[2]*plnorm(y, parms$mu_l, parms$sig_l, lower.tail=T),1)
  a4 <- ifelse(D==3, p[1]*(pgamma(y, shape=parms$alpha, scale=parms$beta, lower.tail=T) - pgamma(yf, shape=parms$alpha, scale=parms$beta, lower.tail=T)) + p[2]*(plnorm(y, parms$mu_l, parms$sig_l, lower.tail=T) - plnorm(y, parms$mu_l, parms$sig_l, lower.tail=T)),1)
  LL <- sum(log(a1) + log(a2) + log(a3) + log(a4) )
  
  ### Log-Posterior
  LP <- LL + sum(beta.prior) + sigma.prior + p.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rglw(length(mu.gamma), c(p[1], p[2], 0), mu.gamma, sigma), parm=parm)
  return(Modelout)
}

Model.gw <- function(parm, Data) {
  
  ### Parameters
  beta <- parm[Data$pos.beta]
  p <- interval(parm[Data$pos.p], 0, 1)
  p <- p/sum(p)
  parm[Data$pos.p] <- p
  
  sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.sigma] <- sigma
  
  ### Log-Prior
  beta.prior  <- ifelse(is.null(Data[[1]]$beta), sum(dnormv(beta, 0, 1000, log=TRUE)), sum(eval(parse(1, text=Data[[1]]$beta))) )
  sigma.prior <- ifelse(is.null(Data[[1]]$sigma), sum(dexp(sigma, 1, log=TRUE)), sum(eval(parse(1, text=Data[[1]]$sigma))) )
  p.prior     <- ifelse(is.null(Data[[1]]$alpha), ddirichlet(p, c(1,1), log=TRUE), ddirichlet(p, Data[[1]]$alpha, log=TRUE) )
  
  ### Log-Likelihood
  X <- Data$X
  y <- Data$y
  D <- Data$D
  yf <- Data$yf
  J  <- Data$J
  
  xbeta <- 0
  if(Data$rint==T){
    if(ncol(X)>1) xbeta <- X[,2:ncol(X)]%*%as.matrix(beta[3:J])
    mu.gamma <- exp(beta[1] + xbeta)
    mu.lnorm <- mu.gamma
    mu.weibu <- exp(beta[2] + xbeta)
  }
  if(Data$rint==F) {
    if(ncol(X)>1) xbeta <- X[,2:ncol(X)]%*%as.matrix(beta[2:J])
    mu.gamma <- exp(beta[1] + xbeta)
    mu.lnorm <- mu.gamma
    mu.weibu <- mu.gamma
  }
  parms <- Parm3(mu.gamma, mu.lnorm, mu.weibu, sigma)
  
  a1 <- ifelse(D==1, p[1]*dgamma(y, shape=parms$alpha, scale=parms$beta) + p[2]*dweibull(y, shape=parms$gama, scale=parms$delta),1)
  a2 <- ifelse(D==0, p[1]*pgamma(y, shape=parms$alpha, scale=parms$beta, lower.tail=F) + p[2]*pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=F),1)
  a3 <- ifelse(D==2, p[1]*pgamma(y, shape=parms$alpha, scale=parms$beta, lower.tail=T) + p[2]*pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T),1)
  a4 <- ifelse(D==3, p[1]*(pgamma(y, shape=parms$alpha, scale=parms$beta, lower.tail=T) - pgamma(yf, shape=parms$alpha, scale=parms$beta, lower.tail=T)) + p[2]*(pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T) - pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T)),1)
  LL <- sum(log(a1) + log(a2) + log(a3) + log(a4) )
  
  ### Log-Posterior
  LP <- LL + sum(beta.prior) + sigma.prior + p.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rglw(length(mu.gamma), c(p[1], 0, p[2]), mu.gamma, sigma), parm=parm)
  return(Modelout)
}


Model.lw <- function(parm, Data) {
  
  ### Parameters
  beta <- parm[Data$pos.beta]
  
  p <- interval(parm[Data$pos.p], 0, 1)
  p <- p/sum(p)
  parm[Data$pos.p] <- p
  
  sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.sigma] <- sigma
  
  ### Log-Prior
  beta.prior  <- ifelse(is.null(Data[[1]]$beta), sum(dnormv(beta, 0, 1000, log=TRUE)), sum(eval(parse(1, text=Data[[1]]$beta))) )
  sigma.prior <- ifelse(is.null(Data[[1]]$sigma), sum(dexp(sigma, 1, log=TRUE)), sum(eval(parse(1, text=Data[[1]]$sigma))) )
  p.prior     <- ifelse(is.null(Data[[1]]$alpha), ddirichlet(p, c(1,1), log=TRUE), ddirichlet(p, Data[[1]]$alpha, log=TRUE) )
  
  ### Log-Likelihood
    X <- Data$X
  y <- Data$y
  D <- Data$D
  yf <- Data$yf
  J  <- Data$J
  
  xbeta <- 0
  if(Data$rint==T){
    if(ncol(X)>1) xbeta <- X[,2:ncol(X)]%*%as.matrix(beta[3:J])
    mu.gamma <- exp(xbeta)
    mu.lnorm <- exp(beta[1] + xbeta)
    mu.weibu <- exp(beta[2] + xbeta)
  }
  if(Data$rint==F) {
    if(ncol(X)>1) xbeta <- X[,2:ncol(X)]%*%as.matrix(beta[2:J])
    mu.gamma <- exp(beta[1] + xbeta)
    mu.lnorm <- mu.gamma
    mu.weibu <- mu.gamma
  }
  parms <- Parm3(mu.gamma, mu.lnorm, mu.weibu, sigma)
  
  a1 <- ifelse(D==1, p[1]*dlnorm(y, parms$mu_l, parms$sig_l) + p[2]*dweibull(y, shape=parms$gama, scale=parms$delta),1)
  a2 <- ifelse(D==0, p[1]*plnorm(y, parms$mu_l, parms$sig_l, lower.tail=F) + p[2]*pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=F),1)
  a3 <- ifelse(D==2, p[1]*plnorm(y, parms$mu_l, parms$sig_l, lower.tail=T) + p[2]*pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T),1)
  a4 <- ifelse(D==3, p[1]*(plnorm(y, parms$mu_l, parms$sig_l, lower.tail=T) - plnorm(y, parms$mu_l, parms$sig_l, lower.tail=T)) + p[2]*(pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T) - pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T)),1)
  LL <- sum(log(a1) + log(a2) + log(a3) + log(a4) )
  
  ### Log-Posterior
  LP <- LL + sum(beta.prior) + sigma.prior + p.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rglw(length(mu.gamma), c(0, p[1], p[2]), mu.gamma, sigma), parm=parm)
  return(Modelout)
}

Model.glw <- function(parm, Data) {
  
  ### Parameters
  beta <- parm[Data$pos.beta]
  
  p <- interval(parm[Data$pos.p], 0, 1)
  p <- p/sum(p)
  parm[Data$pos.p] <- p
  
  sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.sigma] <- sigma
  
  ### Log-Prior
  beta.prior  <- ifelse(is.null(Data[[1]]$beta), sum(dnormv(beta, 0, 1000, log=TRUE)), sum(eval(parse(1, text=Data[[1]]$beta))) )
  sigma.prior <- ifelse(is.null(Data[[1]]$sigma), sum(dexp(sigma, 1, log=TRUE)), sum(eval(parse(1, text=Data[[1]]$sigma))) )
  p.prior     <- ifelse(is.null(Data[[1]]$alpha), ddirichlet(p, c(1,1,1), log=TRUE), ddirichlet(p, Data[[1]]$alpha, log=TRUE) )
  
  ### Log-Likelihood
  X <- Data$X
  y <- Data$y
  D <- Data$D
  yf <- Data$yf
  J <- Data$J
  
  xbeta <- 0
  if(Data$rint==T){
    if(ncol(X)>1) xbeta <- X[,2:ncol(X)]%*%as.matrix(beta[4:J])
    mu.gamma <- exp(beta[1] + xbeta)
    mu.lnorm <- exp(beta[2] + xbeta)
    mu.weibu <- exp(beta[3] + xbeta)
  }
  if(Data$rint==F) {
    if(ncol(X)>1) xbeta <- X[,2:ncol(X)]%*%as.matrix(beta[2:J])
    mu.gamma <- exp(beta[1] + xbeta)
    mu.lnorm <- mu.gamma
    mu.weibu <- mu.gamma
  }
  parms <- Parm3(mu.gamma, mu.lnorm, mu.weibu, sigma)
  
  a1 <- ifelse(D==1, p[1]*dgamma(y, shape=parms$alpha, scale=parms$beta) + p[2]*dlnorm(y, parms$mu_l, parms$sig_l) + p[3]*dweibull(y, shape=parms$gama, scale=parms$delta),1)
  a2 <- ifelse(D==0, p[1]*pgamma(y, shape=parms$alpha, scale=parms$beta, lower.tail=F) + p[2]*plnorm(y, parms$mu_l, parms$sig_l, lower.tail=F) + p[3]*pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=F),1)
  a3 <- ifelse(D==2, p[1]*pgamma(y, shape=parms$alpha, scale=parms$beta, lower.tail=T) + p[2]*plnorm(y, parms$mu_l, parms$sig_l, lower.tail=T) + p[3]*pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T),1)
  a4 <- ifelse(D==3, p[1]*(pgamma(y, shape=parms$alpha, scale=parms$beta, lower.tail=T) - pgamma(yf, shape=parms$alpha, scale=parms$beta, lower.tail=T)) + p[2]*(plnorm(y, parms$mu_l, parms$sig_l, lower.tail=T) - plnorm(y, parms$mu_l, parms$sig_l, lower.tail=T)) + p[3]*(pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T) - pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T)),1)
  LL <- sum(log(a1) + log(a2) + log(a3) + log(a4) )
  
  ### Log-Posterior
  LP <- LL + sum(beta.prior) + sigma.prior + p.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rglw(length(mu.gamma), p, mu.gamma, sigma), parm=parm)
  return(Modelout)
}
