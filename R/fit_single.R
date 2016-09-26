#----------------------------------------------------#
#--------------GAMMA---------------------------------#
#----------------------------------------------------#
#' Fits the Gamma Model
#' 
#' \code{gamma_glw} generates a MCMC sample considering a gama distribution for the response variable
#' 
#' @param parm0: vector of initial values, in the order betas, variance, weights.
#' @param Data: list of the data to be analyzed. Must contain: \cr
#'          y: survival times. \cr
#'          d: censorship indicator. The code is 0=right censored, 1=event at time, 2=left censored, 3=interval censored. \cr
#'          X: matrix of covariates (model matrix) or vector of 1's when tere are no covariates.\cr
#'          rint: True or False indicating when to consider different intercept for each component of the finite mixture.
#' @param Specs: list of simulation parameters. Must contain: (Same as Laplaces Demon). \cr
#'          Iterations: number of iterations to be simulated. \cr
#'          Status: number of simulations before showing the status. \cr
#'          Thinning: thinning of the generated sample. \cr
#'          Algorithm: algorith used to generate the posterior samples. \cr
#'          LogFile: path to print the simulation status.
#' @param Specs2: list of simulation parameters specific to the Algorithm chosen (Same as Laplaces Demon).
#' @param prior: list of log of prior distribution for beta an sigma and vector of Dirichlet parameters for the weight. 
#'        Default: list(beta='dnormv(beta, 0, 100, log=T)', sigma='dexp(sigma, 1, log=T)', alpha==c(1,1,1)).
#' 
#' @return A LaplacesDemon object which contains also a MCMC sample of the posterior distribution of the GLW parameters.
#' 
#' @export
gamma_glw <- function(parm0, Data, Specs, Specs2=list(Adaptive=1000, Periodicity=100), prior=NULL) {
  
  y <- Data$y
  N <- length(y)
  X <- as.matrix(Data$X)
  d <- Data$d
  yf <- Data$yf ; if(is.null(Data$yf)) yf <- rep(0, N)
  J <- ncol(X)
  
  mon.names  <- "LP"
  parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
  
  pos.beta  <- grep("beta", parm.names)
  pos.sigma <- grep("sigma", parm.names)
  
  MyData <- list(prior, J=J, X=X, D=d, y=y, yf=yf, mon.names=mon.names, parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma)
  
  if(is.null(parm0)) parm0 <- c(rep(0,J), 1)
  
  Model <- function(parm, Data) {
    
    ### Parameters
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    
    ### Log-Prior
    beta.prior  <- ifelse(is.null(Data[[1]]$beta), sum(dnormv(beta, 0, 1000, log=TRUE)), sum(eval(parse(1, text=Data[[1]]$beta))) )
    sigma.prior <- ifelse(is.null(Data[[1]]$sigma), sum(dexp(sigma, 1, log=TRUE)), sum(eval(parse(1, text=Data[[1]]$sigma))) )
    
    ### Log-Likelihood
    X <- Data$X
    y <- Data$y
    D <- Data$D
    yf <- Data$yf
    J  <- Data$J
    
    mu <- exp(X%*%as.matrix(beta))
    al <- mu^2/sigma
    be <- sigma/mu
    
    a1 <- ifelse(D==1, dgamma(y, shape=al, scale=be), 1)
    a2 <- ifelse(D==0, pgamma(y, shape=al, scale=be, lower.tail=F), 1)
    a3 <- ifelse(D==2, pgamma(y, shape=al, scale=be, lower.tail=T), 1)
    a4 <- ifelse(D==3, pgamma(y, shape=al, scale=be, lower.tail=T) - pgamma(yf, shape=al, scale=be, lower.tail=T), 1)
    
    LL <- sum(log(a1) + log(a2) + log(a3) + log(a4) )
    
    ### Log-Posterior
    LP <- LL + sum(beta.prior) + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rgamma(length(mu), shape=mu^2/sigma, scale=sigma/mu), parm=parm)
    return(Modelout)
  }
  Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values=parm0, Covar=Specs$Covar, Iterations=Specs$Iterations, Status=Specs$Status, Thinning=Specs$Thinning,  Algorithm=Specs$Algorithm, Specs=Specs2, LogFile=Specs$LogFile)
  return(Fit)
}

#----------------------------------------------------#
#--------------WEIBULL-------------------------------#
#----------------------------------------------------#
#' Fits the Weibull Model
#' 
#' \code{weibull_glw} generates a MCMC sample considering a weibull distribution for the response variable
#' 
#' @param parm0: vector of initial values, in the order betas, variance, weights.
#' @param Data: list of the data to be analyzed. Must contain: \cr
#'          y: survival times. \cr
#'          d: censorship indicator. The code is 0=right censored, 1=event at time, 2=left censored, 3=interval censored. \cr
#'          X: matrix of covariates (model matrix) or vector of 1's when tere are no covariates.\cr
#'          rint: True or False indicating when to consider different intercept for each component of the finite mixture.
#' @param Specs: list of simulation parameters. Must contain: (Same as Laplaces Demon). \cr
#'          Iterations: number of iterations to be simulated. \cr
#'          Status: number of simulations before showing the status. \cr
#'          Thinning: thinning of the generated sample. \cr
#'          Algorithm: algorith used to generate the posterior samples. \cr
#'          LogFile: path to print the simulation status.
#' @param Specs2: list of simulation parameters specific to the Algorithm chosen (Same as Laplaces Demon).
#' @param prior: list of log of prior distribution for beta an sigma and vector of Dirichlet parameters for the weight. 
#'        Default: list(beta='dnormv(beta, 0, 100, log=T)', sigma='dexp(sigma, 1, log=T)', alpha==c(1,1,1)).
#' 
#' @return A LaplacesDemon object which contains also a MCMC sample of the posterior distribution of the GLW parameters.
#' 
#' @export
weibull_glw <- function(parm0, Data, Specs, Specs2=list(Adaptive=1000, Periodicity=100), prior=NULL) {
  
  y <- Data$y
  N <- length(y)
  X <- as.matrix(Data$X)
  d <- Data$d
  yf <- Data$yf ; if(is.null(Data$yf)) yf <- rep(0, N)
  J <- ncol(X)
  
  mon.names  <- "LP"
  parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
  
  pos.beta  <- grep("beta", parm.names)
  pos.sigma <- grep("sigma", parm.names)
  
  MyData <- list(prior, J=J, X=X, D=d, y=y, yf=yf, mon.names=mon.names, parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma)
  
  if(is.null(parm0)) parm0 <- c(rep(0,J), 1)
  
  Model <- function(parm, Data) {
    
    ### Parameters
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    
    ### Log-Prior
    beta.prior  <- ifelse(is.null(Data[[1]]$beta), sum(dnormv(beta, 0, 1000, log=TRUE)), sum(eval(parse(1, text=Data[[1]]$beta))) )
    sigma.prior <- ifelse(is.null(Data[[1]]$sigma), sum(dexp(sigma, 1, log=TRUE)), sum(eval(parse(1, text=Data[[1]]$sigma))) )
    
    ### Log-Likelihood
    X <- Data$X
    y <- Data$y
    D <- Data$D
    yf <- Data$yf
    J  <- Data$J
    
    mu <- exp(X%*%as.matrix(beta))
    parms <- Parm3(mu, mu, mu, sigma)
    
    a1 <- ifelse(D==1, dweibull(y, shape=parms$gama, scale=parms$delta), 1)
    a2 <- ifelse(D==0, pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=F), 1)
    a3 <- ifelse(D==2, pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T), 1)
    a4 <- ifelse(D==3, pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T) - pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T), 1)
    
    LL <- sum(log(a1) + log(a2) + log(a3) + log(a4) )
    
    ### Log-Posterior
    LP <- LL + sum(beta.prior) + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rep(0,length(mu)), parm=parm)
    return(Modelout)
  }
  Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values=parm0, Covar=Specs$Covar, Iterations=Specs$Iterations, Status=Specs$Status, Thinning=Specs$Thinning,  Algorithm=Specs$Algorithm, Specs=Specs2, LogFile=Specs$LogFile)
  return(Fit)
}

#----------------------------------------------------#
#------------LOG-NORMAL------------------------------#
#----------------------------------------------------#
#' Fits the Lognormal Model
#' 
#' \code{lnorm_glw} generates a MCMC sample considering a lognormal distribution for the response variable
#' 
#' @param parm0: vector of initial values, in the order betas, variance, weights.
#' @param Data: list of the data to be analyzed. Must contain: \cr
#'          y: survival times. \cr
#'          d: censorship indicator. The code is 0=right censored, 1=event at time, 2=left censored, 3=interval censored. \cr
#'          X: matrix of covariates (model matrix) or vector of 1's when tere are no covariates.\cr
#'          rint: True or False indicating when to consider different intercept for each component of the finite mixture.
#' @param Specs: list of simulation parameters. Must contain: (Same as Laplaces Demon). \cr
#'          Iterations: number of iterations to be simulated. \cr
#'          Status: number of simulations before showing the status. \cr
#'          Thinning: thinning of the generated sample. \cr
#'          Algorithm: algorith used to generate the posterior samples. \cr
#'          LogFile: path to print the simulation status.
#' @param Specs2: list of simulation parameters specific to the Algorithm chosen (Same as Laplaces Demon).
#' @param prior: list of log of prior distribution for beta an sigma and vector of Dirichlet parameters for the weight. 
#'        Default: list(beta='dnormv(beta, 0, 100, log=T)', sigma='dexp(sigma, 1, log=T)', alpha==c(1,1,1)).
#' 
#' @return A LaplacesDemon object which contains also a MCMC sample of the posterior distribution of the GLW parameters.
#' 
#' @export
lnorm_glw <- function(parm0, Data, Specs, Specs2=list(Adaptive=1000, Periodicity=100), prior=NULL) {
  
  y <- Data$y
  N <- length(y)
  X <- as.matrix(Data$X)
  d <- Data$d
  yf <- Data$yf ; if(is.null(Data$yf)) yf <- rep(0, N)
  J <- ncol(X)
  
  mon.names  <- "LP"
  parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
  
  pos.beta  <- grep("beta", parm.names)
  pos.sigma <- grep("sigma", parm.names)
  
  MyData <- list(prior, J=J, X=X, D=d, y=y, yf=yf, mon.names=mon.names, parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma)
  
  if(is.null(parm0)) parm0 <- c(rep(0,J), 1)
  
  Model <- function(parm, Data) {
    
    ### Parameters
    beta <- parm[Data$pos.beta]
    sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
    parm[Data$pos.sigma] <- sigma
    
    ### Log-Prior
    beta.prior  <- ifelse(is.null(Data[[1]]$beta), sum(dnormv(beta, 0, 1000, log=TRUE)), sum(eval(parse(1, text=Data[[1]]$beta))) )
    sigma.prior <- ifelse(is.null(Data[[1]]$sigma), sum(dexp(sigma, 1, log=TRUE)), sum(eval(parse(1, text=Data[[1]]$sigma))) )
    
    ### Log-Likelihood
    X <- Data$X
    y <- Data$y
    D <- Data$D
    yf <- Data$yf
    J  <- Data$J
    
    mu <- exp(X%*%as.matrix(beta))
    mu_l  <- log(mu^2/sqrt(mu^2+sigma))
    sig_l <- sqrt(log((mu^2+sigma)/mu^2))
    
    a1 <- ifelse(D==1, dlnorm(y, mu_l, sig_l), 1)
    a2 <- ifelse(D==0, plnorm(y, mu_l, sig_l, lower.tail=F), 1)
    a3 <- ifelse(D==2, plnorm(y, mu_l, sig_l, lower.tail=T), 1)
    a4 <- ifelse(D==3, plnorm(y, mu_l, sig_l, lower.tail=T) - plnorm(y, mu_l, sig_l, lower.tail=T), 1)
    
    LL <- sum(log(a1) + log(a2) + log(a3) + log(a4) )
    
    ### Log-Posterior
    LP <- LL + sum(beta.prior) + sigma.prior
    Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rlnorm(length(mu), log(mu^2/sqrt(mu^2+sigma)), sqrt(log((mu^2+sigma)/mu^2))), parm=parm)
    return(Modelout)
  }
  Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values=parm0, Covar=Specs$Covar, Iterations=Specs$Iterations, Status=Specs$Status, Thinning=Specs$Thinning,  Algorithm=Specs$Algorithm, Specs=Specs2, LogFile=Specs$LogFile)
  return(Fit)
}

