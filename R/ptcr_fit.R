#------------------------------------------------------------------------------------------------------------------------#
#---------------Functions to generate a MCMC sample from the GLW finite mixture model------------------------------------#
#------------------------------------------------------------------------------------------------------------------------#
#' Fits the Promotion Time Cure Rate model with covariates
#' 
#' \code{glw_ptcr} generates a MCMC sample from the Standard Cure Rate model, with covariates on the cure rate and the GLW FM distribution for the susceptible individuals
#' 
#' @param parm0: vector of initial values, in the order betas, variance, weights, cure rate.
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
glw_ptcr <- function (parm0=NULL, Data, Specs, Specs2=list(Adaptive=1000, Periodicity=100), prior=NULL, mixture="glw") {
  
  y <- Data$y
  N <- length(y)
  X <- as.matrix(Data$X)
  d <- Data$d
  K <- nchar(mixture)
  J <- ncol(X)
  
  mon.names <- "LP"
  if(K>1) {
    parm.names <- as.parm.names(list(mu=0, sigma=0, p=rep(0, K), beta=rep(0, J)))
    
    pos.mu   <- grep("mu", parm.names)
    pos.sigma <- grep("sigma", parm.names)
    pos.p   <- grep("p", parm.names)
    pos.beta 	<- grep("beta", parm.names)
    
    MyData <- list(prior, J=J, K=K, X=X, D=d, y=y, mon.names=mon.names, parm.names=parm.names, pos.beta=pos.beta, pos.p=pos.p, pos.sigma=pos.sigma, pos.mu=pos.mu)
    
    if (is.null(parm0)) parm0 <- c(1, 1, rep(1/K, K), rep(0,J))
  }
  if(K==1) {
    parm.names <- as.parm.names(list(mu=0, sigma=0, beta=rep(0, J)))
    
    pos.mu   <- grep("mu", parm.names)
    pos.sigma <- grep("sigma", parm.names)
    pos.beta 	<- grep("beta", parm.names)
    
    MyData <- list(prior, J=J, K=K, X=X, D=d, y=y, mon.names=mon.names, parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma, pos.mu=pos.mu)
    
    if (is.null(parm0)) parm0 <- c(1, 1, rep(0,J))
  }
  
  Model <- eval(parse(text=paste('PTCR.', mixture, sep='')))
  
  Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values=parm0, Covar=Specs$Covar, Iterations=Specs$Iterations, Status=Specs$Status, 
                       Thinning=Specs$Thinning, Algorithm=Specs$Algorithm, Specs=Specs2, LogFile=Specs$LogFile)
  return(Fit)
}


PTCR.glw <- function(parm, Data) {
  
  ### Parameters
  mu <- interval(parm[Data$pos.mu], 1e-100, Inf)
  parm[Data$pos.mu] <- mu
  
  sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.sigma] <- sigma  
  
  p <- interval(parm[Data$pos.p], 0, 1)
  p <- p/sum(p)
  parm[Data$pos.p] <- p
  
  beta <- parm[Data$pos.beta]
  
  ### Log-Prior
  mu.prior	  <- ifelse(is.null(Data[[1]]$mu),    dgamma(mu, 0.01, 0.01, log=TRUE),    sum(eval(parse(1, text=Data[[1]]$mu))))
  sigma.prior <- ifelse(is.null(Data[[1]]$sigma), dgamma(sigma, 0.01, 0.01, log=TRUE), sum(eval(parse(1, text=Data[[1]]$sigma))) )
  p.prior     <- ifelse(is.null(Data[[1]]$alpha), ddirichlet(p, c(1,1,1), log=TRUE),   ddirichlet(p, Data[[1]]$alpha, log=TRUE))
  beta.prior  <- ifelse(is.null(Data[[1]]$beta),  sum(dnormv(beta, 0, 1000, log=TRUE)),sum(eval(parse(1, text=Data[[1]]$beta))))
  
  ### Log-Likelihood
  X <- Data$X
  y <- Data$y
  D <- Data$D
  J <- Data$J
  
  xbeta <- X%*%as.matrix(beta)
  eta <- exp(xbeta)
  
  a1 <- D*(log(dglw(y, p, mu, sigma)) + xbeta)
  a2 <- -eta*pglw(y, p, mu, sigma)
  
  LL <- sum(a1 + a2)
  
  ### Log-Posterior
  LP <- LL + sum(beta.prior) + sigma.prior + p.prior + mu.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=1, parm=parm)
  return(Modelout)
}

PTCR.gl <- function(parm, Data) {
  
  ### Parameters
  mu <- interval(parm[Data$pos.mu], 1e-100, Inf)
  parm[Data$pos.mu] <- mu
  
  sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.sigma] <- sigma  
  
  p <- interval(parm[Data$pos.p], 0, 1)
  p <- p/sum(p)
  parm[Data$pos.p] <- p
  
  beta <- parm[Data$pos.beta]
  
  ### Log-Prior
  mu.prior	  <- ifelse(is.null(Data[[1]]$mu),    dgamma(mu, 0.01, 0.01, log=TRUE),    sum(eval(parse(1, text=Data[[1]]$mu))))
  sigma.prior <- ifelse(is.null(Data[[1]]$sigma), dgamma(sigma, 0.01, 0.01, log=TRUE), sum(eval(parse(1, text=Data[[1]]$sigma))) )
  p.prior     <- ifelse(is.null(Data[[1]]$alpha), ddirichlet(p, c(1,1), log=TRUE),     ddirichlet(p, Data[[1]]$alpha, log=TRUE))
  beta.prior  <- ifelse(is.null(Data[[1]]$beta),  sum(dnormv(beta, 0, 1000, log=TRUE)),sum(eval(parse(1, text=Data[[1]]$beta))))
  
  ### Log-Likelihood
  X <- Data$X
  y <- Data$y
  D <- Data$D
  J <- Data$J
  
  xbeta <- X%*%as.matrix(beta)
  eta <- exp(xbeta)
  
  parms <- Par(mu, sigma)
  a1 <- D*(log(p[1]*dgamma(y, shape=parms$alpha, scale=parms$beta) + p[2]*dlnorm(y, parms$mu_l, parms$sig_l)) + xbeta)
  a2 <- -eta*(p[1]*pgamma(y, shape=parms$alpha, scale=parms$beta) + p[2]*plnorm(y, parms$mu_l, parms$sig_l))
  
  LL <- sum(a1 + a2)
  
  ### Log-Posterior
  LP <- LL + sum(beta.prior) + sigma.prior + p.prior + mu.prior
  
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=1, parm=parm)
  return(Modelout)
}

PTCR.gw <- function(parm, Data) {
  
  ### Parameters
  mu <- interval(parm[Data$pos.mu], 1e-100, Inf)
  parm[Data$pos.mu] <- mu
  
  sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.sigma] <- sigma  
  
  p <- interval(parm[Data$pos.p], 0, 1)
  p <- p/sum(p)
  parm[Data$pos.p] <- p
  
  beta <- parm[Data$pos.beta]
  
  ### Log-Prior
  mu.prior	  <- ifelse(is.null(Data[[1]]$mu),    dgamma(mu, 0.01, 0.01, log=TRUE),    sum(eval(parse(1, text=Data[[1]]$mu))))
  sigma.prior <- ifelse(is.null(Data[[1]]$sigma), dgamma(sigma, 0.01, 0.01, log=TRUE), sum(eval(parse(1, text=Data[[1]]$sigma))) )
  p.prior     <- ifelse(is.null(Data[[1]]$alpha), ddirichlet(p, c(1,1), log=TRUE),     ddirichlet(p, Data[[1]]$alpha, log=TRUE))
  beta.prior  <- ifelse(is.null(Data[[1]]$beta),  sum(dnormv(beta, 0, 1000, log=TRUE)),sum(eval(parse(1, text=Data[[1]]$beta))))
  
  ### Log-Likelihood
  X <- Data$X
  y <- Data$y
  D <- Data$D
  J <- Data$J
  
  xbeta <- X%*%as.matrix(beta)
  eta <- exp(xbeta)
  
  parms <- Par(mu, sigma)
  a1 <- D*(log(p[1]*dgamma(y, shape=parms$alpha, scale=parms$beta) + p[2]*dweibull(y, shape=parms$gama, scale=parms$delta)) + xbeta)
  a2 <- -eta*(p[1]*pgamma(y, shape=parms$alpha, scale=parms$beta) + p[2]*pweibull(y, shape=parms$gama, scale=parms$delta))
  
  LL <- sum(a1 + a2)
  
  ### Log-Posterior
  LP <- LL + sum(beta.prior) + sigma.prior + p.prior + mu.prior
  
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=1, parm=parm)
  return(Modelout)
}

PTCR.lw <- function(parm, Data) {
  
  ### Parameters
  mu <- interval(parm[Data$pos.mu], 1e-100, Inf)
  parm[Data$pos.mu] <- mu
  
  sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.sigma] <- sigma  
  
  p <- interval(parm[Data$pos.p], 0, 1)
  p <- p/sum(p)
  parm[Data$pos.p] <- p
  
  beta <- parm[Data$pos.beta]
  
  ### Log-Prior
  mu.prior	  <- ifelse(is.null(Data[[1]]$mu),    dgamma(mu, 0.01, 0.01, log=TRUE),    sum(eval(parse(1, text=Data[[1]]$mu))))
  sigma.prior <- ifelse(is.null(Data[[1]]$sigma), dgamma(sigma, 0.01, 0.01, log=TRUE), sum(eval(parse(1, text=Data[[1]]$sigma))) )
  p.prior     <- ifelse(is.null(Data[[1]]$alpha), ddirichlet(p, c(1,1), log=TRUE),     ddirichlet(p, Data[[1]]$alpha, log=TRUE))
  beta.prior  <- ifelse(is.null(Data[[1]]$beta),  sum(dnormv(beta, 0, 1000, log=TRUE)),sum(eval(parse(1, text=Data[[1]]$beta))))
  
  ### Log-Likelihood
  X <- Data$X
  y <- Data$y
  D <- Data$D
  J <- Data$J
  
  xbeta <- X%*%as.matrix(beta)
  eta <- exp(xbeta)
  
  parms <- Par(mu, sigma)
  a1 <- D*(log(p[1]*dlnorm(y, parms$mu_l, parms$sig_l) + p[2]*dweibull(y, shape=parms$gama, scale=parms$delta)) + xbeta)
  a2 <- -eta*(p[1]*plnorm(y, parms$mu_l, parms$sig_l) + p[2]*pweibull(y, shape=parms$gama, scale=parms$delta))
  
  LL <- sum(a1 + a2)
  
  ### Log-Posterior
  LP <- LL + sum(beta.prior) + sigma.prior + p.prior + mu.prior
  
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=1, parm=parm)
  return(Modelout)
}

PTCR.g <- function(parm, Data) {
  
  ### Parameters
  mu <- interval(parm[Data$pos.mu], 1e-100, Inf)
  parm[Data$pos.mu] <- mu
  
  sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.sigma] <- sigma  
  
  beta <- parm[Data$pos.beta]
  
  ### Log-Prior
  mu.prior    <- ifelse(is.null(Data[[1]]$mu),    dgamma(mu, 0.01, 0.01, log=TRUE),    sum(eval(parse(1, text=Data[[1]]$mu))))
  sigma.prior <- ifelse(is.null(Data[[1]]$sigma), dgamma(sigma, 0.01, 0.01, log=TRUE), sum(eval(parse(1, text=Data[[1]]$sigma))) )
  beta.prior  <- ifelse(is.null(Data[[1]]$beta),  sum(dnormv(beta, 0, 1000, log=TRUE)),sum(eval(parse(1, text=Data[[1]]$beta))))
  
  ### Log-Likelihood
  X <- Data$X
  y <- Data$y
  D <- Data$D
  J <- Data$J
  
  xbeta <- X%*%as.matrix(beta)
  eta <- exp(xbeta)
  
  al <- mu^2/sigma
  be  <- sigma/mu
  
  a1 <- D*(dgamma(y, shape=al, scale=be, log=T) + xbeta)
  a2 <- -eta*pgamma(y, shape=al, scale=be)
  
  LL <- sum(a1 + a2)
  
  ### Log-Posterior
  LP <- LL + sum(beta.prior) + sigma.prior + mu.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=1, parm=parm)
  return(Modelout)
}

PTCR.l <- function(parm, Data) {
  
  ### Parameters
  mu <- interval(parm[Data$pos.mu], 1e-100, Inf)
  parm[Data$pos.mu] <- mu
  
  sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.sigma] <- sigma  
  
  beta <- parm[Data$pos.beta]
  
  ### Log-Prior
  mu.prior	  <- ifelse(is.null(Data[[1]]$mu),    dgamma(mu, 0.01, 0.01, log=TRUE),    sum(eval(parse(1, text=Data[[1]]$mu))))
  sigma.prior <- ifelse(is.null(Data[[1]]$sigma), dgamma(sigma, 0.01, 0.01, log=TRUE), sum(eval(parse(1, text=Data[[1]]$sigma))) )
  beta.prior  <- ifelse(is.null(Data[[1]]$beta),  sum(dnormv(beta, 0, 1000, log=TRUE)),sum(eval(parse(1, text=Data[[1]]$beta))))
  
  ### Log-Likelihood
  X <- Data$X
  y <- Data$y
  D <- Data$D
  J <- Data$J
  
  xbeta <- X%*%as.matrix(beta)
  eta <- exp(xbeta)
  
  mu_l  <- log(mu^2/sqrt(mu^2+sigma))
  sig_l <- sqrt(log((mu^2+sigma)/mu^2))
  
  a1 <- D*(dlnorm(y, mu_l, sig_l, log=T) + xbeta)
  a2 <- -eta*plnorm(y, mu_l, sig_l)
  
  LL <- sum(a1 + a2)
  
  ### Log-Posterior
  LP <- LL + sum(beta.prior) + sigma.prior + mu.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=1, parm=parm)
  return(Modelout)
}

PTCR.w <- function(parm, Data) {
  
  ### Parameters
  mu <- interval(parm[Data$pos.mu], 1e-100, Inf)
  parm[Data$pos.mu] <- mu
  
  sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.sigma] <- sigma  
  
  beta <- parm[Data$pos.beta]
  
  ### Log-Prior
  mu.prior	  <- ifelse(is.null(Data[[1]]$mu),    dgamma(mu, 0.01, 0.01, log=TRUE),    sum(eval(parse(1, text=Data[[1]]$mu))))
  sigma.prior <- ifelse(is.null(Data[[1]]$sigma), dgamma(sigma, 0.01, 0.01, log=TRUE), sum(eval(parse(1, text=Data[[1]]$sigma))) )
  beta.prior  <- ifelse(is.null(Data[[1]]$beta),  sum(dnormv(beta, 0, 1000, log=TRUE)),sum(eval(parse(1, text=Data[[1]]$beta))))
  
  ### Log-Likelihood
  X <- Data$X
  y <- Data$y
  D <- Data$D
  J <- Data$J
  
  xbeta <- X%*%as.matrix(beta)
  eta <- exp(xbeta)
  
  parms <- Par(mu, sigma)
  
  a1 <- D*(dweibull(y, shape=parms$gama, scale=parms$delta, log=T) + xbeta)
  a2 <- -eta*pweibull(y, shape=parms$gama, scale=parms$delta)
  
  LL <- sum(a1 + a2)
  
  ### Log-Posterior
  LP <- LL + sum(beta.prior) + sigma.prior + mu.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=1, parm=parm)
  return(Modelout)
}
