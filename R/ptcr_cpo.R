#---------------------------------------------------------------------#
#------CPO - Conditional Predictive Ordinate--------------------------#
#---------------------------------------------------------------------#
#'CPO of the GLW FM Promotion Time Cure Rate Model
#'
#' \code{cpo_ptcr} estimates the CPO based on a MCMC sample of the posterior distribuiton.
#' 
#' @param parm: sample of posterior parameters.
#' @param Data: the same data as in \code{\link{glw_scr}}.
#' @param mixture: indicates the finite mixture to be fitted to the data. Options are: 'glw', 'gl', 'gw', 'lw'.
#' 
#' @return Returns the Conditional Predictive Ordinate of the GLW finite mixture model for every observation.
#' 
#' @export
cpo_ptcr <- function (parm, Data, mixture) {
  cpo <- matrix(nrow = nrow(parm), ncol = nrow(Data$X))
  if(nchar(mixture)>1) {
    for (i in 1:nrow(parm)) cpo[i, ] <- 1/exp(loglik_ptcr(parm = parm[i,], Data = Data, mixture = mixture))
  }
  if(nchar(mixture)==1) {
    for (i in 1:nrow(parm)) {
      aux <- eval(parse(text=paste('loglik_ptcr.', mixture, '(parm = parm[i,], Data = Data)', sep='')))
      cpo[i, ] <- 1/exp(aux)    
    }}
  out <- 1/colMeans(cpo)
  return(out)
}

loglik_ptcr <- function (parm, Data, mixture) {
  
  if (is.null(Data$X)) Data$X <- matrix(0, nrow = length(Data$y), ncol = 1)
  
  y <- Data$y
  D <- Data$d
  X <- as.matrix(Data$X)
  K <- nchar(mixture)
  J <- ncol(X)
  
  parm.names <- as.parm.names(list(mu=0, sigma=0, p=rep(0, K), beta=rep(0, J)))
  
  pos.mu <- grep("mu", parm.names)
  pos.sigma <- grep("sigma", parm.names)
  pos.p <- grep("p", parm.names)
  pos.beta <- grep("beta", parm.names)
  
  mu <- parm[pos.mu]
  sigma <- parm[pos.sigma]
  p <- parm[pos.p]
  beta <- parm[pos.beta]
  
  if (K == 2) {
    if (length(grep("gl", mixture)) == 1) {
      p[3] <- 0
    }
    if (length(grep("gw", mixture)) == 1) {
      p[3] <- p[2]
      p[2] <- 0
    }
    if (length(grep("lw", mixture)) == 1) {
      p[2:3] <- p[1:2]
      p[1] <- 0
    }
    names(p) <- as.parm.names(list(p=rep(0,3)))
  }
  
  xbeta <- X%*%as.matrix(beta)
  eta <- exp(xbeta)
  
  a1 <- D*(log(dglw(y, p, mu, sigma)) + xbeta)
  a2 <- -eta*pglw(y, p, mu, sigma)
  
  LL <- a1 + a2
  
  return(matrix(LL, nrow = 1))
}


loglik_ptcr.g <- function (parm, Data) {
  
  if (is.null(Data$X)) Data$X <- matrix(0, nrow = length(Data$y), ncol = 1)
  
  y <- Data$y
  D <- Data$d
  X <- as.matrix(Data$X)
  J <- ncol(X)
  
  parm.names <- as.parm.names(list(mu=0, sigma=0, beta=rep(0, J)))
  
  pos.mu <- grep("mu", parm.names)
  pos.sigma <- grep("sigma", parm.names)
  pos.beta <- grep("beta", parm.names)
  
  mu <- parm[pos.mu]
  sigma <- parm[pos.sigma]
  beta <- parm[pos.beta]
  
  xbeta <- X%*%as.matrix(beta)
  eta <- exp(xbeta)
  
  al <- mu^2/sigma
  be  <- sigma/mu
  
  a1 <- D*(dgamma(y, shape=al, scale=be, log=T) + xbeta)
  a2 <- -eta*pgamma(y, shape=al, scale=be)
  
  LL <- a1 + a2
  
  return(matrix(LL, nrow = 1))
}

loglik_ptcr.l <- function (parm, Data) {
  
  if (is.null(Data$X)) Data$X <- matrix(0, nrow = length(Data$y), ncol = 1)
  
  y <- Data$y
  D <- Data$d
  X <- as.matrix(Data$X)
  J <- ncol(X)
  
  parm.names <- as.parm.names(list(mu=0, sigma=0, beta=rep(0, J)))
  
  pos.mu <- grep("mu", parm.names)
  pos.sigma <- grep("sigma", parm.names)
  pos.beta <- grep("beta", parm.names)
  
  mu <- parm[pos.mu]
  sigma <- parm[pos.sigma]
  beta <- parm[pos.beta]
  
  xbeta <- X%*%as.matrix(beta)
  eta <- exp(xbeta)
  
  mu_l  <- log(mu^2/sqrt(mu^2+sigma))
  sig_l <- sqrt(log((mu^2+sigma)/mu^2))
  
  a1 <- D*(dlnorm(y, mu_l, sig_l, log=T) + xbeta)
  a2 <- -eta*plnorm(y, mu_l, sig_l)
  
  LL <- a1 + a2
  
  return(matrix(LL, nrow = 1))
}


loglik_ptcr.w <- function (parm, Data) {
  
  if (is.null(Data$X)) Data$X <- matrix(0, nrow = length(Data$y), ncol = 1)
  
  y <- Data$y
  D <- Data$d
  X <- as.matrix(Data$X)
  J <- ncol(X)
  
  parm.names <- as.parm.names(list(mu=0, sigma=0, beta=rep(0, J)))
  
  pos.mu <- grep("mu", parm.names)
  pos.sigma <- grep("sigma", parm.names)
  pos.beta <- grep("beta", parm.names)
  
  mu <- parm[pos.mu]
  sigma <- parm[pos.sigma]
  beta <- parm[pos.beta]
  
  xbeta <- X%*%as.matrix(beta)
  eta <- exp(xbeta)
  
  parms <- Par(mu, sigma)
  
  a1 <- D*(dweibull(y, shape=parms$gama, scale=parms$delta, log=T) + xbeta)
  a2 <- -eta*pweibull(y, shape=parms$gama, scale=parms$delta)
  
  LL <- a1 + a2
  
  return(matrix(LL, nrow = 1))
}
