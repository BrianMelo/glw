#---------------------------------------------------------------------#
#------CPO - Conditional Predictive Ordinate--------------------------#
#---------------------------------------------------------------------#
#'CPO of the GLW FM Standard Cure Rate Model
#'
#' \code{cpo_scr} estimates the CPO based on a MCMC sample of the posterior distribuiton.
#' 
#' @param parm: sample of posterior parameters.
#' @param Data: the same data as in \code{\link{glw_scr}}.
#' @param mixture: indicates the finite mixture to be fitted to the data. Options are: 'glw', 'gl', 'gw', 'lw'.
#' 
#' @return Returns the Conditional Predictive Ordinate of the GLW finite mixture model for every observation.
#' 
#' @export
cpo_scr <- function (parm, Data, mixture) {
  cpo <- matrix(nrow = nrow(parm), ncol = nrow(Data$X))
  if(nchar(mixture)>1) {
    for (i in 1:nrow(parm)) cpo[i, ] <- 1/exp(loglik_scr(parm = parm[i,], Data = Data, mixture = mixture))
  }
  if(nchar(mixture)==1) {
    for (i in 1:nrow(parm)) {
      aux <- eval(parse(text=paste('loglik_scr.', mixture, '(parm = parm[i,], Data = Data)', sep='')))
      cpo[i, ] <- 1/exp(aux)    
  }}
  out <- 1/colMeans(cpo)
  return(out)
}

loglik_scr <- function (parm, Data, mixture) {
  
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
  pi <- exp(xbeta)/(1+exp(xbeta))
  
  a1 <- ifelse(D==1, (1-pi)*dglw(y, p, mu, sigma),1)
  a2 <- ifelse(D==0, pi + (1-pi)*pglw(y, p, mu, sigma, lower.tail=F),1)
  
  LL <- log(a1) + log(a2)
  
  return(matrix(LL, nrow = 1))
}


loglik_scr.g <- function (parm, Data) {
  
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
  pi <- exp(xbeta)/(1+exp(xbeta))
  
  al <- mu^2/sigma
  be  <- sigma/mu
  
  a1 <- ifelse(D==1, (1-pi)*dgamma(y, shape=al, scale=be),1)
  a2 <- ifelse(D==0, pi + (1-pi)*pgamma(y, shape=al, scale=be, lower.tail=F),1)
  
  LL <- log(a1) + log(a2)
  
  return(matrix(LL, nrow = 1))
}

loglik_scr.l <- function (parm, Data) {
  
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
  pi <- exp(xbeta)/(1+exp(xbeta))
  
  mu_l  <- log(mu^2/sqrt(mu^2+sigma))
  sig_l <- sqrt(log((mu^2+sigma)/mu^2))
  
  a1 <- ifelse(D==1, (1-pi)*dlnorm(y, mu_l, sig_l),1)
  a2 <- ifelse(D==0, pi + (1-pi)*plnorm(y, mu_l, sig_l, lower.tail=F),1)
  
  LL <- log(a1) + log(a2)
  
  return(matrix(LL, nrow = 1))
}

loglik_scr.w <- function (parm, Data) {
  
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
  pi <- exp(xbeta)/(1+exp(xbeta))
  
  parms <- Par(mu, sigma)
  
  a1 <- ifelse(D==1, (1-pi)*dweibull(y, shape=parms$gama, scale=parms$delta),1)
  a2 <- ifelse(D==0, pi + (1-pi)*pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=F),1)
  
  LL <- log(a1) + log(a2)
  
  return(matrix(LL, nrow = 1))
}
