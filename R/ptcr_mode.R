#-----------------------------------------------------------------------------------------#
#-----Computes the posterior mode---------------------------------------------------------#
#-----------------------------------------------------------------------------------------#
#'Computes the Posterior Mode
#'
#'\code{mode_pt} computes the posterior mode of the GLW FM Promotion Time Cure Rate model
#'
#' @param Data: data as in \code{\link{glwfm}}.
#' @param prior: priori distributions as in \code{\link{glwfm}}.
#' @param mixture: indicates the finite mixture to be fitted to the data. Options are: 'glw', 'gl', 'gw', 'lw', 'g', 'l', 'w'
#' 
#' @return Returns a list containing the mode of the posterior distribution of the GLW finite mixture regression model, 
#'      the log-posterior maximum, the convergence status (as in optim) and a message (as in optim).
#'      
#' @export
mode_pt <- function (Data, prior, mixture) {
  if (nchar(mixture) > 1) 
    out <- mode_pt_mix(Data = Data, prior = prior, mixture = mixture)
  if (nchar(mixture) == 1) 
    out <- mode_pt_single(Data = Data, prior = prior, mixture = mixture)
  return(out)
}

#Funcao que retorna a moda a posteriori da mistura GLW modelo mistura padrão, usada tambem para o FBST 
#Data: base de dados incluindo vetor initial.par de valores iniciais para a otimizacao (optim)
#prior: distribuicoes a priori
mode_pt_mix <- function(Data, prior, mixture){
  Model <- function (parm, Data, prior, mixture) {
    
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
    
    mu <- exp(parm[pos.mu])
    sigma <- exp(parm[pos.sigma])
    p <- exp(parm[pos.p])/(1+exp(parm[pos.p]))
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
    
    a1 <- ifelse(D==1, (log(dglw(y, p, mu, sigma)) + xbeta), 0)
    a2 <- -eta*pglw(y, p, mu, sigma)
    
    LL <- sum(a1 + a2)
    
    ### Log-Prior
    mu.prior    <- ifelse(is.null(prior$mu),    dgamma(mu, 0.01, 0.01, log=TRUE),    sum(eval(parse(1, text=prior$mu))))
    sigma.prior <- ifelse(is.null(prior$sigma), dgamma(sigma, 0.01, 0.01, log=TRUE), sum(eval(parse(1, text=prior$sigma))) )
    p.prior     <- ifelse(is.null(prior$alpha), ifelse(K<3, 0, ddirichlet(p, c(1,1,1), log=TRUE)),   ddirichlet(p, prior$alpha, log=TRUE))
    beta.prior  <- ifelse(is.null(prior$beta),  sum(dnormv(beta, 0, 1000, log=TRUE)), sum(eval(parse(1, text=prior$beta))))
    
    ### Log-Posterior
    LP <- LL + sum(beta.prior) + sigma.prior + p.prior + mu.prior
    
    return(LP/nrow(X))
  }
  K <- nchar(mixture)
  r <- ncol(as.matrix(Data$X))
  ini <- c(log(Data$initial.par[1:(2+K)]), Data$initial.par[(3+K):(3+K+r-1)])
  op <- optim(par=ini, fn=Model, Data=Data, prior=prior, mixture=mixture, control=list(fnscale=-1, reltol=1e-10, maxit=10**4) )
  out <- list(Posterior.mode=round(c(exp(op$par[1:2]), exp(op$par[3:(3+K-1)])/sum(exp(op$par[3:(3+K-1)])), op$par[(3+K):(3+K+r-1)]), 4), 
              value=op$value*length(Data$y), Convergence=op$convergence, message=op$message)
  return(out)
}

mode_pt_single <- function(Data, prior, mixture){
  Model <- function (parm, Data, prior, mixture) {
    
    if (is.null(Data$X)) Data$X <- matrix(0, nrow = length(Data$y), ncol = 1)
    
    y <- Data$y
    D <- Data$d
    X <- as.matrix(Data$X)
    J <- ncol(X)
    
    parm.names <- as.parm.names(list(mu=0, sigma=0, beta=rep(0, J)))
    
    pos.mu <- grep("mu", parm.names)
    pos.sigma <- grep("sigma", parm.names)
    pos.beta <- grep("beta", parm.names)
    
    mu <- exp(parm[pos.mu])
    sigma <- exp(parm[pos.sigma])
    beta <- parm[pos.beta]
    
    xbeta <- X%*%as.matrix(beta)
    eta <- exp(xbeta)
    
    if(mixture=='g') p <- c(1, 0, 0)
    if(mixture=='l') p <- c(0, 1, 0)
    if(mixture=='w') p <- c(0, 0, 1)
    
    a1 <- ifelse(D==1, (log(dglw(y, p, mu, sigma)) + xbeta), 0)
    a2 <- -eta*pglw(y, p, mu, sigma)
    
    LL <- sum(a1 + a2)
    
    ### Log-Prior
    mu.prior    <- ifelse(is.null(prior$mu),    dgamma(mu, 0.01, 0.01, log=TRUE),    sum(eval(parse(1, text=prior$mu))))
    sigma.prior <- ifelse(is.null(prior$sigma), dgamma(sigma, 0.01, 0.01, log=TRUE), sum(eval(parse(1, text=prior$sigma))) )
    beta.prior  <- ifelse(is.null(prior$beta),  sum(dnormv(beta, 0, 1000, log=TRUE)), sum(eval(parse(1, text=prior$beta))))
    
    ### Log-Posterior
    LP <- LL + sum(beta.prior) + sigma.prior + mu.prior
    
    return(LP/nrow(X))
  } 
  r <- ncol(as.matrix(Data$X))
  ini <- c(log(Data$initial.par[1:2]), Data$initial.par[3:(3+r-1)])
  op <- optim(par=ini, fn=Model, Data=Data, prior=prior, mixture=mixture, control=list(fnscale=-1, reltol=1e-10, maxit=10**4) )
  out <- list(Posterior.mode=round(c(exp(op$par[1:2]), op$par[3:(3+r-1)]),4), value=op$value*length(Data$y), Convergence=op$convergence, message=op$message)
  return(out)
}