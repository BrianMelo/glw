#-----------------------------------------------------------------------------------------#
#-----Computes the posterior mode---------------------------------------------------------#
#-----------------------------------------------------------------------------------------#
#'Computes the Posterior Mode
#'
#'\code{Mode_glw} computes the posterior mode of the GLW finite mixture model
#'
#' @param Data: data as in \code{\link{glwfm}}.
#' @param prior: priori distributions as in \code{\link{glwfm}}.
#' @param mixture: indicates the finite mixture to be fitted to the data. Options are: 'glw', 'gl', 'gw', 'lw', 'g', 'l', 'w'
#' 
#' @return Returns a list containing the mode of the posterior distribution of the GLW finite mixture regression model, 
#'      the log-posterior maximum, the convergence status (as in optim) and a message (as in optim).
#'      
#' @export
mode_glw <- function(Data, prior, mixture) {
  
  if(mixture=='glw') out <- Mode.glw(Data=Data, prior=prior)
  
  if(mixture=='gl')  out <- Mode.gl(Data=Data, prior=prior)
  if(mixture=='gw')  out <- Mode.gw(Data=Data, prior=prior)
  if(mixture=='lw')  out <- Mode.lw(Data=Data, prior=prior)
  
  if(mixture=='g') out <- Mode.gamma(Data=Data, prior=prior)
  if(mixture=='l') out <- Mode.lnorm(Data=Data, prior=prior)
  if(mixture=='w') out <- Mode.weibull(Data=Data, prior=prior)
  
  return(out)
}

#Funcao que retorna a moda a posteriori da mistura GLW, (usada tambem para o FBST)
#Data: base de dados incluindo vetor initial.par de valores iniciais para a otimizacao (optim)
#prior: distribuicoes a priori
Mode.glw <- function(Data, prior) {
  
  J <- ifelse(Data$rint==T, 2, 0) + ncol(Data$X)
  
  Model <- function(parm, Data, prior) {
    y <- Data$y
    yf <- Data$yf
    D  <- Data$d
    X <- Data$X
    rint <- Data$rint
    
    J <- ifelse(rint==T, 2, 0)+ncol(X)
    parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0, p=rep(0,3)))
    
    pos.beta  <- grep("beta", parm.names)
    pos.p     <- grep("p", parm.names)
    pos.sigma <- grep("sigma", parm.names)      
    
    ### Parameters
    beta <- parm[pos.beta]
    p <- exp(parm[pos.p])/sum(exp(parm[pos.p]))
    sigma <- exp(parm[pos.sigma])
    
    ### Log-Prior
    beta.prior  <- ifelse(is.null(prior$beta), sum(dnormv(beta, 0, 1000, log=TRUE)), sum(eval(parse(1, text=prior$beta))) )
    sigma.prior <- ifelse(is.null(prior$sigma), sum(dexp(sigma, 1, log=TRUE)), sum(eval(parse(1, text=prior$sigma))) )
    p.prior     <- ifelse(is.null(prior$alpha), ddirichlet(p, c(1,1,1), log=TRUE), ddirichlet(p, prior$alpha, log=TRUE) )
    
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
    
    LL <- sum(log(a1) + log(a2) + log(a3) + log(a4) )
    
    ### Log-Posterior
    LP <- LL + sum(beta.prior) + sigma.prior + p.prior
    
    return(as.numeric(LP/length(Data$y)))
  }
  ini <- c(Data$initial.par[1:J], log(Data$initial.par[(J+1):(J+4)]))
  op <- optim(par=ini, fn=Model, Data=Data, prior=prior, control=list(fnscale=-1, reltol=1e-10, maxit=10**4) )
  out <- list(Posterior.mode=c(op$par[1:J], exp(op$par[J+1]), exp(op$par[(J+2):(J+4)])/sum(exp(op$par[(J+2):(J+4)]))), 
              value=op$value*length(Data$y), Convergence=op$convergence, message=op$message)
  return(out)
}

Mode.lw <- function(Data, prior) {
  J <- ifelse(Data$rint==T, 1, 0) + ncol(Data$X)
  
  Modellw <- function(parm, Data, prior) {
    
    y <- Data$y
    yf <- Data$yf
    D  <- Data$d
    X <- Data$X
    rint <- Data$rint
    
    J <- ifelse(rint==T, 1, 0) + ncol(X)
    parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0, p=rep(0,2)))
    
    pos.beta  <- grep("beta", parm.names)
    pos.p     <- grep("p", parm.names)
    pos.sigma <- grep("sigma", parm.names)      
    
    ### Parameters
    beta <- parm[pos.beta]
    p <- exp(parm[pos.p])/sum(exp(parm[pos.p]))
    sigma <- exp(parm[pos.sigma])
    
    ### Log-Prior
    beta.prior  <- ifelse(is.null(prior$beta), sum(dnormv(beta, 0, 1000, log=TRUE)), sum(eval(parse(1, text=prior$beta))) )
    sigma.prior <- ifelse(is.null(prior$sigma), sum(dexp(sigma, 1, log=TRUE)), sum(eval(parse(1, text=prior$sigma))) )
    p.prior     <- ifelse(is.null(prior$alpha), ddirichlet(p, c(1,1), log=TRUE), ddirichlet(p, prior$alpha, log=TRUE) )
    
    ### Log-Likelihood
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
    
    a1 <- ifelse(D==1, p[1]*dlnorm(y, parms$mu_l, parms$sig_l) + 
                   p[2]*dweibull(y, shape=parms$gama, scale=parms$delta),1)
    
    a2 <- ifelse(D==0, p[1]*plnorm(y, parms$mu_l, parms$sig_l, lower.tail=F) + 
                   p[2]*pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=F),1)
    
    a3 <- ifelse(D==2, p[1]*plnorm(y, parms$mu_l, parms$sig_l, lower.tail=T) + 
                   p[2]*pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T),1)
    
    a4 <- ifelse(D==3, p[1]*(plnorm(y, parms$mu_l, parms$sig_l, lower.tail=T) - plnorm(y, parms$mu_l, parms$sig_l, lower.tail=T)) +
                   p[2]*(pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T) - pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T)),1)
    
    LL <- sum(log(a1) + log(a2) + log(a3) + log(a4) )
    
    ### Log-Posterior
    LP <- LL + sum(beta.prior) + sigma.prior + p.prior
    return(LP/length(y))
  }
  ini <- c(Data$initial.par[1:J], log(Data$initial.par[(J+1):(J+3)]))
  op <- optim(par=ini, fn=Modellw, Data=Data, prior=prior, control=list(fnscale=-1, reltol=1e-10, maxit=10**4) )
  out <- list(Posterior.mode=c(op$par[1:J], exp(op$par[J+1]), exp(op$par[(J+2):(J+3)])/sum(exp(op$par[(J+2):(J+3)]))), 
              value=op$value*length(Data$y), Convergence=op$convergence, message=op$message)
  return(out)
}

Mode.gl <- function(Data, prior) {
  J <- ifelse(Data$rint==T, 1, 0) + ncol(Data$X)
  
  Modelgl <- function(parm, Data, prior) {
    
    y <- Data$y
    yf <- Data$yf
    D  <- Data$d
    X <- Data$X
    rint <- Data$rint
    
    J <- ifelse(rint==T, 1, 0) + ncol(X)
    parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0, p=rep(0,2)))
    
    pos.beta  <- grep("beta", parm.names)
    pos.p     <- grep("p", parm.names)
    pos.sigma <- grep("sigma", parm.names)      
    
    ### Parameters
    beta <- parm[pos.beta]
    p <- exp(parm[pos.p])/sum(exp(parm[pos.p]))
    sigma <- exp(parm[pos.sigma])
    
    ### Log-Prior
    beta.prior  <- ifelse(is.null(prior$beta), sum(dnormv(beta, 0, 1000, log=TRUE)), sum(eval(parse(1, text=prior$beta))) )
    sigma.prior <- ifelse(is.null(prior$sigma), sum(dexp(sigma, 1, log=TRUE)), sum(eval(parse(1, text=prior$sigma))) )
    p.prior     <- ifelse(is.null(prior$alpha), ddirichlet(p, c(1,1), log=TRUE), ddirichlet(p, prior$alpha, log=TRUE) )
    
    ### Log-Likelihood
    xbeta <- 0
    if(Data$rint==T){
      if(ncol(X)>1) xbeta <- X[,2:ncol(X)]%*%as.matrix(beta[3:J])
      mu.gamma <- exp(beta[1] + xbeta)
      mu.lnorm <- exp(beta[2] + xbeta)
      mu.weibu <- exp(xbeta)
    }
    if(Data$rint==F) {
      if(ncol(X)>1) xbeta <- X[,2:ncol(X)]%*%as.matrix(beta[2:J])
      mu.gamma <- exp(beta[1] + xbeta)
      mu.lnorm <- mu.gamma
      mu.weibu <- mu.gamma
    }
    parms <- Parm3(mu.gamma, mu.lnorm, mu.weibu, sigma)
    
    a1 <- ifelse(D==1, p[1]*dgamma(y, shape=parms$alpha, scale=parms$beta) + 
                   p[2]*dlnorm(y, parms$mu_l, parms$sig_l),1)
    
    a2 <- ifelse(D==0, p[1]*pgamma(y, shape=parms$alpha, scale=parms$beta, lower.tail=F) + 
                   p[2]*plnorm(y, parms$mu_l, parms$sig_l, lower.tail=F),1)
    
    a3 <- ifelse(D==2, p[1]*pgamma(y, shape=parms$alpha, scale=parms$beta, lower.tail=T) + 
                   p[2]*plnorm(y, parms$mu_l, parms$sig_l, lower.tail=T),1)
    
    a4 <- ifelse(D==3, p[1]*(pgamma(y, shape=parms$alpha, scale=parms$beta, lower.tail=T) - pgamma(yf, shape=parms$alpha, scale=parms$beta, lower.tail=T)) + 
                   p[2]*(plnorm(y, parms$mu_l, parms$sig_l, lower.tail=T) - plnorm(y, parms$mu_l, parms$sig_l, lower.tail=T)),1)
    
    LL <- sum(log(a1) + log(a2) + log(a3) + log(a4) )
    
    ### Log-Posterior
    LP <- LL + sum(beta.prior) + sigma.prior + p.prior
    return(LP/length(y))
  }
  ini <- c(Data$initial.par[1:J], log(Data$initial.par[(J+1):(J+3)]))
  op <- optim(par=ini, fn=Modelgl, Data=Data, prior=prior, control=list(fnscale=-1, reltol=1e-10, maxit=10**4) )
  out <- list(Posterior.mode=c(op$par[1:J], exp(op$par[J+1]), exp(op$par[(J+2):(J+3)])/sum(exp(op$par[(J+2):(J+3)]))), 
              value=op$value*length(Data$y), Convergence=op$convergence, message=op$message)
  return(out)
}


Mode.gw <- function(Data, prior) {
  J <- ifelse(Data$rint==T, 1, 0) + ncol(Data$X)
  
  Modelgw <- function(parm, Data, prior) {
    
    y <- Data$y
    yf <- Data$yf
    D  <- Data$d
    X <- Data$X
    rint <- Data$rint
    
    J <- ifelse(rint==T, 1, 0) + ncol(X)
    parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0, p=rep(0,2)))
    
    pos.beta  <- grep("beta", parm.names)
    pos.p     <- grep("p", parm.names)
    pos.sigma <- grep("sigma", parm.names)      
    
    ### Parameters
    beta <- parm[pos.beta]
    p <- exp(parm[pos.p])/sum(exp(parm[pos.p]))
    sigma <- exp(parm[pos.sigma])
    
    ### Log-Prior
    beta.prior  <- ifelse(is.null(prior$beta), sum(dnormv(beta, 0, 1000, log=TRUE)), sum(eval(parse(1, text=prior$beta))) )
    sigma.prior <- ifelse(is.null(prior$sigma), sum(dexp(sigma, 1, log=TRUE)), sum(eval(parse(1, text=prior$sigma))) )
    p.prior     <- ifelse(is.null(prior$alpha), ddirichlet(p, c(1,1), log=TRUE), ddirichlet(p, prior$alpha, log=TRUE) )
    
    ### Log-Likelihood
    xbeta <- 0
    if(Data$rint==T){
      if(ncol(X)>1) xbeta <- X[,2:ncol(X)]%*%as.matrix(beta[3:J])
      mu.gamma <- exp(beta[1] + xbeta)
      mu.lnorm <- exp(xbeta)
      mu.weibu <- exp(beta[2] + xbeta)
    }
    if(Data$rint==F) {
      if(ncol(X)>1) xbeta <- X[,2:ncol(X)]%*%as.matrix(beta[2:J])
      mu.gamma <- exp(beta[1] + xbeta)
      mu.lnorm <- mu.gamma
      mu.weibu <- mu.gamma
    }
    parms <- Parm3(mu.gamma, mu.lnorm, mu.weibu, sigma)
    
    a1 <- ifelse(D==1, p[1]*dgamma(y, shape=parms$alpha, scale=parms$beta) + 
                   p[2]*dweibull(y, shape=parms$gama, scale=parms$delta),1)
    
    a2 <- ifelse(D==0, p[1]*pgamma(y, shape=parms$alpha, scale=parms$beta, lower.tail=F) + 
                   p[2]*pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=F),1)
    
    a3 <- ifelse(D==2, p[1]*pgamma(y, shape=parms$alpha, scale=parms$beta, lower.tail=T) + 
                   p[2]*pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T),1)
    
    a4 <- ifelse(D==3, p[1]*(pgamma(y, shape=parms$alpha, scale=parms$beta, lower.tail=T) - pgamma(yf, shape=parms$alpha, scale=parms$beta, lower.tail=T)) + 
                   p[2]*(pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T) - pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T)),1)
    
    LL <- sum(log(a1) + log(a2) + log(a3) + log(a4) )
    
    ### Log-Posterior
    LP <- LL + sum(beta.prior) + sigma.prior + p.prior
    return(LP/length(y))
  }
  ini <- c(Data$initial.par[1:J], log(Data$initial.par[(J+1):(J+3)]))
  op <- optim(par=ini, fn=Modelgw, Data=Data, prior=prior, control=list(fnscale=-1, reltol=1e-10, maxit=10**4) )
  out <- list(Posterior.mode=c(op$par[1:J], exp(op$par[J+1]), exp(op$par[(J+2):(J+3)])/sum(exp(op$par[(J+2):(J+3)]))), 
              value=op$value*length(Data$y), Convergence=op$convergence, message=op$message)
  return(out)
}

Mode.weibull <- function(Data, prior) {
  
  J <- ncol(as.matrix(Data$X))
  
  Modelweib <- function(parm, Data, prior) {
    
    X <- as.matrix(Data$X)
    y <- Data$y
    D <- Data$d
    yf <- Data$yf
    J  <- ncol(X)
    
    parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
    
    pos.beta  <- grep("beta", parm.names)
    pos.sigma <- grep("sigma", parm.names)
    
    ### Parameters
    beta <- parm[pos.beta]
    sigma <- exp(parm[pos.sigma])
    
    ### Log-Prior
    beta.prior  <- ifelse(is.null(prior$beta), sum(dnormv(beta, 0, 1000, log=TRUE)), sum(eval(parse(1, text=prior$beta))) )
    sigma.prior <- ifelse(is.null(prior$sigma), sum(dexp(sigma, 1, log=TRUE)), sum(eval(parse(1, text=prior$sigma))) )
    
    ### Log-Likelihood
    mu <- exp(X%*%as.matrix(beta))
    parms <- Parm3(mu, mu, mu, sigma)
    
    a1 <- ifelse(D==1, dweibull(y, shape=parms$gama, scale=parms$delta), 1)
    a2 <- ifelse(D==0, pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=F), 1)
    a3 <- ifelse(D==2, pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T), 1)
    a4 <- ifelse(D==3, pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T) - pweibull(y, shape=parms$gama, scale=parms$delta, lower.tail=T), 1)
    
    LL <- sum(log(a1) + log(a2) + log(a3) + log(a4) )
    
    ### Log-Posterior
    LP <- LL + sum(beta.prior) + sigma.prior
    return(LP/length(y))
  }
  ini <- c(Data$initial.par[1:J], log(Data$initial.par[J+1]))
  op <- optim(par=ini, fn=Modelweib, Data=Data, prior=prior, control=list(fnscale=-1, reltol=1e-10, maxit=10**4) )
  out <- list(Posterior.mode=c(op$par[1:J], exp(op$par[J+1]) ),value=op$value*length(Data$y), Convergence=op$convergence, message=op$message)
  return(out)
}

Mode.gamma <- function(Data, prior) {
  
  J <- ncol(as.matrix(Data$X))
  
  Modelgama <- function(parm, Data, prior) {
    
    X <- as.matrix(Data$X)
    y <- Data$y
    D <- Data$d
    yf <- Data$yf
    J  <- ncol(X)
    
    parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
    
    pos.beta  <- grep("beta", parm.names)
    pos.sigma <- grep("sigma", parm.names)
    
    ### Parameters
    beta <- parm[pos.beta]
    sigma <- exp(parm[pos.sigma])
    
    ### Log-Prior
    beta.prior  <- ifelse(is.null(prior$beta), sum(dnormv(beta, 0, 1000, log=TRUE)), sum(eval(parse(1, text=prior$beta))) )
    sigma.prior <- ifelse(is.null(prior$sigma), sum(dexp(sigma, 1, log=TRUE)), sum(eval(parse(1, text=prior$sigma))) )
    
    ### Log-Likelihood
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
    return(LP/length(y))
  }
  ini <- c(Data$initial.par[1:J], log(Data$initial.par[J+1]))
  op <- optim(par=ini, fn=Modelgama, Data=Data, prior=prior, control=list(fnscale=-1, reltol=1e-10, maxit=10**4) )
  out <- list(Posterior.mode=c(op$par[1:J], exp(op$par[J+1]) ),value=op$value*length(Data$y), Convergence=op$convergence, message=op$message)
  return(out)
}

Mode.lnorm <- function(Data, prior) {
  
  J <- ncol(as.matrix(Data$X))
  
  Modellnorm <- function(parm, Data, prior) {
    
    X <- as.matrix(Data$X)
    y <- Data$y
    D <- Data$d
    yf <- Data$yf
    J  <- ncol(X)
    
    parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
    
    pos.beta  <- grep("beta", parm.names)
    pos.sigma <- grep("sigma", parm.names)
    
    ### Parameters
    beta <- parm[pos.beta]
    sigma <- exp(parm[pos.sigma])
    
    ### Log-Prior
    beta.prior  <- ifelse(is.null(prior$beta), sum(dnormv(beta, 0, 1000, log=TRUE)), sum(eval(parse(1, text=prior$beta))) )
    sigma.prior <- ifelse(is.null(prior$sigma), sum(dexp(sigma, 1, log=TRUE)), sum(eval(parse(1, text=prior$sigma))) )
    
    ### Log-Likelihood
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
    return(LP/length(y))
  }
  ini <- c(Data$initial.par[1:J], log(Data$initial.par[J+1]))
  op <- optim(par=ini, fn=Modellnorm, Data=Data, prior=prior, control=list(fnscale=-1, reltol=1e-10, maxit=10**4) )
  out <- list(Posterior.mode=c(op$par[1:J], exp(op$par[J+1]) ),value=op$value*length(Data$y), Convergence=op$convergence, message=op$message)
  return(out)
}