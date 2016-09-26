#---------------------------------------------------------------------------------------------------------------------#
### Funcao que cria replicas da amostra original baseada em uma observacao da posteriori dos parametros do modelo
### parm: resultado MCMC
### X: Matriz Modelo original
### Resultado: Cada linha da matriz representa um Yrep !!!
#---------------------------------------------------------------------------------------------------------------------#
Yrep <- function(parm, X) {
  
  J <- ncol(parm)-3
  r <- ncol(X)
  N <- nrow(X)
  yrep <- matrix(0, nrow=nrow(parm), ncol=N)
  
  for(i in 1:N) {
    xbeta <- 0
    Betas <- parm[,1:J]
    Sigma <- parm[,(J+1)]
    Pesos <- cbind(0, parm[,(J+2):(J+3)])
    if( J>r ) {
      if(r>1) xbeta <- t(X[i,2:r]%*%t(Betas[,3:J]))
      mu.lnorm <- exp(Betas[,1] + xbeta)
      mu.weibu <- exp(Betas[,2] + xbeta)
      mu.gamma <- mu.lnorm
    } #end if
    if(J==r) {
      if(r>1) xbeta <- X[i,2:r]%*%t(Betas[,2:J])
      mu.gamma <- exp(Betas[,1] + xbeta)
      mu.lnorm <- mu.gamma
      mu.weibu <- mu.gamma
    } #end if
    for(j in 1:nrow(parm)) yrep[j,i] <- rglw2(n=1, P=Pesos[j,], mu.gamma[j], mu.lnorm[j], mu.weibu[j], Sigma=Sigma[j])
  }  # end for i in 1:N
  return(yrep)
}


#' Predictive Distribution of GLW 
#' 
#' \code{pred} estimates the predictive distribution for a new observation
#' 
#' @param parm: matrix of parameters generated from the posterior distribution.
#' @param X: vector containing the new observation covariates.
#' @param mixture: indicates the finite mixture to be fitted to the data. Options are: 'glw', 'gl', 'gw', 'lw'.
#' 
#' @return A generated sample of the predictive distribution for the new observation.
#' 
#' @export
pred <- function(parm, X, mixture) {
  
  if(mixture=='glw') out <- pred.glw(parm=parm,X=X)
  if(mixture=='gl')  out <- pred.gl(parm=parm, X=X)
  if(mixture=='gw')  out <- pred.gw(parm=parm, X=X)
  if(mixture=='lw')  out <- pred.lw(parm=parm, X=X)
  return(out)
}

pred.lw <-function(parm, X){
  
  y <- c()
  J <- ncol(parm)-3
  r <- length(X)
  
  xbeta <- 0
  Betas <- parm[,1:J]
  Sigma <- parm[,(J+1)]
  Pesos <- cbind(0, parm[,(J+2):(J+3)])
  if( J>r ) {
    if(r>1) xbeta <- t(X[2:r]%*%t(Betas[,3:J]))
    mu.lnorm <- exp(Betas[,1] + xbeta)
    mu.weibu <- exp(Betas[,2] + xbeta)
    mu.gamma <- mu.lnorm
  }
  if(J==r) {
    if(r>1) xbeta <- X[2:r]%*%t(Betas[,2:J])
    mu.gamma <- exp(Betas[,1] + xbeta)
    mu.lnorm <- mu.gamma
    mu.weibu <- mu.gamma
  }
  
  for(j in 1:nrow(Betas)) y[j] <- rglw2(n=1, P=Pesos[j,], mu.gamma[j], mu.lnorm[j], mu.weibu[j], Sigma=Sigma[j])
  return(y)
}

pred.gw <-function(parm, X){
  
  y  <- c()
  J <- ncol(parm)-3
  r <- length(X)
  
  xbeta <- 0
  Betas <- parm[,1:J]
  Sigma <- parm[,(J+1)]
  Pesos <- cbind(parm[,(J+2)], 0, parm[,(J+3)])
  if( J>r ) {
    if(r>1) xbeta <- t(X[2:r]%*%t(Betas[,3:J]))
    mu.gamma <- exp(Betas[,1] + xbeta)
    mu.lnorm <- mu.gamma
    mu.weibu <- exp(Betas[,2] + xbeta)
  }
  if(J==r) {
    if(r>1) xbeta <- X[2:r]%*%t(Betas[,2:J])
    mu.gamma <- exp(Betas[,1] + xbeta)
    mu.lnorm <- mu.gamma
    mu.weibu <- mu.gamma
  }
  
  for(j in 1:nrow(Betas)) y[j] <- rglw2(n=1, P=Pesos[j,], mu.gamma[j], mu.lnorm[j], mu.weibu[j], Sigma=Sigma[j])
  return(y)
}


pred.gl <-function(parm, X){
  
  y  <- c()
  J <- ncol(parm)-3
  r <- length(X)
  
  xbeta <- 0
  Betas <- parm[,1:J]
  Sigma <- parm[,(J+1)]
  Pesos <- cbind(parm[,(J+2):(J+3)],0)
  if( J>r ) {
    if(r>1) xbeta <- t(X[2:r]%*%t(Betas[,3:J]))
    mu.gamma <- exp(Betas[,1] + xbeta)
    mu.lnorm <- exp(Betas[,2] + xbeta)
    mu.weibu <- mu.gamma
  }
  if(J==r) {
    if(r>1) xbeta <- X[2:r]%*%t(Betas[,2:J])
    mu.gamma <- exp(Betas[,1] + xbeta)
    mu.lnorm <- mu.gamma
    mu.weibu <- mu.gamma
  }
  
  for(j in 1:nrow(Betas)) y[j] <- rglw2(n=1, P=Pesos[j,], mu.gamma[j], mu.lnorm[j], mu.weibu[j], Sigma=Sigma[j])
  return(y)
}


pred.glw <-function(parm, X){
  
  y  <- c()
  J <- ncol(parm)-4
  r <- length(X)
  
  xbeta <- 0
  Betas <- parm[,1:J]
  Sigma <- parm[,(J+1)]
  Pesos <- parm[,(J+2):(J+4)]
  if( J>r ) {
    if(r>1) xbeta <- t(X[2:r]%*%t(Betas[,4:J]))
    mu.gamma <- exp(Betas[,1] + xbeta)
    mu.lnorm <- exp(Betas[,2] + xbeta)
    mu.weibu <- exp(Betas[,3] + xbeta)
  }
  if(J==r) {
    if(r>1) xbeta <- X[2:r]%*%t(Betas[,2:J])
    mu.gamma <- exp(Betas[,1] + xbeta)
    mu.lnorm <- mu.gamma
    mu.weibu <- mu.gamma
  }
  
  for(j in 1:nrow(Betas)) y[j] <- rglw2(n=1, P=Pesos[j,], mu.gamma[j], mu.lnorm[j], mu.weibu[j], Sigma=Sigma[j])
  return(y)
}




