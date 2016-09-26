#' Generate Random Sample
#' 
#' Generates a random sample from the GLW finite mixture with different means for the different distributions.
#' 
#' @param n: number of observations.
#' @param P: weights of the finite mixture.
#' @param Mu.gamma,Mu.lnorm,Mu.weibull: means of the Gamma, Lognormal and Weibull distribution.
#' @param Sigma: variance 
#' 
#' @export
rglw2 <- function(n, P, Mu.gamma, Mu.lnorm, Mu.weibull, Sigma) {
  
  a <- Parm3(Mu.gamma, Mu.lnorm, Mu.weibull, Sigma)
  z <- rmultinom(n, 1, P)
  d <- z[1,]*rgamma(n, shape=a$alpha, scale=a$beta) + z[2,]*rlnorm(n, a$mu_l, a$sig_l) + z[3,]*rweibull(n, shape=a$gama, scale=a$delta)
  return(d)
}

#----------------------------------------------------------#
#-----Computes the right censored times paramater--------- #
#----------------------------------------------------------#
Tau <- function(PC, p, Mu, S2) {
  
  Pc <- function(tau, PC, p, Mu, S2) {
    
    aux <- Par(Mu, S2)
    
    al <- aux$alpha
    be <- aux$beta
    ga <- aux$gama
    de <- aux$delta
    ml <- aux$mu_l
    sl <- aux$sig_l
    
    f <- integrate(function(c) dexp(c, 1/tau)*(p[1]*pgamma(c, shape=al, scale=be) + p[2]*plnorm(c, ml, sl) + p[3]*pweibull(c, shape=ga, scale=de)) , lower=0, upper=Inf)$value-(1-PC)
    return(f)
  }
  
  g = uniroot(Pc, PC=PC, p=p, Mu=Mu, S2=S2, lower=10^-10, upper=10^4)$root
  return(g)
}

#----------------------------------------------------------#
#-----Computes the left censored times paramater---------- #
#----------------------------------------------------------#
Phi <- function(PC, p, Mu, S2) {
  
  Pc <- function(phi, PC, p, Mu, S2) {
    
    aux <- Par(Mu, S2)
    
    al <- aux$alpha
    be <- aux$beta
    ga <- aux$gama
    de <- aux$delta
    ml <- aux$mu_l
    sl <- aux$sig_l
    
    f <- integrate(function(c) dexp(c, 1/phi)*(p[1]*pgamma(c, shape=al, scale=be) + p[2]*plnorm(c, ml, sl) + p[3]*pweibull(c, shape=ga, scale=de)) , lower=0, upper=Inf)$value-PC
    return(f)
  }
  
  g <- uniroot(Pc, PC=PC, p=p, Mu=Mu, S2=S2, lower=10^-10, upper=10^4)$root
  return(g)
}

#-----------------------------------------------------------------------------------------------------------#
#-----Generates a sample of the GLW finite mixture with right or left censored times and covariates---------#
#-----------------------------------------------------------------------------------------------------------#
#' Random Generation From The GLW FM Regression Model
#' 
#' \code{rglw3} generates a sample from the glw finite mixture regression model 
#' with right OR left censored times given a model matrix and the regression coefficients.
#' 
#' @param n: number of observations.
#' @param Betas: regression coefficients.
#' @param X: matrix of covariates (with column of 1's).
#' @param P: weights of the finite mixture.
#' @param Sigma: variance.
#' @param Pc: vector with probabiliity of right and left censored times.
#' 
#' @return A matrix with first column being the times, the second column the censorship 
#' indicator and the remaining columns are the covariates matrix \code{X}.
#' 
#' @export
rglw3 <- function(n, Betas, X, P, Sigma, Pc=c(0,0)) {
  
  k <- length(Betas)
  Mu <- exp(X%*%Betas)
  aux <- Par3(Mu, Sigma)
  Alpha <- aux$alpha ; Beta  <- aux$beta ; Gama  <- aux$gama
  Delta <- aux$delta ; Mu_l  <- aux$mu_l ; Sig_l <- aux$sig_l
  
  if(Pc[1] > 0) {  #Right censored times
    tau   <- rep(0,n)
    for(i in 1:n) tau[i] <- 1/Tau(Pc[1], P, Mu[i], Sigma)
    C <- rexp(n, tau)
  }
  
  if(Pc[2] > 0) {	#Left censored times
    phi   <- rep(0,n)
    for(i in 1:n) phi[i] <- 1/Phi(Pc[2], P, Mu[i], Sigma)
    E <- rexp(n, phi)
  }
  
  Z <- rmultinom(n, 1, P)
  
  T <- Z[1,]*rgamma(n, shape=Alpha, scale=Beta) + Z[2,]*rlnorm(n, Mu_l, Sig_l) + Z[3,]*rweibull(n, shape=Gama, scale=Delta) #Tempo de sobrevida simulado
  
  if(Pc[1]==0) C <- rep(10**6,n)
  if(Pc[2]==0) E <- rep(0,n)
  
  d <- rep(0,n)
  Y <- apply(cbind(apply(cbind(C,T),1,min),E), 1, max) #Observed time
  d[Y==T] <- 1 #Failure
  d[Y==C] <- 0 #Right Censored	
  d[Y==E] <- 2 #Left Censored
  
  return(as.matrix(cbind(Y, d, X)))
}


#-----------------------------------------------------------------------------------------------------------#
#-----Generates a sample of the GLW Standard Cure Rate Model with covariates on the cure rate parameter-----#
#-----------------------------------------------------------------------------------------------------------#
#' Random Generation From The GLW FM Standard Cure Rate Model
#' 
#' \code{rscr_glw} generates a sample from the glw fm standard cure rate regression model 
#' with right censored times given a model matrix and the regression coefficients.
#' 
#' @param n: number of observations.
#' @param Mu: mean
#' @param Sigma: variance.
#' @param P: weights of the finite mixture.
#' @param Betas: regression coefficients associated with the long term parameter.
#' @param X: matrix of covariates (with column of 1's).
#' @param pc: vector with probability of right censored times.
#' 
#' @return A matrix with first column being the times, the second column the censorship 
#' indicator and the remaining columns are the covariates matrix \code{X}.
#' 
#' @export
rscr_glw <- function(n, Mu, Sigma, P, Betas, X, pc) {
  
  X <- as.matrix(X)
  P <- P/sum(P)
  
  Pi <- exp(X%*%Betas)/(1+exp(X%*%Betas))
  
  W <- rep(0, n)
  Y <- rep(0, n)
  D <- rep(0, n)
  
  aux <- Par(Mu, Sigma)
  tau <- 1/Tau(pc, P, Mu, Sigma) 
  
  for(i in 1:n) {
    W[i] <- rbinom(1, 1, Pi[i])
    if(W[i]==1) {
      #tau.cr <- 1/Tau(Pi[i], P, Mu, Sigma) 
      Y[i] <- Inf #rexp(1, tau)
      D[i] <- 0
    } else {
      C <- rexp(1, tau)
      
      Z <- rmultinom(1, 1, P)
      T <- Z[1,]*rgamma(1, shape=aux$alpha, scale=aux$beta) + Z[2,]*rlnorm(1, aux$mu_l, aux$sig_l) + Z[3,]*rweibull(1, shape=aux$gama, scale=aux$delta)
      
      Y[i] <- min(T,C)
      D[i] <- ifelse(Y[i]==T, 1, 0)
    }
  }
  Xnames <- as.parm.names(list(X=rep(0,ncol(X))))
  out <- data.frame(cbind(Y, D, X)) ; colnames(out) <- c('Y', 'D', Xnames)
  return(out)
}


#-----------------------------------------------------------------------------------------------------------#
#-----Generates a sample of the GLW Promotion Time Cure Rate model with covariates -------------------------#
#-----------------------------------------------------------------------------------------------------------#
#' Random Generation From The GLW FM Promotion Time Cure Rate Model
#' 
#' \code{rpt_glw} generates a sample from the glw fm promotion time cure rate regression model 
#' with right censored times given a model matrix and the regression coefficients.
#' 
#' @param n: number of observations.
#' @param Mu: mean
#' @param Sigma: variance.
#' @param P: weights of the finite mixture.
#' @param Betas: regression coefficients associated with the long term parameter.
#' @param X: matrix of covariates (with column of 1's).
#' @param pc: vector with probability of right censored times.
#' 
#' @return A data frame with first column being the times, the second column the censorship 
#' indicator and the remaining columns are the covariates matrix \code{X}.
#' 
#' @export
rpt_glw <- function(n, Mu, Sigma, P, Betas, X, pc) {
  
  Y <- rep(0, n)
  D <- rep(0, n)
  N <- rep(0, n)
  X <- as.matrix(X)
  P <- P/sum(P)
  
  eta <- exp(X%*%Betas)
  tau <- 1/Tau(pc, P, Mu, Sigma)
  
  for(i in 1:n) {
    N[i] <- rpois(1, eta[i])
    
    if(N[i]==0) {
      #tau.pt <- 1/Tau(exp(-eta[i]), P, Mu, Sigma)
      Y[i] <- Inf #rexp(1, tau)
      D[i] <- 0
    } 
    else {
      R <- min(rglw(n=N[i], P=P, Mu=Mu, Sigma=Sigma))
      C <- rexp(1, tau)
      Y[i] <- min(R,C)
      D[i] <- ifelse(Y[i] == R, 1, 0)
    }
  }
  Xnames <- as.parm.names(list(X = rep(0, ncol(X))))
  out <- data.frame(cbind(Y, D, X))
  colnames(out) <- c("Y", "D", Xnames)
  return(out)
}
