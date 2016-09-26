#---------------------------------------------------------------------------------------------------------------------#
##---Inversao mais eficiente para k covariáveis (Mesma media nas 3 distribuicoes)-------------------------------------#
#---------------------------------------------------------------------------------------------------------------------#
Par <- function(mu, sig2) {
 
 n <- length(mu)
 
 alpha <- mu^2/sig2               #Shape
 beta  <- sig2/mu         #Scale  
 
 mu_l  <- log(mu^2/sqrt(mu^2+sig2))               #Mean
 sig_l <- sqrt(log((mu^2+sig2)/mu^2))             #StdDev
 
 gama <- c()
 for(j in 1:n){
   f <- function(k) {log((sig2+mu[j]^2)/mu[j]^2) - lgamma(1+2/k) + 2*lgamma(1+1/k)}
   gama[j] <- uniroot(f,interval=c(10^-200, 10^10))$root      #Shape
 }
 delta <- mu/gamma(1+(1/gama))                                    #Scale
 
 aux <- as.data.frame(cbind(mu, alpha, beta, mu_l, sig_l, gama, delta))
 colnames(aux) <- c('mu', 'alpha', 'beta', 'mu_l', 'sig_l', 'gama', 'delta')
 return(aux)
}

Par3 <- function(Mu, S2) {
  
 Mu <- round(Mu, 10) ; tabela = table(Mu)
 mu.tab <- as.numeric(names(tabela))
  
 aux.par <- round(Par(mu.tab, S2),10)
 param <- matrix(nrow=nrow(as.matrix(Mu)), ncol=6)	
  
 for(i in 1:length(mu.tab)) param[Mu==aux.par[i,1],] <- matrix(as.numeric(aux.par[i,2:7]), ncol=6, nrow=as.numeric(tabela[i]), byrow=T)
 
 colnames(param) <- c('alpha', 'beta', 'mu_l', 'sig_l', 'gama', 'delta')
 return(as.data.frame(param))
}

#----------------------------------------------------------------------------------------------------------------------#
#--Inversao mais eficiente para k covariaveis  para o modelo com intercepto variavel-----------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
## Mu.gamma, Mu.lnorm, Mu.weibull: medias de cada componente da mistura;
## S2: variancia comum a todos os individuos e componentes;
Parm <- function(Mu.gamma, Mu.lnorm, Mu.weibull, S2){
  
  n <- length(Mu.gamma)
  
  alpha <- Mu.gamma^2/S2    #Shape
  beta  <- S2/Mu.gamma       #Scale  
  
  mu_l  <- log(Mu.lnorm^2/sqrt(Mu.lnorm^2+S2))       #Mean
  sig_l <- sqrt(log((Mu.lnorm^2+S2)/Mu.lnorm^2))     #StdDev
  
  gama <- c()
  for(j in 1:n){
    f <- function(k) {log((S2+Mu.weibull[j]^2)/Mu.weibull[j]^2) - lgamma(1+2/k) + 2*lgamma(1+1/k)}
    gama[j] <- uniroot(f,interval=c(10^-200, 10^10))$root      #Shape
  }
  delta <- Mu.weibull/gamma(1+(1/gama))    #Scale
  
  aux <- as.data.frame(cbind(Mu.gamma, alpha, beta, mu_l, sig_l, gama, delta))
  colnames(aux) <- c('Mu.gamma', 'alpha', 'beta', 'mu_l', 'sig_l', 'gama', 'delta')
  return(aux)
}
Parm3 <- function(Mu.gamma, Mu.lnorm, Mu.weibull, S2) {
  
 Mu.ga <- round(Mu.gamma, 10) ; tabela.ga <- table(Mu.ga)
 mu.tab.ga <- as.numeric(names(tabela.ga))
  
 Mu.ln <- round(Mu.lnorm, 10) ; tabela <- table(Mu.ln)
 mu.tab.ln <- as.numeric(names(tabela))
  
 Mu.we <- round(Mu.weibull, 10) ; tabela <- table(Mu.we)
 mu.tab.we <- as.numeric(names(tabela))
  
 aux.par <- round(Parm(mu.tab.ga, mu.tab.ln, mu.tab.we, S2),10)
 param <- matrix(nrow=nrow(as.matrix(Mu.gamma)), ncol=6)
  
 for(i in 1:length(mu.tab.ga)) param[Mu.ga==aux.par[i,1],] <- matrix(as.numeric(aux.par[i,2:7]), ncol=6, nrow=as.numeric(tabela.ga[i]), byrow=T)
  
 colnames(param) <- c('alpha', 'beta', 'mu_l', 'sig_l', 'gama', 'delta')
 return(as.data.frame(param))
 }

#----------------------------------------------------------#
#------Funcoes do modelo GLW-------------------------------#
#----------------------------------------------------------#
#x: quantil
#n: tamanho da amostra
#P: vetor de pesos da mistura (normalizado para somar 1)
#Mu: Media
#S2: variancia

#' Density of GLW Finite Mixture
#' 
#' \code{dglw} returns the value of the density of the GLW finite mixtue.
#' @param x quantile
#' @param P weights
#' @param Mu mean
#' @param S2 variance
#' 
#' @return Density of the Gamma Lognormal Weibull finite mixture with mean Mu, variance S2, and vector of weights P.
#' 
#' @family pglw rglw
#' 
#' @export
dglw <- function(x, P, Mu, Sigma) {
 P <- P/sum(P)
 a <- Par3(Mu, Sigma)
 d <- P[1]*dgamma(x, shape=a$alpha, scale=a$beta) + P[2]*dlnorm(x, a$mu_l, a$sig_l) + P[3]*dweibull(x, shape=a$gama, scale=a$delta)
 return(d)
 }

#' Distribution Function of GLW Finite Mixture
#' 
#' \code{pglw} returns the cdf of the GLW finite mixtue.
#' @param x vector of quantiles.
#' @param P vector of mixture weights.
#' @param Mu mean.
#' @param Sigma variance.
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' @return CDF of the Gamma Lognormal Weibull finite mixture with mean Mu, variance Sigma, and vector of weights P.
#' @family dglw rglw
#' 
#' @export
pglw <- function(x, P, Mu, Sigma, lower.tail=T) {
 P <- P/sum(P)
 a <- Par3(Mu, Sigma)
 d <- P[1]*pgamma(x, shape=a$alpha, scale=a$beta, lower.tail=lower.tail) + P[2]*plnorm(x, a$mu_l, a$sig_l, lower.tail=lower.tail) + P[3]*pweibull(x, shape=a$gama, scale=a$delta, lower.tail=lower.tail)
 return(d)
 }

#' Random Generation From the GLW Finite Mixture
#' 
#' \code{rglw} generates a sample from the GLW finite mixtue.
#' 
#' @param n number of observations.
#' @param P vector of mixture weights.
#' @param Mu mean.
#' @param Sigma variance.
#' 
#' @return Generates a sample with \code{n} observations from the Gamma Lognormal Weibull finite mixture model with mean Mu, variance Sigma, and vector of weights P.
#' @seealso \code{\link{pglw}} for cdf, \code{\link{dglw}} for density
#' 
#' @export
rglw <- function(n, P, Mu, Sigma) {
 a <- Par3(Mu, Sigma)
 z <- rmultinom(n, 1, P)
 d <- z[1,]*rgamma(n, shape=a$alpha, scale=a$beta) + z[2,]*rlnorm(n, a$mu_l, a$sig_l) + z[3,]*rweibull(n, shape=a$gama, scale=a$delta)
 return(d)
 }

