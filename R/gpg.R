#' Generalized Poisson Type-G distribution (GP-G)
#'
#' \code{dgpg} - Generalized Poisson type-G distribution probability mass function
#'
#' @param x Integer, vector of integers or numeric that can be coerced into integers.
#' @param n Integer or numeric that can be coerced into integer. Number of random deviations
#' to be generated
#' @param mu Positive numeric. Mean parameter of the distribution
#' @param phi Positive numeric. Dispersion parameter.
#' @param g Positive numeric <= 1. Functional parameter of the distribution.
#' @param log Logical. If TRUE, return the probabilities in log scale. The default is FALSE.
#' @param use.gp12 Logical. If TRUE (default) recovers the GP-1 parametrization for g=1 or,
#' for g=2, recovers the GP-2 parametrization. This is recommended due otimization purposes.
#'
#' @author Davi Butturi-Gomes
#'
#' Silvio S. Zocchi
#'
#' @seealso \code{\link{dgp1}}, \code{\link{rgp1}}, \code{\link{dgp2}}, \code{\link{rgp2}}
#'
#' @examples
#' x <- rgpg(1000,mu=mu<-20,phi=phi<-0.01,g=g<-2) # Same as using rgp2.
#' mu; mu*((1+phi*(mu^(g-1)))^2) # nominal mean and variance
#' mean(x);var(x) # "observed" mean and variance
#' #
#' dx <- (dgpg(10:30,mu=mu,phi=phi,g=g))
#' plot(10:30,dx,type='h')
#'
#' @rdname rgpg
#' @export

dgpg <- function(x,mu,phi,g,log=F,use.gp12=T){
  if(is.integer(x)==F){
    warning('Non-integer x')
    px<-0
  }
  else{
    if(g==1 && use.gp12) phi <- sqrt(phi)-1
    px <- log(mu) + (x-1)*(log(mu+phi*x*mu^(g-1))) - x*log(1+phi*mu^(g-1)) -
      lfactorial(x) - (mu+phi*x*mu^(g-1))/(1+phi*mu^(g-1))
    if(!log) px <- exp(px)
  }
  return(px)
}

#' Generalized Poisson Type-G distribution (GP-G)
#'
#' \code{rgpg} - Generalized Poisson type-G distribution random number generator
#'
#' @rdname rgpg
#' @export

rgpg<-function(n,mu,phi,g,use.gp12=T){
  if(g==1 && use.gp12){
    y <- rgp1(n,mu,phi)
  }else if(g==2 && use.gp12){
    y <- rgp2(n,mu,phi)
  }else{
    if(length(mu)==1){
      VY <- mu*((1+phi*(mu^(g-1)))^2)
      CV <- VY/mu
      minmas <- max( c(0, floor(mu - ifelse(CV < 5 , 3*VY, VY )  )) )
      maxmas <- ifelse(CV < 5, ceiling( mu + ifelse(mu<5,10*VY,3*VY) ), ceiling( mu + ifelse(mu<5,3*VY,VY)) )
      probs <- dgpg(minmas:maxmas,mu=mu,phi=phi,g=g)
      GP    <- DiscreteDistribution(supp=minmas:maxmas,prob=probs)
      y <- r(GP)(n)
    }
    if(length(mu)>1){
      VY <- mu*((1+phi*(mu^(g-1)))^2)
      CV <- VY/mu
      maxmas <- do.call(c,lapply(1:length(mu),
                                 function(i) ifelse(CV[i] < 5, ceiling( mu[i] + ifelse(mu[i]<5,10*VY[i],3*VY[i]) ),
                                                    ceiling( mu[i] + ifelse(mu[i]<5,3*VY[i],VY[i]))) ))
      minmas <- do.call(c,lapply(1:length(mu),
                                 function(i) max( c(0, floor(mu[i] - ifelse(CV[i] < 5 , 3*VY[i], VY[i] )))) ))

      supps <- do.call(list,lapply(1:length(mu),
                                   function(i) minmas[i]:maxmas[i]))

      probs <- do.call(list,lapply(1:length(mu),
                                   function(i) dgpg(x=minmas[i]:maxmas[i],mu=mu[i],phi=phi[i],g=g[i])))

      GPs <- do.call(list,lapply(1:length(mu),function(i) DiscreteDistribution(supp=supps[[i]],prob=probs[[i]])))

      y <- do.call(c,lapply(1:length(mu),function(i) r(GPs[[i]])(1)))
      if(n>length(mu)) warning('n forced to have same length as mu.')
    }
  }
  return(y)
}
