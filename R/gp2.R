#' Generalized Poisson Type-2 distribution (GP-2)
#'
#' \code{dgp2} - Generalized Poisson type-2 distribution probability mass function
#'
#' @param x Integer, vector of integers or numeric that can be coerced into integers.
#' @param n Integer or numeric that can be coerced into integer. Number of random deviations
#' to be generated
#' @param mu Positive numeric. Mean parameter of the distribution
#' @param phi Positive numeric. Dispersion parameter.
#' @param log Logical. If TRUE, return the probabilities in log scale. The default is FALSE.
#'
#' @author Davi Butturi-Gomes
#'
#' Silvio S. Zocchi
#'
#' @seealso \code{\link{dgp1}}, \code{\link{rgp1}}, \code{\link{dgpg}}, \code{\link{rgpg}}
#'
#' @examples
#' x <- rgp2(1000,mu=mu<-20,phi=phi<-0.01)
#' mu; mu*((1+phi*mu)^2) # nominal mean and variance
#' mean(x);var(x) # "observed" mean and variance
#' #
#' dx <- (dgp2(10:30,mu=mu,phi=phi))
#' plot(10:30,dx,type='h')
#'
#' @rdname rgp2
#' @export

dgp2 <- function(x,mu,phi,log=F){
  if(is.integer(x)==F){
    warning('Non-integer x')
    px<-0
  }
  else{
    px <- x*log(mu/(1+phi*mu)) + (x-1)*log(1+phi*x) - lfactorial(x) - mu*(1+phi*x)/(1+phi*mu)
    if(!log) px <- exp(px)
  }
  return(px)
}

#' Generalized Poisson Type-2 distribution (GP-2)
#'
#' \code{rgp2} - Generalized Poisson type-2 distribution random number generator
#'
#' @rdname rgp2
#' @export

rgp2<-function(n,mu,phi){
  if(length(mu)==1){
    VY <- mu*((1+phi*mu)^2)
    CV <- VY/mu
    minmas <- max( c(0, floor(mu - ifelse(CV < 5 , 3*VY, VY )  )) )
    maxmas <- ifelse(CV < 5, ceiling( mu + ifelse(mu<5,10*VY,3*VY) ), ceiling( mu + ifelse(mu<5,3*VY,VY)) )
    probs <- dgp2(minmas:maxmas,mu=mu,phi=phi)
    GP    <- DiscreteDistribution(supp=minmas:maxmas,prob=probs)
    y <- r(GP)(n)
  }
  if(length(mu)>1){
    VY <- mu*((1+phi*mu)^2)
    CV <- VY/mu
    maxmas <- do.call(c,lapply(1:length(mu),
                               function(i) ifelse(CV[i] < 5, ceiling( mu[i] + ifelse(mu[i]<5,10*VY[i],3*VY[i]) ),
                                                  ceiling( mu[i] + ifelse(mu[i]<5,3*VY[i],VY[i]))) ))
    minmas <- do.call(c,lapply(1:length(mu),
                               function(i) max( c(0, floor(mu[i] - ifelse(CV[i] < 5 , 3*VY[i], VY[i] )))) ))

    supps <- do.call(list,lapply(1:length(mu),
                                 function(i) minmas[i]:maxmas[i]))

    probs <- do.call(list,lapply(1:length(mu),
                                 function(i) dgp2(x=minmas[i]:maxmas[i],mu=mu[i],phi=phi[i])))

    GPs <- do.call(list,lapply(1:length(mu),function(i) DiscreteDistribution(supp=supps[[i]],prob=probs[[i]])))

    y <- do.call(c,lapply(1:length(mu),function(i) r(GPs[[i]])(1)))
    if(n>length(mu)) warning('n forced to have same length as mu.')
  }
  return(y)
}
