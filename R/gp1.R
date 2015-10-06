#' Generalized Poisson Type-1 distribution (GP-1)
#'
#' \code{dgp1} - Generalized Poisson type-1 distribution probability mass function
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
#' @seealso \code{\link{dgp2}}, \code{\link{rgp2}}, \code{\link{dgpg}}, \code{\link{rgpg}}
#'
#' @examples
#' x <- rgp1(1000,mu=mu<-20,phi=phi<-2)
#' mu; mu*phi # Nominal mean and variance
#' mean(x);var(x) # "Observed" mean and variance
#' #
#' dx <- (dgp1(15:25,mu=mu,phi=phi))
#' plot(15:25,dx,type='h')
#'
#' @rdname rgp1
#' @export

dgp1 <- function(x,mu,phi,log=F){
  if(is.integer(x)==F){
    warning('Non-integer x')
    px<-0
  }
  else{
    px <- (x-1)*log((phi^(-1/2))*(mu-x) + x) + log(mu) - log(phi)/2 -
      ((phi^(-1/2))*(mu-x) + x + lfactorial(x))
    if(!log) px <- exp(px)
  }
  return(px)
}

#' Generalized Poisson Type-1 distribution (GP-1)
#'
#' \code{rgp1} - Generalized Poisson type-1 distribution random number generator
#'
#' @rdname rgp1
#' @export

rgp1 <-
  function(n,mu,phi){
    if(length(mu)==1){
      VY <- mu*phi
      CV <- VY/mu
      minmas <- max( c(0, floor(mu - ifelse(CV < 10 , 3*VY, VY ))))
      maxmas <- ifelse(CV < 10, ceiling( mu + ifelse(mu<5,10*VY,3*VY) ),
                       ceiling( mu + ifelse(mu<5, 3*VY,VY)))
      probs <- dgp1(minmas:maxmas,mu=mu,phi=phi)
      GP    <- DiscreteDistribution(supp=minmas:maxmas,prob=probs)
      y <- r(GP)(n)
    }
    if(length(mu)>1){
      VY <- mu*phi
      CV <- VY/mu
      maxmas <- do.call(c,lapply(1:length(mu),
                                 function(i) ifelse(CV[i] < 10, ceiling( mu[i] + ifelse(mu[i]<5,10*VY[i],3*VY[i]) ),
                                                    ceiling( mu[i] + ifelse(mu[i]<5, 3*VY[i],VY[i]))) ))
      minmas <- do.call(c,lapply(1:length(mu),
                                 function(i) max( c(0, floor(mu[i] - ifelse(CV[i] < 10 , 3*VY[i], VY[i] )))) ))

      supps <- do.call(list,lapply(1:length(mu),
                                   function(i) minmas[i]:maxmas[i]))

      probs <- do.call(list,lapply(1:length(mu),
                                   function(i) dgp1(x=minmas[i]:maxmas[i],mu=mu[i],phi=phi[i])))

      GPs <- do.call(list,lapply(1:length(mu),function(i) DiscreteDistribution(supp=supps[[i]],prob=probs[[i]])))

      y <- do.call(c,lapply(1:length(mu),function(i) r(GPs[[i]])(1)))
      if(n>length(mu)) warning('n forced to have same length as mu.')
    }
    return(y)
  }
