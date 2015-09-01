#' Type-I generalized Poisson distribution
#'
#'
#' @param ...
#'
#' @author Davi Butturi-Gomes
#'
#' Silvio S. Zocchi
#'
#' @examples
#' fit.ames()
#'
#'
#' @export


rgenpois<-function(n,mu,phi){
  if(length(mu)==1){
    VY <- mu*phi
    minmas <- max( c(0, floor(mu - 3*VY )) )
    maxmas <- ceiling( mu + ifelse(mu<5,20*VY,4*VY) )
    probs <- dgenpois(minmas:maxmas,mu=mu,phi=phi)
    GP    <- DiscreteDistribution(supp=minmas:maxmas,prob=probs)
    y <- r(GP)(n)
  }
  if(length(mu)>1){
    VY <- mu*phi
    maxmas <- do.call(c,lapply(1:length(mu),function(i) ceiling( mu[i] + ifelse(mu[i]<5,20*VY[i],4*VY[i]) ) ))
    minmas <- do.call(c,lapply(1:length(mu),function(i) max(c(0, floor(mu[i] - 3*VY[i]))) ))

    supps <- do.call(list,lapply(1:length(mu),
                                 function(i) minmas[i]:maxmas[i]))

    probs <- do.call(list,lapply(1:length(mu),
                                 function(i) dgenpois(x=minmas[i]:maxmas[i],mu=mu[i],phi=phi)))

    GPs <- do.call(list,lapply(1:length(mu),function(i) DiscreteDistribution(supp=supps[[i]],prob=probs[[i]])))

    y <- do.call(c,lapply(1:length(mu),function(i) r(GPs[[i]])(1)))
    if(n>length(mu)) warning('n forced to have same length as mu.')
  }
  return(y)
}
