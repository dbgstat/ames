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
#' @export


dgenpois <- function(x,mu,phi){
  if(is.integer(x)==F){
    warning('Non-integer x')
    px<-0
  }
  else{
    lk <- max(c(-1,-mu/4))
    k  <- 1-sqrt(1/phi)
    px <- (x-1)*log((1-k)*mu + k*x) + log( (1-k)*mu ) -
      lfactorial(x) -(1-k)*mu - k*x
    px<-exp(px)
  }
  return(px)
}
