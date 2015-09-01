#' Ames test data random generator
#'
#'
#' @param ...
#'
#'
#' @author Davi Butturi-Gomes
#'
#' Silvio S. Zocchi
#'
#'
#'
#' @examples
#' iv.ames()
#'
#'
#' @export

rames <- function(dose,theta,dispersion=1,mu.only=F,
                  lh=c("poisson","nb2","gp1"),
                  predictor=c("krewski","myers","bernstein","breslow","stead","margolin","svetliza"),
                  link=log,
                  controls=control.fitames(),
                  quiet=T){

  data      <- data.frame(x=dose)
  lh        <- match.arg(lh)
  predictor <- match.arg(predictor)
  theta     <- abs(theta)

  lambda <- 1
  if(lambda!=1&lambda!=0){
    data$x1 <- data$x^lambda
  }else if(lambda==0){
    data$x1 <- log(data$x)
  }else if(lambda==1){
    data$x1<-data$x
  }
  if(length(which(data$x==-Inf | data$x==Inf ) )!=0 ){
    data$x1 <- log(data$x+min(data$x[which(data$x>0)]))
  }
  data$x <- data$x1

  switch(predictor,
         krewski = {
           alpha <- controls$alpha
           Eta <- function(x,theta) (theta[1]+theta[2]*x)*exp(-theta[3]*(x^alpha))
         },
         myers = {
           Eta <- function(x,theta) (theta[1]+theta[2]*x)*exp(-theta[3]*x)
         },
         bernstein = {
           Eta <- function(x,theta) theta[1]+theta[2]*x
         },
         breslow = {
           delta <- controls$delta
           Eta <- function(x,theta) theta[1] + theta[2]*log((x+delta)) -theta[3]*x
         },
         stead = {
           Eta <- function(x,theta) (theta[1]+theta[2]*(x^theta[3]))*exp(-theta[4]*x)
         },
         svetliza = {
           Eta <- function(x,theta) exp(theta[1]-exp(theta[2]-theta[3]*x))
         },
         margolin = {
           m <- link(controls$cfu)
           Eta <- function(x,theta) m*(1-exp(-(theta[1]+theta[2]*x)))*exp(-theta[3]*x)
         }
  )

  mu <- exp(Eta(data$x,theta))
  if(mu.only==T){
    return(mu)
  }

  switch(lh,
         poisson={
           rsample <- function(mu,phi) rpois(length(mu),lambda=mu)
         },
         nb2={
           rsample <- function(mu,phi) rnbinom(length(mu),mu=mu,size=phi)
         },
         gp1={
           rsample <- function(mu,phi) rgenpois(1,mu=mu,phi=phi)
         }
  )

  data$y <- rsample(mu=mu,phi=dispersion)
  data <- data.frame(count=data$y,dose=data$x)
  rets <- list(data=data,theta=theta,dispersion=dispersion,lh=lh,predictor=predictor,link=link,controls=controls)

  class(rets) <- "sim.ames"
  if(quiet==F) print(rets)
  return(rets)

}

