#' Dean score test for overdispersion
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


dtest.ames<-function(data,count,dose,
                     alternative.var=c("gp1","nb2"),
                     alternative=c("greater","less","two.sided"),
                     predictor=c("krewski","myers","bernstein","breslow","stead","margolin","svetliza"),
                     link=log,
                     theta.start="smart",
                     controls=control.fitames(),
                     quiet=F
){

  if(missing(data)==T){
    data<-data.frame(y=count,x=dose)
  }else{
    if(class(data)=="sim.ames"){
      if(is.null(which.y)==T) which.y=1
      data<-data.frame(x=data[[1]],y=data[[(which.y+1)]] )
    }else if(class(data)=="fitames"){
      if(as.character( data$model$Likelihood )!='poisson') stop('Use a model with Poisson random component.')
      controls <- data$controls
      mu <- data$fitted.values
      y  <- data$data$count
      if(missing(alternative)==T){
        alternative <- "greater"
      }else alternative <- match.arg(alternative)
      alternative.var <- match.arg(alternative.var)

      switch(alternative.var,
             nb2 = {
               Vhat <- (1/2)*sum( mu^2 )
               Ti   <- (1/2)*((y - mu )^2 - y)
             },
             gp1 = {
               Vhat <- length(y)/2
               Ti   <- (1/2)*((y - mu )^2 - y)/mu
             }
      )
      Sstat <- sum(Ti)/sqrt(Vhat)
      switch(alternative,
             greater = {
               pvalue <- 1-pnorm(Sstat)
             },
             less = {
               pvalue <- pnorm(Sstat)
               warning('p-value probably doesn\'t make sense')
             },
             two.sided = {
               pvalue <- (1-pnorm(abs(Sstat)))*2
               warning('p-value probably doesn\'t make sense')
             }
      )

      Test <- data.frame(Var=alternative.var,Z.statistic=Sstat,pvalue=pvalue)
      if(quiet==F){
        print(Test)
      }
      return(invisible(list(Test=Test,Tis=Ti,Tis.Var=Vhat)))
    }else{
      data<-data.frame(y=data[,1],x=data[,2])
    }
  }

  if(missing(alternative)==T){
    alternative <- "greater"
  }else alternative <- match.arg(alternative)
  alternative.var <- match.arg(alternative.var)
  predictor       <- match.arg(predictor)

  lambda<-1
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
  if(length(which(data$y==0))>0) data$y[which(data$y==0)] <- controls$zerocounts

  if(theta.start[1]=="smart"){
    theta.start <- iv.ames(data,predictor=predictor,)
  }

  switch(predictor,
         krewski = {
           alpha <- controls$alpha
           Eta <- function(x,theta) (theta[1]+theta[2]*x)*exp(-theta[3]*(x^alpha))
           coefnam <- c("beta0","beta1","gamma")
         },
         myers = {
           Eta <- function(x,theta) (theta[1]+theta[2]*x)*exp(-theta[3]*x)
           coefnam <- c("beta0","beta1","gamma")
         },
         bernstein = {
           Eta <- function(x,theta) theta[1]+theta[2]*x
           coefnam <- c("beta0","beta1")
         },
         breslow = {
           delta <- controls$delta
           Eta <- function(x,theta) theta[1] + theta[2]*log((x+delta)) -theta[3]*x
           coefnam <- c("beta0","beta1","gamma")
         },
         stead = {
           if(length(which(data$x==0))>0) data$x[which(data$x==0)] <- 1e-3
           Eta <- function(x,theta) (theta[1]+theta[2]*(x^theta[3]))*exp(-theta[4]*x)
           coefnam <- c("beta0","beta1","beta2","gamma")
         },
         svetliza = {
           Eta <- function(x,theta) exp(theta[1]-exp(theta[2]-theta[3]*x))
           coefnam <- c("alpha","beta","gamma")
         },
         margolin = {
           m <- link(controls$cfu)
           Eta <- function(x,theta) m*(1-exp(-(theta[1]+theta[2]*x)))*exp(-theta[3]*x)
           coefnam <- c("beta0","beta1","gamma")
         }
  )

  invlink <- exp
  Mu  <- function(x,theta) invlink(Eta(x,theta))


  poissonmod <- suppressWarnings(
    fit.ames(count=data$y,dose=data$x,
             lh='poisson',predictor=predictor,theta.start=theta.start,
             controls=controls,quiet=T)
  )

  tries <- 0
  while( is.null(poissonmod$coef)==T && tries < 1000){
    poissonmod <- suppressWarnings(
      fit.ames(count=data$y,dose=data$x,
               lh='poisson',predictor=predictor,theta.start=jitter(theta.start),
               controls=controls,quiet=T)
    )
    tries <- tries + 1
  }

  if(is.null(poissonmod$coef)==F){
    theta.hat <- poissonmod$coef$Estimate
  }else return(NULL)

  switch(alternative.var,
         nb2 = {
           Vhat <- (1/2)*sum( Mu(data$x,theta.hat)^2 )
           Ti   <- (1/2)* ((data$y - Mu(data$x,theta.hat) )^2 - data$y)
         },
         gp1 = {
           Vhat <- length(data$y)/2
           Ti   <- (1/2)*((data$y - Mu(data$x,theta.hat) )^2 - data$y)/Mu(data$x,theta.hat)
         }
  )
  Sstat <- sum(Ti)/sqrt(Vhat)

  switch(alternative,
         greater = {
           pvalue <- 1-pnorm(Sstat)
         },
         less = {
           pvalue <- pnorm(Sstat)
           warning('p-value probably doesn\'t make sense')
         },
         two.sided = {
           pvalue <- (1-pnorm(abs(Sstat)))*2
           warning('p-value probably doesn\'t make sense')
         }
  )

  Test <- data.frame(Var=alternative.var,Z.statistic=Sstat,pvalue=pvalue)
  if(quiet==F){
    print(Test)
  }
  return(invisible(list(Test=Test,Tis=Ti,Tis.Var=Vhat)))
}
