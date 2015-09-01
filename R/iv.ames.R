#' Compute inital values for fitting models to Ames test data
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

iv.ames <- function(data,count,dose,
                    predictor=c("krewski","myers","margolin","stead"),
                    link=log,controls=control.fitames(),
                    interval=c(-5,5),submod=F,...){

  predictor <- match.arg(predictor)
  #
  b2.stead <- seq(1,5,length.out=1000)
  which.x  <- 'simple'
  ufcs     <- controls$cfu
  alpha    <- controls$alpha
  delta    <- controls$delta
  #
  if(missing(data)==T){
    data<-data.frame(y=count,x=dose)
  }
  y <- data[,1]
  x <- data[,2]

  lambda <- 1
  if(lambda!=1&lambda!=0){
    x1 <- x^lambda
  }
  else if(lambda==0){
    x1 <- log(x)
  }
  else if(lambda==1){
    x1<-x
  }

  if(length(which(x1==-Inf | x1==Inf ) )!=0 ){
    x1 <- log(x+min(x[which(x>0)]))
  }
  x<-x1;rm(x1)

  if(which.x=='simple'){
    xfac <- as.factor(x)
    rep  <- unname(table(xfac))
    x.0  <- c(x[1+rep[1]],x[1+rep[1]+rep[2]]  )
    y.0  <- c(link(mean(y[1:rep[1]] )),link(mean(y[(rep[1]+1):(rep[1]+rep[2])])),link(mean(y[(rep[1]+rep[2]+1):(sum(rep[1:3]))] )))
  }

  if(which.x!='simple') stop('This method is not implemented')

  if(predictor=='myers'){
    if(submod){
      i.maxmut <- 1:max( which( x == x[ which( y == max(y) ) ] ) )
      i.minmut <- min( which( x == x[ which( y == max(y) ) ] ) ):length(x)

      x.maxmut <- x[i.maxmut]
      x.minmut <- x[i.minmut]
      y.maxmut <- y[i.maxmut]
      y.minmut <- y[i.minmut]

      if(length(which(x==0))>0){
        b0.0 <- link( mean(y[which(x==0)]) )
      }else{
        maxmod <- glm(y.maxmut ~ x.maxmut, family=poisson )
        b0.0 <- unname( maxmod$coef[1] )
      }

      minmod <- glm(y.minmut ~ x.minmut, family=poisson )
      g.0    <- -unname( minmod$coef[2] )

      offb0.0 <- rep(b0.0,length(y.maxmut))
      offg.0 <- -exp(g.0*x.maxmut)

      medmod <- glm(y.maxmut ~ -1 + offset(offb0.0) + offset(offg.0) + x.maxmut, family=poisson )

      b1.0.1 <- unname( medmod$coef[1] )
      b1.0.2 <- unname( glm(y.maxmut ~ -1 + x.maxmut, family=poisson )$coef[1] - glm(y.maxmut ~  x.maxmut, family=poisson )$coef[2] )
      b1.0  <- max(b1.0.1,b1.0.2)

      st <- data.frame(b0.0=b0.0, b1.0=b1.0, g.0=g.0 )
    }else{
      z.myers <- function(z,y0,x0){
        y0[1]*x0[1]*exp(z+z*x0[1]/x0[2])-x0[2]*y0[1]*exp(z+z*x0[1]/x0[2])+exp(z)*x0[2]*y0[2]-y0[3]*exp(z*x0[1]/x0[2])*x0[1]
      }

      b1.0fun <- function(zroots,y0,x0){
        obj <- (-exp(zroots*x0[1]/x0[2])*y0[1]+y0[2])/(exp(zroots*x0[1]/x0[2])*x0[1])
        if(length(obj[which(obj>0)])==0)obj2<-0
        else obj2 <- obj[which(obj>0)]
        return(obj2)
      }
      g.0fun <- function(zroots,x0){
        obj <- zroots/x0[2]
        if(length(obj[which(obj<0)])==0)obj2<-0
        else obj2 <- obj[which(obj<0) ]
        return(obj2)
      }
      zs <- uniroot.all(z.myers,interval=interval,y0=y.0,x0=x.0)
      st <- data.frame( b0.0=y.0[1] , b1.0=b1.0fun(zs,y.0,x.0) , g.0=-g.0fun(zs,x.0) )
    }
  }

  else if(predictor=='margolin'){
    if(submod){
      i.maxmut <- 1:max( which( x == x[ which( y == max(y) ) ] ) )
      i.minmut <- min( which( x == x[ which( y == max(y) ) ] ) ):length(x)

      x.maxmut <- x[i.maxmut]
      x.minmut <- x[i.minmut]
      y.maxmut <- y[i.maxmut]
      y.minmut <- y[i.minmut]

      if(length(which(x==0))>0){
        b0.0 <- -log( 1 - link( mean(y[which(x==0)]))/link(ufcs) )
        offb0.0 <- rep(link( mean(y[which(x==0)])) ,length(y.maxmut) )
      }else{
        maxmod  <- glm(y.maxmut ~ x.maxmut, family=poisson )
        b0.0    <- unname( -log( 1 - maxmod$coef[1]/link(ufcs) ) )
        offb0.0 <- unname( rep( maxmod$coef[1] ,length(y.maxmut) ) )

      }
      b1.00 <- unname( glm(y.maxmut ~ -1 + offset(offb0.0) + x.maxmut, family=poisson)$coef )
      b1.0  <- -( log( 1- (b1.00/link(ufcs)) ) - b0.0 )
      offb0b1 <- -link(ufcs)*(1-exp(-b0.0 - b1.0*x.minmut))

      minmod  <- glm(y.minmut ~  x.minmut + offset(offb0b1) , family=poisson)
      g.0 <- unname( minmod$coef[2]/2 )

      st <- data.frame(b0.0=b0.0, b1.0=b1.0, g.0=g.0 )
    }else{
      z.margo <- function(z,nn,y0,x0){
        nn*exp(z+ (log( (nn-y0[1])/(nn*exp(z*x0[1]/x0[2])-y0[2]) )*x0[2]+z*x0[1]) /x0[1])-
          exp(z)*nn+
          exp(z)*y0[1]-
          y0[3]*exp( (log((nn-y0[1])/(nn*exp(z*x0[1]/x0[2])-y0[2]))*x0[2]+z*x0[1])/x0[1] )
      }
      b0.0fun <- function(y0,n=link(ufcs)){
        obj  <- -log((n-y.0[1])/n)
        obj2 <- obj[which(obj>0)]
        return(obj2)
      }
      b1.0fun <- function(zroots,y0,x0,nn=link(ufcs)){
        obj <- (1/(x.0[1]*x.0[2]))*(log((nn-y0[1])/(nn*exp(zroots*x0[1]/x0[2])-y0[2]))*x0[2]+zroots*x0[1])
        obj2 <- obj[which(obj>0)]
        return(obj2)
      }
      g.0fun <- function(zroots,x0){
        obj <- zroots/x0[2]
        obj2 <- obj[which(obj<0) ]
        return(obj2)
      }
      zs <- suppressWarnings( uniroot.all(z.margo,interval=interval,y0=y.0,x0=x.0,nn=link(ufcs)) )
      st <- data.frame(b0.0 = b0.0fun(y.0), b1.0 = b1.0fun(zs,y.0,x.0), g.0 = -g.0fun(zs,x.0) )
    }
  }

  else if(predictor=='krewski'){
    if(submod){
      alpha <- controls$alpha

      i.maxmut <- 1:max( which( x == x[ which( y == max(y) ) ] ) )
      i.minmut <- min( which( x == x[ which( y == max(y) ) ] ) ):length(x)

      x.maxmut <- x[i.maxmut]
      x.minmut <- x[i.minmut]
      y.maxmut <- y[i.maxmut]
      y.minmut <- y[i.minmut]

      if(length(which(x==0))>0){
        b0.0 <- link( mean(y[which(x==0)]) )
      }else{
        maxmod <- glm(y.maxmut ~ x.maxmut, family=poisson )
        b0.0   <- unname( maxmod$coef[1] )
      }

      minmod <- glm(y.minmut ~ I(alpha*x.minmut), family=poisson)
      g.0  <- -unname( minmod$coef[2]/10 )

      offb0.0 <- rep(b0.0,length(y.maxmut))
      offg.0 <- -exp(g.0*(x.maxmut^alpha))

      medmod1 <- glm(y.maxmut ~ -1 + offset(offb0.0) + offset(offg.0) + x.maxmut,family=poisson)
      medmod2 <- glm(y.maxmut ~ offset(offg.0) + x.maxmut,family=poisson )

      b1.0 <- unname(2*medmod1$coef[1] - medmod2$coef[2])

      st <- data.frame(b0.0=b0.0, b1.0=b1.0, g.0=g.0 )
    }else{
      z.krewski <- function(z,y0,x0,alpha){
        if(missing(alpha)){
          alpha <- 2
        }
        -y0[1]*x0[2]*exp(z*(1+x0[1]^(-alpha)*x0[2]^alpha))+x0[1]*y0[1]*exp(z*(1+x0[1]^(-alpha)*x0[2]^alpha))-exp(z)*x0[1]*y0[3]+y0[2]*exp(z*x0[1]^(-alpha)*x0[2]^alpha)*x0[2]
      }
      b1.0fun <- function(zroots,y0,x0,alpha){
        if(missing(alpha)){
          alpha <- 2
        }
        obj <- (-y0[1]+exp(-zroots*x0[2]^alpha*x0[1]^(-alpha))*y0[3])/x0[2]
        if(length(obj[which(obj>0)])==0)obj2<-0
        else obj2 <- obj[which(obj>0)]
        return(obj2)
      }
      g.0fun <- function(zroots,x0,alpha){
        if(missing(alpha)){
          alpha <- 2
        }
        obj <- zroots*x0[1]^(-alpha)
        if(length(obj[which(obj<0)])==0)obj2<-0
        else obj2 <- obj[which(obj<0) ]
        return(obj2)
      }
      zs <- uniroot.all(z.krewski,interval=interval,y0=y.0,x0=x.0,alpha=controls$alpha)
      st <- data.frame( b0.0=y.0[1] , b1.0=b1.0fun(zs,y.0,x.0,alpha=controls$alpha) , g.0=-g.0fun(zs,x.0,alpha=controls$alpha) )
    }
  }
  else if(predictor=='stead'){
    if(submod){
      i.maxmut <- 1:max( which( x == x[ which( y == max(y) ) ] ) )
      i.minmut <- min( which( x == x[ which( y == max(y) ) ] ) ):length(x)

      x.maxmut <- x[i.maxmut]
      x.minmut <- x[i.minmut]
      y.maxmut <- y[i.maxmut]
      y.minmut <- y[i.minmut]

      g.0     <- -unname( glm(y.minmut~x.minmut,family=poisson)$coef[2])
      b0.0off <- rep(b0.0<-log( mean(y[1:3]) ), length(y.maxmut))

      b2.00 <- seq(0,5,length.out=1000)[-1]
      lb2.p <- c()
      for(i in 1:length(b2.00)){
        lb2.p <- c(lb2.p, as.numeric(logLik( glm(y.maxmut~-1+offset(b0.0off)+I(x.maxmut^b2.00[i]),family=poisson))))
      }
      b2.0 <- b2.00[ which(lb2.p==max(lb2.p)) ]
      b1.0 <- unname( log( glm(y.maxmut~-1+offset(b0.0off) + I(x.maxmut^b2.0),family=poisson)$coef^-b2.0 ) )

      st <- data.frame(b0.0=b0.0, b1.0=b1.0, b2.0=b2.0, g.0=g.0 )
    }else{
      b2.00   <- b2.stead
      z.stead <- function(z,y0,x0,b2.00){
        -y0[1]*(x0[1]^b2.00)*exp(z*(x0[2]+x0[1])/x0[2])+(x0[2]^b2.00)*y0[1]*exp(z*(x0[2]+x0[1])/x0[2])-exp(z)*(x0[2]^b2.00)*y0[2]+y0[3]*exp(z*x0[1]/x0[2])*(x0[1]^b2.00)
      }
      b1.0fun <- function(zroots,y0,x0,b2.00) (1/(x0[1]^b2.00))*(y0[2]*exp(-zroots) - y0[1])
      g.0fun  <- function(zroots,x0) zroots/x0[2]

      b1.0<-NULL;g.0<-NULL;b2.0<-0;kb1<-matrix(0,ncol=1,nrow=1);kg<-matrix(0,ncol=1,nrow=1)
      for(i in 1:length(b2.00)){
        zs   <- uniroot.all(z.stead,interval=interval,y0=y.0,x0=x.0,b2.00=b2.00[i])
        b1.0[[i]] <- matrix(b1.0fun(zs,y.0,x.0,b2.00[i]),nrow=length(zs),ncol=1)
        g.0[[i]]  <- matrix(g.0fun(zs,x.0),nrow=length(zs),ncol=1)
        b2.0 <- c(b2.0,rep(b2.00[i],length(zs)))
      }

      for(i in 1:length(b1.0)){
        kb1 <- rbind(kb1,b1.0[[i]])
        kg  <- rbind(kg,g.0[[i]])
      }

      b0.0<-rep(y.0[1],length(b2.0))

      st.0 <- data.frame(b0.0=b0.0[-1],b1.0=kb1[-1,],b2.0=b2.0[-1],g.0=kg[-1,] )
      st.1 <- st.0[which(st.0$b1.0>0&st.0$b1.0<10&st.0$g.0<0&st.0$g.0>-10),]
      st   <- data.frame(b0.0=st.1$b0.0[1],
                         b1.0=min(st.1[st.1$b2.0==min(st.1$b2.0),]$b1.0),
                         b2.0=st.1[st.1$b2.0==min(st.1$b2.0),]$b2.0[1],
                         g.0=-min(st.1[st.1$b2.0==min(st.1$b2.0),]$g.0))
    }
  }
  return(st)
}
