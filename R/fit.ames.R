#' Fitting models for Ames test data
#'
#' This is the main function of the "ames" package -- potentially, fits several different models
#' to the data. Its output is an object of class "fitames", used by most of the other functions
#' available in the package.
#'
#' @param data a matrix or a data frame. The standard use of this structure inside the function is
#' revernt count in the first column and doses placed in the second column. If not supplied,
#' the function will searc for the arguments \code{count} and \code{dose}.
#'
#' @param count Integer. A vector of length equal to the number of observations, describing the
#' number of revertent colonies per plate. If both \code{data} and \code{count} are supplied, then
#' \code{count} is ignored.
#'
#' @param dose Numeric. A vector of same length as \code{count}, describind the doses in each plate.
#'
#' @param lh Character. A sequence of characters describing which random component should used
#' for model fitting. For GNLMs, the valid sequences are "poisson", "nb2", "quasipoisson" and
#' "quasipower". For VGNLMs, the valid sequences are "skellam" and "nlvnb2". See \bold{Details}.
#'
#'
#' @author Davi Butturi-Gomes
#'
#' Silvio S. Zocchi
#'
#' @details The families described by the argument \code{lh} differ either by the complete
#' especification of likelihood function, or by the variance function formula. If "poisson"
#' is chosen, then a complete Poisson likelihood is used, whth \emph{phi} fixed (=1) and
#' variance function \emph{V(mu)=mu}. If "nb2", the assumptions depend on the argument
#' \code{dispersion.method}. Thus, for "ml", then a complete negative binomial type-II
#' likelihood is used, \emph{phi} is fixed (=1), and \emph{V(mu)=mu(1+mu/k)}; for "mpl", instead,
#' \emph{phi=mu(1+mu/k)} and \emph{V(mu)=1} and the estimation of \emph{k} is obtained by
#' maximum pseudo-likelihood. If "quasipoisson", no likelihood function is completely especified,
#' and the estimation procceeds using the maximum quasi-likelihood with \emph{V(mu)=mu},
#' \emph{Var(Y)=phi*V(mu)} and \emph{phi} an unknown constant, estimated by a moments method
#' (generalized Pearson statistic). If "quasipower".
#'
#'
#' @examples
#' fit.ames()
#'
#'
#' @export

fit.ames<- function(data,count,dose,
                    lh=c("poisson","nb2","quasipoisson","quasipower","skellam","nlvnb2"),
                    predictor=c("krewski","myers","bernstein","breslow","stead","margolin"),
                    link=log,
                    form1,form2,
                    theta.start="smart",
                    dispersion.start=1,
                    dispersion.method=c("ml","mpl","profile"),
                    controls=control.fitames(),
                    quiet=F,...){

  if(missing(data)==T){
    data<-data.frame(y=count,x=dose)
  }else{
    if(class(data)=="sim.ames"){
      lh        <- ifelse(data$lh!="gp1",data$lh,"quasipoisson")
      predictor <- data$predictor
      controls  <- data$controls
      if(missing(theta.start)) theta.start <- data$theta
      if(missing(dispersion.start)) dispersion.start <- data$dispersion
      data      <- data.frame(y=data$data$count,x=data$data$dose)
    }else{
      dataF  <- data
      data$y <- data[,1]
      data$x <- data[,2]
    }
  }
  if(length(dispersion.start)==1 & lh=='quasipower') stop('There is more than one dispersion parameter for Quasi-Power family')
  if(length(dispersion.start)==1 & lh=='hd') stop('There is more than one dispersion parameter for HD family')
  if(lh=='poisson'|lh=='nb2'|lh=='quasipoisson') dispersion.start <- dispersion.start[1]
  lambda<-1

  lh <- match.arg(lh)
  if(lh!="poisson" & lh!="nb2" & lh!="quasipoisson" & lh!="quasipower" & lh!="hd" & lh!="nlvnb2"){
    predictor <- NULL
    if(missing(form1)==T | missing(form2)==T){
      stop('Formulas are required for \'skellam\', \'pskellam\', \'vnb2\' and \'nlvnb2\' ')
    }
  }else if(missing(predictor)){
    stop('Choosing a predictor is required for non-linear models')
  }else{
    predictor <- match.arg(predictor)
  }

  if(lh=="nlvnb2"){
    if(missing(form1)==T&missing(form2)==T){
      stop("At least one formula is required for \'nlvnb2\' ")
    }else if(missing(form1)==F&missing(form2)==T){
      form <- form1
    }else if(missing(form1)==T&missing(form2)==F){
      form <- form2
    }else if(missing(form1)==F&missing(form2)==F){
      if(class(form1)=="formula"&class(form2)=="formula"){
        if(!quiet) warning("Two formulas were provided for \'nlvnb2\', unsing \'vnb2\' instead.")
        lh <- "vnb2"
        predictor <- NULL
      }else if(class(form1)=="formula"&class(form2)!="formula"){
        form <- form1
      }else if(class(form1)!="formula"&class(form2)=="formula"){
        form <- form2
      }else{
        stop("Only formula classes can be used.")
      }
    }
  }

  if(missing(dispersion.method)==T){
    dispersion.method <- "default"
  }else dispersion.method <- match.arg(dispersion.method)


  if(lambda!=1&lambda!=0){
    data$xTEMP <- data$x^lambda
  }else if(lambda==0){
    data$xTEMP <- log(data$x)
  }else if(lambda==1){
    data$xTEMP<-data$x
  }
  if(length(which(data$x==-Inf | data$x==Inf ) )!=0 ){
    data$xTEMP <- log(data$x+min(data$x[which(data$x>0)]))
  }
  data$x <- data$xTEMP

  if( lh!="skellam" & lh!="pskellam"){
    if(length(which(data$y==0))>0  ) data$y[which(data$y==0)] <- controls$zerocounts
  }

  if(theta.start[1]=="smart"){
    if(controls$allow.jitter==F){
      theta.start <- as.numeric( unname( iv.ames(data=data,predictor=predictor,link=link,
                                                 dose.power=dose.power,controls=controls,...) ) )
      return(fit.ames(data=data,
                      lh=lh,
                      predictor=predictor,
                      link=log,
                      dose.power=dose.power,
                      theta.start=theta.start,
                      dispersion.start=dispersion.start,
                      dispersion.method=dispersion.method,
                      controls=controls,
                      quiet=quiet))
    }else{
      theta.start <- as.numeric( unname( iv.ames(data=data,predictor=predictor,link=link,
                                                 dose.power=dose.power,controls=controls,...)))
      temp.mod <- try( fit.ames(data=data,
                                lh=lh,
                                predictor=predictor,
                                link=log,
                                dose.power=dose.power,
                                theta.start=theta.start,
                                dispersion.start=dispersion.start,
                                dispersion.method=dispersion.method,
                                controls=controls,
                                quiet=T),silent=T
      )
      cnd <- ifelse(class(temp.mod)=="try-error",T,is.null(temp.mod$coef) )
      max.jitter <- 0

      while(cnd==T && max.jitter <= controls$max.jitter){
        temp.mod <- try( fit.ames(data=data,
                                  lh=lh,
                                  predictor=predictor,
                                  link=log,
                                  dose.power=dose.power,
                                  theta.start=theta.start,
                                  dispersion.start=dispersion.start,
                                  dispersion.method=dispersion.method,
                                  controls=controls,
                                  quiet=T),silent=T )

        cnd <- ifelse(class(temp.mod)=="try-error",T,is.null(temp.mod$coef) )
        max.jitter <- max.jitter + 1
      }
      #cat(paste("Number of initial iterations:",max.jitter,sep=" "),sep="\n")

      return(temp.mod)
    }
  }

  if(lh=="skellam"|lh=="pskellam"|lh=="vnb2"){
    fit.method <- "vglm"
  }else if(lh=="nlvnb2"){
    fit.method <- "vgnlm"
  }else fit.method <- "gnlm"
  if(fit.method=="gnlm" | fit.method=="vgnlm"){
    switch(predictor,
           krewski = {
             alpha <- controls$alpha
             Eta <- function(x,theta) (theta[1]+theta[2]*x)*exp(-theta[3]*(x^alpha))
             J   <- function(x,theta){
               c1 <- exp(-theta[3]*(x^alpha))
               c2 <- x*exp(-theta[3]*(x^alpha))
               c3 <- -(x^alpha)*Eta(x,theta)
               return(cbind(c1,c2,c3))
             }
             coefnam <- c("beta0","beta1","gamma")
           },
           myers = {
             Eta <- function(x,theta) (theta[1]+theta[2]*x)*exp(-theta[3]*x)
             J   <- function(x,theta){
               c1 <- exp(-theta[3]*x)
               c2 <- x*exp(-theta[3]*x)
               c3 <- -x*(theta[1]+theta[2]*x)*exp(-theta[3]*x)
               return(cbind(c1,c2,c3))
             }
             coefnam <- c("beta0","beta1","gamma")

           },
           bernstein = {
             Eta <- function(x,theta) theta[1]+theta[2]*x
             J   <- function(x,theta){
               c1 <- rep(1,length(x))
               c2 <- x
               return(cbind(c1,c2))
             }
             coefnam <- c("beta0","beta1")
           },
           breslow = {
             delta <- controls$delta
             Eta <- function(x,theta) theta[1] + theta[2]*log((x+delta)) -theta[3]*x
             J   <- function(x,theta){
               c1 <- rep(1,length(x))
               c2 <- log(x+delta)
               c3 <- -x
               return(cbind(c1,c2,c3))
             }
             coefnam <- c("beta0","beta1","gamma")
           },
           stead = {
             if(length(which(data$x==0))>0) data$x[which(data$x==0)] <- 1e-3
             Eta <- function(x,theta) (theta[1]+theta[2]*(x^theta[3]))*exp(-theta[4]*x)
             J   <- function(x,theta){
               c1 <- exp(-theta[4]*x)
               c2 <- (x^theta[3])*exp(-theta[4]*x)
               c3 <- log(x)*theta[1]*(x^theta[3])*exp(-theta[4]*x)
               c4 <- -x*(theta[1]+theta[2]*(x^theta[3]))*exp(-theta[4]*x)
               return(cbind(c1,c2,c3,c4))
             }
             coefnam <- c("beta0","beta1","beta2","gamma")
           },
           margolin = {
             m <- link(controls$cfu)
             Eta <- function(x,theta) m*(1-exp(-(theta[1]+theta[2]*x)))*exp(-theta[3]*x)
             J   <- function(x,theta){
               c1 <- m*exp(-(theta[1]+theta[2]*x)-theta[3]*x)
               c2 <- m*x*exp(-(theta[1]+theta[2]*x)-theta[3]*x)
               c3 <- -x*m*(1-exp(-(theta[1]+theta[2]*x)))*exp(-theta[3]*x)
               return(cbind(c1,c2,c3))
             }
             coefnam <- c("beta0","beta1","gamma")
           }
    )
  }
  if(fit.method=="vglm" | fit.method=="vgnlm"){
    switch(lh,
           skellam={
             vary <- function(eta1,eta2,extra=NULL){
               mu1 <- exp(eta1)
               mu2 <- exp(eta2)
               return( mu1 + mu2 )
             }
             dl.detas <- function(y,eta1,eta2,method){
               if(missing(method)) method = 1
               if(method==1){
                 mu1 <- exp(eta1-ifelse(is.null(off1),0,off1))
                 mu2 <- exp(eta2-ifelse(is.null(off2),0,off2))
                 k   <- 2*sqrt(mu1*mu2)
                 num <- suppressWarnings( besselI(x=k,nu=y+1) + besselI(x=k,nu=y-1) )
                 den <- suppressWarnings( k*besselI(x=k,nu=y) )

                 d1  <- (mu2*num/den + (y/(2*mu1)) - 1)*mu1
                 d2  <- (mu1*num/den - (y/(2*mu2)) -1)*mu2
                 return(rbind(d1,d2))
               }
               else{
                 mu1 <- exp(eta1-ifelse(is.null(off1),0,off1))
                 mu2 <- exp(eta2-ifelse(is.null(off2),0,off2))
                 k   <- 2*sqrt(mu1*mu2)
                 num1 <- besselIs2(x=k,nu=abs(y+1),expon=T,log=T)
                 num2 <- besselIs2(x=k,nu=abs(y-1),expon=T,log=T)
                 den  <- besselIs2(x=k,nu=abs(y),expon=T,log=T)
                 div  <- exp(num1-den) + exp(num2-den)

                 d1  <- ( (mu2/k) * div + (y/(2*mu1)) - 1)*mu1
                 d2  <- ( (mu1/k) * div - (y/(2*mu2)) -1)*mu2
                 return(rbind(d1,d2))
               }
             }
             Wfun <- function(y,eta1,eta2,nsim){
               W0sim <- array(0,dim=c(2,2,length(y)))
               njsim <- 1
               while(njsim<=nsim){
                 ysim <- rskellam(length(y),mu1=exp(eta1),mu2=exp(eta2))
                 disim <- NULL
                 Wisim <- NULL
                 for(i in 1:length(y)){
                   disim[[i]] <- dl.detas(ysim[i],
                                          eta1[i],
                                          eta2[i],
                                          method=controls$detas.method)
                   Wisim[[i]] <- disim[[i]] %*% t( disim[[i]] )
                 }
                 cnd <- any(sapply(Wisim,is.nan,simplify=T),
                            sapply(Wisim,is.infinite,simplify=T),
                            sapply(Wisim,is.na,simplify=T)
                 )
                 if(!cnd){
                   for(i in 1:length(y)){
                     W0sim[,,i] <- W0sim[,,i] + Wisim[[i]]
                   }
                   njsim <- njsim + 1
                 }
               }
               Wi <- W0sim/nsim
               return(Wi)
             }
             Lfun <- function(y,eta1,eta2){
               Lret <- suppressWarnings( sum(dskellam(y,mu1=exp(eta1),mu2=exp(eta2),log=T)) )
               if(is.finite(Lret)==F){
                 Lret <- suppressWarnings( sum(dskellam2(y,mu1=exp(eta1),mu2=exp(eta2),log=T)) )
               }
               return(Lret)
             }
           },
           pskellam={
             vary <- function(eta1,eta2,extra){
               mu1 <- exp(eta1)
               mu2 <- exp(eta2)
               return( mu1 + mu2 + extra )
             }
             dl.detas <- function(y,eta1,eta2,method){
               if(missing(method)) method = 1
               if(method==1){
                 mu1 <- exp(eta1-ifelse(is.null(off1),0,off1))
                 mu2 <- exp(eta2-ifelse(is.null(off2),0,off2))
                 k   <- 2*sqrt(mu1*mu2)
                 num <- suppressWarnings( besselI(x=k,nu=y+1) + besselI(x=k,nu=y-1) )
                 den <- suppressWarnings( k*besselI(x=k,nu=y) )
                 d1  <- (mu2*num/den + (y/(2*mu1)) - 1)*mu1
                 d2  <- (mu1*num/den - (y/(2*mu2)) -1)*mu2
                 return(rbind(d1,d2))
               }
               else{
                 mu1 <- exp(eta1-ifelse(is.null(off1),0,off1))
                 mu2 <- exp(eta2-ifelse(is.null(off2),0,off2))
                 k   <- 2*sqrt(mu1*mu2)
                 num1 <- besselIs2(x=k,nu=abs(y+1),expon=T,log=T)
                 num2 <- besselIs2(x=k,nu=abs(y-1),expon=T,log=T)
                 den  <- besselIs2(x=k,nu=abs(y),expon=T,log=T)
                 div  <- exp(num1-den) + exp(num2-den)

                 d1  <- ( (mu2/k) * div + (y/(2*mu1)) - 1)*mu1
                 d2  <- ( (mu1/k) * div - (y/(2*mu2)) -1)*mu2
                 return(rbind(d1,d2))
               }
             }
             Wfun <- function(y,eta1,eta2,nsim){
               W0sim <- array(0,dim=c(2,2,length(y)))
               for(j in 1:nsim){
                 ysim <- rskellam(length(y),mu1=exp(eta1),mu2=exp(eta2))
                 disim <- NULL
                 Wisim <- NULL
                 for(i in 1:length(y)){
                   disim[[i]] <- dl.detas(ysim[i],
                                          eta1[i],
                                          eta2[i],
                                          method=controls$detas.method)
                   Wisim[[i]] <- disim[[i]] %*% t( disim[[i]] )
                   W0sim[,,i] <- W0sim[,,i] + Wisim[[i]]
                 }
               }
               Wi <- W0sim/nsim
               return(Wi)
             }
             Lfun <- function(y,eta1,eta2){
               Lret <- suppressWarnings( sum(dskellam(y,mu1=exp(eta1),mu2=exp(eta2),log=T)) )
               if(is.finite(Lret)==F){
                 Lret <- suppressWarnings( sum(dskellam2(y,mu1=exp(eta1),mu2=exp(eta2),log=T)) )
               }
               return(Lret)
             }
           },
           vnb2={
             vary <- function(eta1,eta2,extra=NULL){
               mu <- exp(eta1)
               k  <- exp(eta2)
               return( mu + (mu^2)/k)
             }
             dl.detas <- function(y,eta1,eta2,method){
               m <- exp(eta1)
               k <- exp(eta2)
               d1 <- y*(1-(m/(m+k))) - (m*k)/(m+k)
               d2 <- -y*(k/(m+k)) + k*(log(k)+1) - k*(log(m+k) + (k/(m+k))) +
                 k*digamma(y+k) - k*digamma(k)

               return(rbind(d1,d2))
             }
             Wfun0 <- function(y,eta1,eta2){
               m <- exp(eta1)
               k <- exp(eta2)

               a1 <- (1-(m/(m+k)))
               b1 <- (m*k)/(m+k)
               d1 <- (a1^2)*( m + (m^2) + (m^2)/k ) - 2*a1*b1*m + b1^2

               Pry <- function(obs,size,mean,q=0:10000){
                 return(
                   sum(
                     (1-pnbinom(q,size=size[obs],mu=mean[obs]) +
                        dnbinom(q,size=size[obs],mu=mean[obs])) /
                       ( (size[obs]+q)^2 )
                   )
                 )
               }
               MuP <- function(obs,size,mean){
                 return(
                   mean[obs]/(size[obs]*(mean[obs]+size[obs]) )
                 )
               }
               kp <- NULL
               for(i in 1:length(y)){
                 kp <- c(kp,(k[i]^2)*(Pry(i,k,m) - MuP(i,k,m)))
               }
               d2 <- sum(kp)

               return( rbind(c(d1,0),c(0,d2)) )
             }
             Wfun <- function(y,eta1,eta2,nsim){
               Wi <- array(0,dim=c(2,2,length(y)))
               for(i in 1:length(y)){
                 Wi[,,i] <- Wfun0(y[i],eta1[i],eta2[i])
               }
               return(Wi)
             }
             Lfun <- function(y,eta1,eta2){
               Lret <- sum(dnbinom(y,mu=exp(eta1),size=exp(eta2),log=T))
               return(Lret)
             }
           },
           nlvnb2={
             vary <- function(eta1,eta2,extra=NULL){
               mu <- exp(eta1)
               k  <- exp(eta2)
               return( mu + (mu^2)/k)
             }
             dl.detas <- function(y,eta1,eta2,method){
               m <- exp(eta1)
               k <- exp(eta2)
               d1 <- y*(1-(m/(m+k))) - (m*k)/(m+k)
               d2 <- -y*(k/(m+k)) + k*(log(k)+1) - k*(log(m+k) + (k/(m+k))) +
                 k*digamma(y+k) - k*digamma(k)

               return(rbind(d1,d2))
             }
             Wfun0 <- function(y,eta1,eta2){
               m <- exp(eta1)
               k <- exp(eta2)

               a1 <- (1-(m/(m+k)))
               b1 <- (m*k)/(m+k)
               d1 <- (a1^2)*( m + (m^2) + (m^2)/k ) - 2*a1*b1*m + b1^2

               Pry <- function(obs,size,mean,q=0:10000){
                 return(
                   sum(
                     (1-pnbinom(q,size=size[obs],mu=mean[obs]) +
                        dnbinom(q,size=size[obs],mu=mean[obs])) /
                       ( (size[obs]+q)^2 )
                   )
                 )
               }
               MuP <- function(obs,size,mean){
                 return(
                   mean[obs]/(size[obs]*(mean[obs]+size[obs]) )
                 )
               }
               kp <- NULL
               for(i in 1:length(y)){
                 kp <- c(kp,(k[i]^2)*(Pry(i,k,m) - MuP(i,k,m)))
               }
               d2 <- sum(kp)

               return( rbind(c(d1,0),c(0,d2)) )
             }
             Wfun <- function(y,eta1,eta2,nsim){
               Wi <- array(0,dim=c(2,2,length(y)))
               for(i in 1:length(y)){
                 Wi[,,i] <- Wfun0(y[i],eta1[i],eta2[i])
               }
               return(Wi)
             }
             Lfun <- function(y,eta1,eta2){
               Lret <- sum(dnbinom(y,mu=exp(eta1),size=exp(eta2),log=T))
               return(Lret)
             }
           }

    )
    thetaIT <- theta.start

    if(lh!="nlvnb2"){
      X1   <- unname( model.matrix(form1,data=data) )
      X2   <- unname( model.matrix(form2,data=data) )
      off1 <- unname( model.offset(model.frame(form1,data=data) ))
      off2 <- unname( model.offset(model.frame(form2,data=data) ))
      eta1I <- X1 %*% thetaIT[1:ncol(X1)] + ifelse(is.null(off1),0,off1)
      eta2I <- X2 %*% thetaIT[(ncol(X1)+1):(ncol(X1)+ncol(X2))] + ifelse(is.null(off1),0,off2)
    }else{
      X1 <- J(data$x,thetaIT)
      X2 <- unname( model.matrix(form,data=data) )
      off2 <- unname( model.offset(model.frame(form,data=data) ))
      eta1I <- Eta(data$x,thetaIT)
      eta2I <- X2 %*% thetaIT[(ncol(X1)+1):(ncol(X1)+ncol(X2))]
    }
    if(lh=="pskellam"){
      y2     <- data$y
      yhat   <- exp(eta1I) - exp(eta2I)
      offI   <- y2 - yhat
      data$y <- y2 - round(offI)
    }

    Xi <- NULL
    for(i in 1:length(data$y)){
      Xi[[i]] <- as.matrix( t( bdiag(list(X1[i,],X2[i,]) ) ) )
    }
    its <- 0
    if(controls$conv.crit=="theta"){
      theta0 <- rep(0,length(thetaIT))
      conv.tol <- sqrt( sum( (thetaIT - theta0)^2 ) )
    }else if(controls$conv.crit=="likelh"){
      l0 <- 0
      l1 <- Lfun(data$y,eta1I,eta2I)
      conv.tol <- abs(l1 - l0)
    }
    conv.chk  <- 0
    step.size <- ifelse(controls$halfstepping,controls$stepsize,1)
    while( (conv.tol > controls$convergence.tol && its<=controls$max.it) || conv.chk==0 ){
      theta0  <- thetaIT
      l0      <- Lfun(data$y,eta1I,eta2I)
      di <- NULL
      for(i in 1:length(data$y)){
        di[[i]] <- dl.detas(data$y[i],eta1I[i],eta2I[i],method=controls$detas.method)
      }
      Wi <- Wfun(data$y,eta1I,eta2I,nsim=controls$nsim)

      for(i in 1:dim(Wi)[3]){
        if( length(which( diag(Wi[,,i]) < controls$diagW.min )) > 0){
          diag(Wi[,,i])[which(diag(Wi[,,i]) < controls$diagW.min)] <- controls$diagW.min
        }
      }

      P1 <- matrix(0,ncol=ncol(Xi[[1]]),nrow=ncol(Xi[[1]]))
      P2 <- matrix(0,ncol=1,nrow=ncol(Xi[[1]]))
      for(i in 1:length(data$y)){
        P1  <- P1 + t(Xi[[i]]) %*% Wi[,,i] %*% Xi[[i]]
        P2  <- P2 + t(Xi[[i]]) %*% Wi[,,i] %*% ( solve(Wi[,,i]) %*% di[[i]] )
      }
      thetaIT <- theta0 + (old.step <- (solve(P1) %*% P2)/step.size )

      if(lh!="nlvnb2"){
        eta1I <- X1 %*% thetaIT[1:ncol(X1)] + ifelse(is.null(off1),0,off1)
        eta2I <- X2 %*% thetaIT[(ncol(X1)+1):(ncol(X1)+ncol(X2))] + ifelse(is.null(off2),0,off2)
      }else{
        X1 <- J(data$x,thetaIT)
        Xi <- NULL
        for(i in 1:length(data$y)){
          Xi[[i]] <- as.matrix( t( bdiag(list(X1[i,],X2[i,]) ) ) )
        }
        eta1I <- Eta(data$x,thetaIT)
        eta2I <- X2 %*% thetaIT[(ncol(X1)+1):(ncol(X1)+ncol(X2))] + ifelse(is.null(off2),0,off2)
      }

      l1     <- Lfun(data$y,eta1I,eta2I)
      if(controls$halfstepping){
        ltemp  <- l1
        if(its==0) new.step <- old.step
        new.step <- min(2*new.step, old.step)
        pits <- 0
        while(ltemp <= l0 && pits<=1000){
          thetaIT  <- theta0 + new.step
          new.step <- new.step/step.size
          if(lh!="nlvnb2"){
            eta1I <- X1 %*% thetaIT[1:ncol(X1)] + ifelse(is.null(off1),0,off1)
            eta2I <- X2 %*% thetaIT[(ncol(X1)+1):(ncol(X1)+ncol(X2))] + ifelse(is.null(off2),0,off2)
          }else{
            X1 <- J(data$x,thetaIT)
            Xi <- NULL
            for(i in 1:length(data$y)){
              Xi[[i]] <- as.matrix( t( bdiag(list(X1[i,],X2[i,]) ) ) )
            }
            eta1I <- Eta(data$x,thetaIT)
            eta2I <- X2 %*% thetaIT[(ncol(X1)+1):(ncol(X1)+ncol(X2))] + ifelse(is.null(off2),0,off2)
          }
          pits <- pits + 1
          ltemp     <- Lfun(data$y,eta1I,eta2I)
        }
        l1 <- ltemp
      }

      if(controls$conv.crit=="theta"){
        conv.tol <- sqrt( sum( (thetaIT - theta0)^2 ) )
      }else if(controls$conv.crit=="likelh"){
        conv.tol <- abs(l1 - l0)
      }
      if(lh=="pskellam"){
        yhatI  <- exp(eta1I) - exp(eta2I)
        offI   <- y2 - yhat
        data$y <- y2 - round(offI)
      }

      if(exists("saveits")) if(saveits < its) conv.chk <- 0
      if(conv.tol < controls$convergence.tol & conv.chk==0){
        conv.chk <- 1
        conv.tol <- 1
        saveits  <- its+1
      }
      its <- its+1
    }
    if(lh=="pskellam"){
      switch(controls$strain,
             ta100={
               ms  <- seq(log(70),log(200),length.out=1000)
             },
             ta98={
               ms  <- seq(log(20),log(50),length.out=1000)
             }
      )
      dms <- c()
      for(i in 1:length(ms)){
        switch(controls$strain,
               ta100={
                 dms <- c(dms, sum(dtrunc(round(offI),spec="pois",a=70,b=200,lambda=exp(ms[i]),log=T)) )
               },
               ta98={
                 dms <- c(dms, sum(dtrunc(round(offI),spec="pois",a=20,b=50,lambda=exp(ms[i]),log=T)) )
               }
        )
      }
      mhat  <- ms[which(dms==max(dms))]
      Lmax  <- dms[which(dms==max(dms))]
      L.cut <- Lmax - qchisq(.95,1)/2
      low.m <- min( ms[which( (1-pchisq(2*(Lmax-dms), df=1)) > 0.05)] )
      upp.m <- max( ms[which( (1-pchisq(2*(Lmax-dms), df=1)) > 0.05)] )
      profiles <- data.frame(Estimate=mhat,MaxLikelh=Lmax,Lower95=low.m,Upper95=upp.m)
      rownames(profiles) <- c("poisson.mean")
      SDy <- sqrt( vary(eta1I,eta2I,extra=exp(mhat)) )
    }else{
      profiles <- NULL
      SDy <- sqrt( vary(eta1I,eta2I) )
    }

    W1 <- NULL
    for(i in 1:dim(Wi)[3]){
      W1 <- rbind(W1,Wi[,,i])
    }

    ids <- matrix(1:(2*dim(Wi)[3]),ncol=2,byrow=T)
    W   <- matrix()
    for(i in 1:nrow(ids)){
      W <- bdiag(W,W1[ids[i,],1:2])
    }
    W <- as.matrix(W); W <- W[-1,]; W <- W[,-1]

    X <- NULL
    for(i in 1:dim(Wi)[3]){
      X <- rbind(X,rbind(c(X1[i,],matrix(0,nrow=1,ncol=ncol(X2)) ),
                         c(matrix(0,nrow=1,ncol=ncol(X1)),X2[i,]))
      )
    }
    fisher   <- t(X) %*% W %*% X
    vartheta <- solve( fisher )
    coefsd   <- sqrt( diag(vartheta) )

    work.resid <- NULL
    for(i in 1:dim(Wi)[3]){
      work.resid <- rbind(work.resid, t( solve( Wi[,,i] ) %*% di[[i]] ))
    }

    pear.resid <- NULL
    for(i in 1:dim(Wi)[3]){
      pear.resid <- rbind(pear.resid, t( sqrtm( Wi[,,i]) %*% work.resid[i,] ) )
    }

    eta <- cbind(eta1I,eta2I)
    zi  <- work.resid + eta
    ResSS <- 0
    for(i in 1:dim(Wi)[3]){
      ResSS <- ResSS + t( (zi - eta)[i,] ) %*% Wi[,,i] %*% (zi - eta)[i,]
    }
    ResDF <- (2*dim(Wi)[3]) - ncol(X)
    H <- chol(W) %*% X %*% solve(t(X)%*%W%*%X )%*%t(X)%*% t(chol(W))

    model <- data.frame(Likelihood=lh,Predictor=ifelse(lh=="nlvnb2",predictor,"model.formulas"),
                        Link=as.character(substitute(link)),
                        Dispersion=NA)

    if(lh!="nlvnb2"){
      model.formulas <- list(formula1=form1,formula2=form2)
    }else{
      model.formulas <- list(formula=form)
    }

    if(lh!="nlvnb2"){
      coefnam <-  c(paste( paste("theta",1:ncol(X1), sep=""), rep("eta1"), sep="."),
                    paste( paste("theta",1:ncol(X2), sep=""), rep("eta2"), sep=".")
      )

    }else{
      coefnam <- c(coefnam,paste( paste("theta",1:ncol(X2), sep=""), rep("eta2"), sep="."))
    }
    coef           <- data.frame(Param=coefnam,Estimate=thetaIT,Sd=coefsd)
    rownames(coef) <- 1:nrow(coef)

    if(exists("dataF")){
      if(length(ncol( as.matrix(dataF[,-c(1:2)])) ) > 0 ){
        nam.dataF <- names(dataF)[-c(1:2)]
        if(lh=="pskellam"){
          data.ans  <- data.frame(count=dataF[,1],dose=data$x,dataF[,-c(1:2)])
        }else{
          data.ans  <- data.frame(count=data$y,dose=data$x,dataF[,-c(1:2)])
        }
        names(data.ans) <- c(names(data.ans)[1:2],nam.dataF)
      }else{
        data.ans  <- data.frame(count=data$y,dose=data$x)
      }
      if(lh=="pskellam"){
        data.ans <- data.frame(data.ans,skellam=data$y,poisson=offI)
      }
    }else{
      data.ans <- data.frame(count=data$y,dose=data$x)
    }


    if(lh=="skellam"){
      fitted.values <- exp(eta[,1]) - exp(eta[,2])
    }else if(lh=="pskellam"){
      fitted.values <- exp(mhat) + exp(eta[,1]) - exp(eta[,2])
    }else{
      fitted.values <- exp(eta[,1])
    }
    resp.resid <- (data.ans$count - fitted.values)/SDy
    ans <- list(model=model,
                coef=coef,
                convergence=data.frame(Iterations=its,ResSS=ResSS,LogLikelihood=l1,ResDF=ResDF),
                predictors=eta,
                fitted.values=fitted.values,
                data=data.ans,
                hat=H,
                residuals=data.frame(working=work.resid,pearson=pear.resid,response=resp.resid),
                fisher=fisher,
                matrices=list(X=X,W=W,Xi=Xi,Wi=Wi,di=di,zi=zi,Wblocks=W1,V=(SDy^2)),
                fit.method=fit.method,
                model.formulas=model.formulas,
                profiles=profiles,
                original.data=dataF,
                controls=controls)
    class(ans) <- "fitames"
    if(quiet==F){
      print.fitames(ans)
    }
    return(invisible(ans))
  }
  if(fit.method=="gnlm"){
    switch(lh,
           poisson={
             if(dispersion.method=="default") dispersion.method <- "ml"
             V  <- function(x,theta) Mu(x,theta)
             dev <- function(y,x,theta){
               core1 <- y*log( y / Mu(x,theta) ) + Mu(x,theta) - y
               core2 <- 0
               core <- abs(core1+core2)
               if(is.finite(sum(c(core1,core2)))==F){resdev <- NA
               }else resdev <- sign(y-Mu(x,theta))*sqrt(2*core)
               return(resdev)
             }
             Up <- function(y,x,theta){return(0)}
             Lp <- function(y,x,theta){return(1)}
             Kp <- function(x,theta){return(NA)}
             f <- function(x,y,d,theta){return(x-1)}
             fall <- uniroot
             llik <- function(y,x,theta) dpois(y,lambda=Mu(x,theta))
             coefnam <- c(coefnam,"phi(fixed)")
           },
           nb2={
             if(dispersion.method=="default"){
               dispersion.method <- "ml"
             }else if(dispersion.method=="profile"){
               dispersion.method <- "ml"
             }
             V  <- function(x,theta) Mu(x,theta) + (Mu(x,theta)^2)/theta[np]
             dev <- function(y,x,theta){
               core1 <- theta[np]*log((Mu(x,theta)+theta[np])/(y+theta[np]))
               core2 <- y*log((y/Mu(x,theta))*((Mu(x,theta)+theta[np])/(y+theta[np])))
               core <- abs(core1+core2)
               if(is.finite(sum(c(core1,core2)))==F){resdev <- NA
               }else resdev <- sign(y-Mu(x,theta))*sqrt(2*core)
               return(resdev)
             }
             Up <- function(y,x,theta){
               return(
                 sum((digamma(theta[np]+y)-digamma(theta[np])-
                        ((y+theta[np])/(theta[np]+Mu(x,theta)))+
                        log(theta[np]/(theta[np]+Mu(x,theta))))+1)
               )
             }

             Lp <- function(y,x,theta){
               return(
                 sum(trigamma(theta[np]+y)+
                       ((y-2*Mu(x,theta)-theta[np])/((theta[np]+Mu(x,theta))^2))+
                       (1-theta[np]*trigamma(theta[np]))/theta[np])
               )
             }
             f <- function(x,y,d,theta){
               return(
                 sum( ((y-Mu(d,theta))^2)/(Mu(d,theta)*(1+Mu(d,theta)/x))) - length(y) + (np-1)
               )
             }
             fall <- uniroot
             llik <- function(y,x,theta) dnbinom(y,mu=Mu(x,theta),size=theta[np])
             Kp <- function(x,theta){
               x=x;theta=theta
               Pry <- function(obs,x,theta,q=0:10000){
                 return(
                   sum(
                     (1-pnbinom(q,size=theta[np],mu=Mu(x,theta)[obs]) +
                        dnbinom(q,size=theta[np],mu=Mu(x,theta)[obs])) /
                       ( (theta[np]+q)^2 )
                   )
                 )
               }
               MuP <- function(obs,x,theta){
                 return(
                   Mu(x,theta)[obs]/(theta[np]*(Mu(x,theta)[obs]+theta[np]) )
                 )
               }
               kp <- NULL
               for(i in 1:length(x)){
                 kp <- c(kp,(Pry(i,x,theta) - MuP(i,x,theta)))
               }
               return(sum(kp))
             }
             coefnam <- c(coefnam,"k(size)")
           },
           quasipoisson={
             if(dispersion.method=="default"){ dispersion.method <- "mpl"
             }else if(dispersion.method=="ml"){
               dispersion.method <- "mpl"
               if(!quiet) warning("ML method for quasi family doesn't work, using MQL.")
             }else if(dispersion.method=="profile"){
               dispersion.method <- "mpl"
             }
             V   <- function(x,theta) Mu(x,theta)
             dev <- function(y,x,theta,core.only=F){
               mu <- Mu(x,theta)
               core1 <- 2*(y*log(ifelse(y == 0, 1, y/mu)) - (y-mu) )
               core2 <- 0
               if(core.only==T) return(core1+core2)
               core <- abs(core1+core2)
               if(is.finite(sum(c(core1,core2)))==F){resdev <- NA
               }else resdev <- sign(y-Mu(x,theta))*sqrt(core)
               return(resdev)
             }
             f <- function(x,y,d,theta){
               return(
                 (1/x)*sum( ((y-Mu(d,theta))^2)/(Mu(d,theta))) - length(y) + (np-1)
               )
             }
             fall <- uniroot.all
             coefnam <- c(coefnam,"phi")
             llik <- function(y,x,theta) exp( -(1/2)*log(2*pi*theta[np]*Mu(x,theta)) -(1/2)*(1/theta[np])*dev(y,x,theta,core.only=T) )

           },
           quasipower={
             if(dispersion.method=="default"){dispersion.method <- "profile"
             }else if(dispersion.method=="ml"){
               dispersion.method <- "mpl"
               if(!quiet) warning("ML method for extended-quasi family doesn't work, using MPL.")
             }
             V   <- function(x,theta) Mu(x,theta)^theta[nq]
             dev <- function(y,x,theta,core.only=F){
               mu <- Mu(x,theta)
               if(theta[nq]==1){
                 core1 <- 2*(y*log(ifelse(y == 0, 1, y/mu)) - (y-mu) )
               }else if(theta[nq]==2){
                 core1 <- 2*( y/mu - log(ifelse(y ==0, 1, y/mu)) -1 )
               }else{
                 core1 <- (2*(y^(2-theta[nq]) - (2-theta[nq])*y*(mu^(1-theta[nq])) + (1-theta[nq])*(mu^(2-theta[nq]))))/( (1-theta[nq])*(2-theta[nq]) )
               }
               core2 <- 0
               if(core.only==T) return(core1+core2)
               core <- abs(core1+core2)
               if(is.finite(sum(c(core1,core2)))==F){resdev <- NA
               }else resdev <- sign(y-Mu(x,theta))*sqrt(core)
               return(resdev)
             }
             f <- function(x,y,d,theta){
               return(
                 (1/x)*sum( ((y-Mu(d,theta))^2)/(Mu(d,theta)^theta[nq])) - length(y) + (np-1)
               )
             }
             fall <- uniroot.all
             llik <- function(y,x,theta) exp( -(1/2)*log(2*pi*theta[np]*Mu(x,theta)^theta[nq]) -(1/2)*(1/theta[np])*dev(y,x,theta,core.only=T) )
             coefnam <- c(coefnam,"phi","tau")
           },
           hd={
             if(dispersion.method=="default"){dispersion.method <- "profile"
             }else if(dispersion.method=="ml"){
               dispersion.method <- "mpl"
               if(!quiet) warning("ML method for Hinde-Demetrio family doesn't work, using MPL instead.")
             }
             V   <- function(x,theta) Mu(x,theta) + (Mu(x,theta)^(theta[nq]+1))/theta[np]
             dev <- function(y,x,theta,core.only=F){
               core1 <- NULL
               mu <- Mu(x,theta)
               for(i in 1:length(y)){
                 core1 <- c(core1,2*integrate(function(t){
                   (y[i]-t)/( t + t^(theta[nq]+1)/theta[np])
                 },lower=y[i],upper=mu[i])$value
                 )
               }
               core2 <- 0
               if(core.only==T) return(core1+core2)
               core <- abs(core1+core2)
               if(is.finite(sum(c(core1,core2)))==F){resdev <- NA
               }else resdev <- sign(y-Mu(x,theta))*sqrt(core)
               return(resdev)
             }
             f <- function(x,y,d,theta){
               return(
                 sum( ((y-Mu(d,theta))^2)/(Mu(d,theta)*(1+ (Mu(d,theta)^theta[nq]) /x))) - length(y) + (np-1)
               )
             }
             fall <- uniroot
             fall2 <- uniroot.all
             llik <- function(y,x,theta) exp( -(1/2)*log(2*pi*(y + (1/theta[np])*y^(theta[nq]+1))) -(1/2)*dev(y,x,theta,core.only=T) )
             coefnam <- c(coefnam,"phi","tau")
           }
    )

    invlink <- exp
    Mu  <- function(x,theta) invlink(Eta(x,theta))
    w   <- function(x,theta) ((Mu(x,theta))^2)/V(x,theta)

    if(lh=="poisson") thetaIT <- c(theta.start,1)
    thetaIT <- c(theta.start,dispersion.start)
    np      <- ifelse(lh=="hd"||lh=="quasipower",length(thetaIT)-1,length(thetaIT) )
    nq      <- np + (length(dispersion.start)-1)

    if(dispersion.method!="profile"){
      devI <- 0
      devF <- dev(data$y,data$x,thetaIT)
      conv.tol <- abs( sum(devF)^2 - sum(devI)^2 )
      if(is.finite(conv.tol)==F){
        conv.tol <- controls$convergence.tol/10
        thetaIT    <- NULL
        if(!quiet) warning('Deviances are not all finite at the start')
      }
      tt <- 0
      while( conv.tol > controls$convergence.tol & tt <= controls$max.it){
        devI  <- dev(data$y,data$x,thetaIT)
        W <- diag(w(data$x,thetaIT))
        X <- J(data$x,thetaIT)
        A <- diag(1/Mu(data$x,thetaIT))
        z <- X%*%thetaIT[1:(np-1)] + A%*%(data$y-Mu(data$x,thetaIT))
        thetaIT[1:(np-1)] <- solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%z
        if(dispersion.method=="ml"){
          thetaIT[np] <- thetaIT[np] - Up(y=data$y,x=data$x,theta=thetaIT)/Lp(y=data$y,x=data$x,theta=thetaIT)
        }
        devF <- dev(data$y,data$x,thetaIT)
        if(is.na(sum(devF))==T || is.na(sum(devI))==T){
          conv.tol <- controls$convergence.tol/10
          thetaIT    <- NULL
          if(!quiet) warning(paste('Deviances are not all finite at interation',tt+1,sep=" "))
        }else conv.tol <-abs( sum(devF)^2 - sum(devI)^2 )
        tt <- tt+1
      }
    }
    else{
      prof <- NULL
      for(i in 1:length( controls$tau.profile )){
        modI <- try(fit.ames(count=data$y,dose=data$x,predictor=predictor,lh=lh,
                             theta.start=theta.start,dispersion.start=c(dispersion.start[1],controls$tau.profile[i]),
                             dispersion.method="mpl",quiet=T),silent=T)
        cnd  <- ifelse(class(modI)=="try-error",T,is.null(modI$coef))
        if(!cnd){
          prof <- c(prof, modI$convergence$LogLikelihood)
        }
      }
      tauh <- controls$tau.profile[which(prof==max(prof))]
      Lmax <- prof[which(prof==max(prof))]
      low.t <- min( controls$tau.profile[which( (1-pchisq(2*(Lmax-prof), df=1)) > 0.05)] )
      upp.t <- max( controls$tau.profile[which( (1-pchisq(2*(Lmax-prof), df=1)) > 0.05)] )
      profiles <- data.frame(Estimate=tauh,MaxLikelh=Lmax,Lower95=low.t,Upper95=upp.t)
      rownames(profiles) <- c("tau")

      modF <- fit.ames(count=data$y,dose=data$x,predictor=predictor,lh=lh,
                       theta.start=theta.start,dispersion.start=c(1,tauh),
                       dispersion.method="mpl",quiet=T)

      modF$model$Dispersion <- "profile"
      modF$profiles <- profiles
      if(!quiet) print(modF)
      return(modF)
    }

    if(is.null(thetaIT)==F){
      if(dispersion.method=="mpl"){
        dispsolve   <- suppressWarnings( try(fall(f,c(0,1e20),y=data$y,d=data$x,theta=thetaIT),silent=T) )
        if( class(dispsolve)=="try-error" ) dispsolve <- suppressWarnings( try(fall2(f,c(0,1e20),y=data$y,d=data$x,theta=thetaIT),silent=T))
        dispsolve   <- ifelse(is.atomic(dispsolve)==T,ifelse(length(dispsolve)>0,dispsolve[1],1),dispsolve$root)
        thetaIT[np] <- dispsolve
      }

      model <- data.frame(Likelihood=lh,Predictor=predictor,Link=as.character(substitute(link)),Dispersion=dispersion.method)
      devf  <- dev(data$y,data$x,thetaIT)
      Wf     <- diag(w(data$x,thetaIT))
      Xf     <- J(data$x,thetaIT)
      Af     <- diag(1/Mu(data$x,thetaIT))
      zf     <- X%*%thetaIT[1:(np-1)] + A%*%(data$y-Mu(data$x,thetaIT))

      ResDev <- sum(devf^2)
      ResDF  <- length(data$y) - nq + ifelse(lh=="poisson",1,0)
      loglik <- sum(log(llik(data$y,data$x,thetaIT)))
      H  <- sqrt(Wf) %*% Xf %*% solve(t(Xf)%*%Wf%*%Xf )%*%t(Xf)%*%sqrt(Wf)
      devfs <- devf/sqrt( 1- diag(H))

      Ftheta <- t(Xf)%*%Wf%*%Xf
      Fdisp  <- ifelse(dispersion.method=="mpl",NA,Kp(data$x,thetaIT))
      if(length(dispersion.start)>1){
        Fdisp <- c(Fdisp,NA)
        Fisher <- unname(rbind(
          cbind(Ftheta,rep(0,(np-1)),rep(0,(np-1))),
          rbind(c(rep(0,(np-1)),Fdisp[1],0),c(rep(0,(nq-1)),Fdisp[2])) ))
      }else{
        Fisher <- unname( rbind( cbind(Ftheta,rep(0,(np-1))),c(rep(0,(np-1)),Fdisp) ) )
      }
      coefsd <- c(sqrt(diag(solve(Ftheta))),sqrt(1/Fdisp))
      if(lh=="quasipoisson"|lh=="quasipower"){
        coefsd <- c(sqrt(thetaIT[np]*diag(solve(Ftheta))),sqrt(1/Fdisp))
      }

      coef   <- data.frame(Param=coefnam,Estimate=thetaIT,Sd=coefsd)
      rownames(coef) <- 1:nrow(coef)

      fitted.values <- Mu(data$x,thetaIT)
      predictors    <- Eta(data$x,thetaIT)
      Vf            <- ifelse(rep(lh=="quasipoisson",length(data$x)),
                              thetaIT[np]*V(data$x,thetaIT),V(data$x,thetaIT))
      PearsonRes    <- (data$y-fitted.values)/sqrt( Vf )
      PearsonStdRes <- (data$y-fitted.values)/sqrt( Vf*(1 - diag(H)) )

      WorkingRes <- (data$y-fitted.values)/fitted.values

      if(exists("dataF")==F){
        dataF <- NULL
      }

      ans <- list(model=model,
                  coef=coef,
                  convergence=data.frame(Iterations=tt,ResidualDeviance=ResDev,LogLikelihood=loglik,ResDF=ResDF),
                  predictors=predictors,
                  fitted.values=fitted.values,
                  data=data.frame(count=data$y,dose=data$x),
                  hat=H,
                  residuals=data.frame(working=WorkingRes,pearson=PearsonRes,pearson.std=PearsonStdRes,deviance=devf,deviance.std=devfs),
                  fisher=Fisher,
                  matrices=list(A=Af,X=Xf,W=Wf,V=Vf),
                  fit.method=fit.method,
                  model.formulas=NULL,
                  profiles=NULL,
                  original.data=dataF,
                  controls=controls
      )
      class(ans) <- "fitames"
      if(quiet==F){
        print.fitames(ans)
      }
      return(invisible(ans))
    }else return(list(coef=NULL))
  }
}
