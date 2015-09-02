#' Simulated envelope in normal plots for fitames objects
#'
#' Used for computing simulated envelopes in ordered residual plots \emph{versus}
#' normal quantiles (normal plots). Although internally called within \code{\link{plot.fitames}}
#' given the right arguments, may be used in an independent fashion.
#'
#' @param obj A model fitted by \code{\link{fit.ames}}
#
#' @param nsim Integer, specifying the number of simulations
#'
#' @param ris.plot See \code{pearson.plots} of \code{\link{plot.fitames}}
#'
#' @param rstar.plot See \code{std.pearson.plots} of \code{\link{plot.fitames}}
#'
#' @param dis.plot See \code{deviance.plots} of \code{\link{plot.fitames}}
#'
#' @param dstar.plot See \code{std.deviance.plots} of \code{\link{plot.fitames}}
#'
#' @param nuisance.plot See \code{\link{diagnostics.fitames}}
#'
#' @param eta.plot See \code{\link{diagnostics.fitames}}
#'
#' @author Davi Butturi-Gomes
#'
#' Silvio S. Zocchi
#'
#' @examples
#' fit.ames()
#'
#' @seealso \code{\link{diagnostics.fitames}} \code{\link{plot.fitames}}
#'
#' @export

envelope.fitames <- function(obj,nsim,ris.plot=F,rstar.plot=F,dis.plot=F,dstar.plot=F,nuisance.plot=F, eta.plot=NULL){
  if(obj$fit.method=="gnlm"){
    if(!ris.plot & !rstar.plot & !dis.plot & !dstar.plot) dstar.plot <- T

    switch(as.character(obj$model$Likelihood),
           poisson={
             rsample <- function(fits,phi) rpois(length(fits),lambda=fits)
             if(as.character(obj$model$Dispersion)=='pearson') obj$model$Dispersion <- as.factor("ml")
           },
           nb2={
             rsample <- function(fits,phi) rnbinom(length(fits),mu=fits,size=phi)
           },
           quasipoisson={
             rsample <- function(fits,phi) rgenpois(1,mu=fits,phi=phi)
           }
    )
    def.par <- par(no.readonly = TRUE)
    thetasim <- obj$coef[,2]
    x.env <- obj$data$dose
    dstar <- obj$residuals$deviance.std
    dis   <- obj$residuals$deviance
    ris   <- obj$residuals$pearson
    rstar <- obj$residuals$pearson.std
    ress.min <- min( c(dstar,dis,ris,rstar) )
    ress.max <- max( c(dstar,dis,ris,rstar) )
    ylimits <- c( min(c(-2.5,ress.min)),max(c(2.5,ress.max)))
    di.env <- NULL
    ds.env <- NULL
    rp.env <- NULL
    rs.env <- NULL
    n.env <- nsim
    probs <- c(0.025,0.975)
    j <- 1

    while(j<=n.env){
      if(as.character(obj$model$Likelihood)=='quasipower'){
        y.env <- round( rtweedie(length(obj$fitted.values),mu=obj$fitted.values,
                                 phi=thetasim[(length(thetasim)-1)],
                                 xi=thetasim[length(thetasim)] ) )
      }else{
        y.env <- rsample(obj$fitted.values,thetasim[length(thetasim)])
      }
      if(length(which(y.env==0))>0) y.env[which(y.env==0)]<-obj$controls$zerocounts
      if(as.character(obj$model$Likelihood)=='quasipower'){
        obj.env <- suppressWarnings( try(fit.ames(count=y.env,dose=x.env,
                                                  lh=as.character(obj$model$Likelihood),
                                                  predictor=as.character(obj$model$Predictor),
                                                  dispersion.method=as.character(obj$model$Dispersion),
                                                  theta.start=thetasim[1:(length(thetasim)-2)],
                                                  dispersion.start=c(thetasim[(length(thetasim)-1)],
                                                                     thetasim[length(thetasim)]),
                                                  quiet=T),silent=T ))
      }else{
        obj.env <- suppressWarnings( try(fit.ames(count=y.env,dose=x.env,
                                                  lh=as.character(obj$model$Likelihood),
                                                  predictor=as.character(obj$model$Predictor),
                                                  dispersion.method=as.character(obj$model$Dispersion),
                                                  theta.start=thetasim[1:(length(thetasim)-1)],
                                                  dispersion.start=thetasim[length(thetasim)],quiet=T),silent=T))
      }
      cnd <- ifelse(class(obj.env)=="try-error",T,is.null(obj.env$coef))
      if(!cnd){
        di.env <- cbind(di.env, obj.env$residuals$deviance)
        ds.env <- cbind(ds.env, obj.env$residuals$deviance.std)
        rp.env <- cbind(rp.env, obj.env$residuals$pearson)
        rs.env <- cbind(rs.env, obj.env$residuals$pearson.std)
        j<-j+1
      }
    }

    env.st <- list(Deviance.envelope=apply(di.env,2,sort),
                   StdDeviance.envelope=apply(ds.env,2,sort),
                   Pearson.envelope=apply(rp.env,2,sort),
                   StdPearson.envelope=apply(rs.env,2,sort))

    env.di <- t(apply(env.st[[1]],1,quantile,probs=probs))
    env.ds <- t(apply(env.st[[2]],1,quantile,probs=probs))
    env.rp <- t(apply(env.st[[3]],1,quantile,probs=probs))
    env.rs <- t(apply(env.st[[4]],1,quantile,probs=probs))

    norm.ords <- 1:length(obj$fitted.values)
    a.normq   <- ifelse(length(obj$fitted.values)<=10,3/8,1/2)
    normq <- qnorm(  (norm.ords-a.normq)/(length(obj$fitted.values) + 1 -2*a.normq)   )
    
    if(dis.plot){
      par(mar=c(4.5,4.5,1,1))
      plot(normq,sort(dis),ylim=ylimits,
           xlab='Normal quantiles',ylab=expression(r[(i)]^'d'),
           cex.axis=0.7,tcl=-0.3)
      lines(normq,env.di[,1],lty=2)
      lines(normq,env.di[,2],lty=2)
      abline(a=0,b=1,lty=2)
    }

    if(dstar.plot){
      par(mar=c(4.5,4.5,1,1))
      plot(normq,sort(dstar),ylim=ylimits,
           xlab='Normal quantiles',ylab=expression(r[(i)]^'d*'),
           cex.axis=0.7,tcl=-0.3)
      lines(normq,env.ds[,1],lty=2)
      lines(normq,env.ds[,2],lty=2)
      abline(a=0,b=1,lty=2)
    }

    if(ris.plot){
      par(mar=c(4.5,4.5,1,1))
      plot(normq,sort(ris),ylim=ylimits,
           xlab='Normal quantiles',ylab=expression(r[(i)]^'p'),
           cex.axis=0.7,tcl=-0.3)
      lines(normq,env.rp[,1],lty=2)
      lines(normq,env.rp[,2],lty=2)
      abline(a=0,b=1,lty=2)
    }

    if(rstar.plot){
      par(mar=c(4.5,4.5,1,1))
      plot(normq,sort(rstar),ylim=ylimits,
           xlab='Normal quantiles',ylab=expression(r[(i)]^'p*'),
           cex.axis=0.7,tcl=-0.3)
      lines(normq,env.rs[,1],lty=2)
      lines(normq,env.rs[,2],lty=2)
      abline(a=0,b=1,lty=2)
    }
    return(invisible(env.st))
  }else{
    if(is.null(eta.plot)){
      eta.plot <- nuisance.plot
    }else{
      nuisance.plot <- eta.plot
    }

    def.par <- par(no.readonly = TRUE)
    lh   <- as.character( obj$model$Likelihood )
    X    <- obj$matrices$X
    x    <- obj$data$dose
    H    <- obj$hat
    Wi   <- obj$matrices$Wi
    nobs <- dim(Wi)[3]
    nij  <- sum(diag(H))
    X1 <- X[seq(1,nrow(X),2),]
    X1 <- as.matrix( X1[,-which( apply(X1,2,sum) == 0)] )
    X2 <- X[seq(2,nrow(X),2),]
    X2 <- as.matrix( X2[,-which( apply(X2,2,sum) == 0)] )
    theta.eta1 <- obj$coef$Estimate[1:ncol(X1)]
    theta.eta2 <- obj$coef$Estimate[(ncol(X1)+1):(ncol(X2)+ncol(X1))]
    switch(lh,
           skellam={
             rsample   <- function(n,eta1,eta2,extra=NULL) rskellam(n,mu1=exp(eta1),mu2=exp(eta2))
             form1.sim <- obj$model.formulas[[1]]
             form2.sim <- obj$model.formulas[[2]]
             pred.sim  <- "none"
             extra.sim <- NULL
             data.sim  <- obj$original.data
           },
           pskellam={
             rsample <- function(n,eta1,eta2,extra) rpois(n,lambda=exp(extra)) + rskellam(n,mu1=exp(eta1),mu2=exp(eta2))
             #rsample <- function(n,eta1,eta2,extra) extra + rskellam(n,mu1=exp(eta1),mu2=exp(eta2))
             form1.sim <- obj$model.formulas[[1]]
             form2.sim <- obj$model.formulas[[2]]
             pred.sim  <- "none"
             extra.sim <- obj$profiles$Estimate
             data.sim  <- obj$original.data
           },
           vnb2={
             rsample <- function(n,eta1,eta2,extra=NULL) rnbinom(n,mu=exp(eta1),size=exp(eta2))
             form1.sim <- obj$model.formulas[[1]]
             form2.sim <- obj$model.formulas[[2]]
             pred.sim  <- "none"
             extra.sim <- NULL
             data.sim  <- obj$original.data
           },
           nlvnb2={
             rsample <- function(n,eta1,eta2,extra=NULL) rnbinom(n,mu=exp(eta1),size=exp(eta2))
             form1.sim <- obj$model.formulas[[1]]
             form2.sim <- "none"
             pred.sim  <- as.character( obj$model$Predictor )
             extra.sim <- NULL
             data.sim  <- obj$original.data
           }
    )

    H     <- obj$hat
    hi1   <- diag(H)[seq(1,2*nobs,2)]
    hi2   <- diag(H)[seq(2,2*nobs,2)]
    rp.eta1 <- obj$residuals$pearson.1
    rp.eta2 <- obj$residuals$pearson.2
    rs.eta1 <- obj$residuals$pearson.1/sqrt(1-hi1)
    rs.eta2 <- obj$residuals$pearson.2/sqrt(1-hi2)
    rz      <- obj$residuals$response
    obj.env     <- NULL
    rp.eta1.env <- NULL
    rp.eta2.env <- NULL
    rs.eta1.env <- NULL
    rs.eta2.env <- NULL
    rz.env      <- NULL
    j <- 1
    while(j<=nsim){
      data.sim[,1]   <- rsample(nobs, eta1=obj$predictors[,1], eta2=obj$predictors[,2], extra=extra.sim)
      obj.env        <- suppressWarnings(try(fit.ames(data=data.sim,lh=lh,predictor=pred.sim,
                                                      form1=form1.sim,form2=form2.sim,
                                                      theta.start=c(theta.eta1,theta.eta2),
                                                      quiet=T,controls=obj$controls),silent=T))
      cnd <- ifelse(class(obj.env)=="try-error",T,is.null(obj.env$coef))
      if(!cnd){
        H.env <- obj.env$hat
        hi1   <- diag(H.env)[seq(1,2*nobs,2)]
        hi2   <- diag(H.env)[seq(2,2*nobs,2)]

        rp.eta1.env <- cbind(rp.eta1.env, obj.env$residuals$pearson.1)
        rp.eta2.env <- cbind(rp.eta2.env, obj.env$residuals$pearson.2)
        rs.eta1.env <- cbind(rs.eta1.env, obj.env$residuals$pearson.1/sqrt(1-hi1))
        rs.eta2.env <- cbind(rs.eta2.env, obj.env$residuals$pearson.2/sqrt(1-hi2))
        rz.env      <- cbind(rz.env, obj.env$residuals$response)
        j<-j+1
      }
    }
    env.st <- list(Pearson.envelope.eta1=apply(rp.eta1.env,2,sort),
                   Pearson.envelope.eta2=apply(rp.eta2.env,2,sort),
                   StdPearson.envelope.eta1=apply(rs.eta1.env,2,sort),
                   StdPearson.envelope.eta2=apply(rs.eta2.env,2,sort),
                   Response.envelope=apply(rz.env,2,sort)
    )

    if(dis.plot | dstar.plot){
      dis.plot   <- F
      dstar.plot <- F
      warning("Deviance residuals are not defined for VG(N)LMs")
    }
    if(!ris.plot & !rstar.plot & !dis.plot & !dstar.plot){
      zis.plot <- T
    }else zis.plot <- F

    probs   <- c(0.025,0.975)
    
    norm.ords <- 1:length(obj$fitted.values)
    a.normq   <- ifelse(length(obj$fitted.values)<=10,3/8,1/2)
    normq <- qnorm(  (norm.ords-a.normq)/(length(obj$fitted.values) + 1 -2*a.normq)   )
    
    env.ri.eta1 <- t(apply(env.st[[1]],1,quantile,probs=probs))
    env.ri.eta2 <- t(apply(env.st[[2]],1,quantile,probs=probs))
    env.rs.eta1 <- t(apply(env.st[[3]],1,quantile,probs=probs))
    env.rs.eta2 <- t(apply(env.st[[4]],1,quantile,probs=probs))
    env.rz <- t(apply(env.st[[5]],1,quantile,probs=probs))

    if((nuisance.plot==0 & !is.logical(nuisance.plot)) || zis.plot || is.na(nuisance.plot)){
#       par(def.par)
      par(mar=c(4.5,4.5,1,1))
      plot(normq,sort(rz),
           xlab='Normal quantiles',ylab=expression(r[(i)]^'z'),
           cex.axis=0.7,tcl=-0.3)
      lines(normq,env.rz[,1],lty=2)
      lines(normq,env.rz[,2],lty=2)
      abline(a=0,b=1,lty=2)
    }

    if(ris.plot){
      if((nuisance.plot==1 & !is.logical(nuisance.plot)) || (!nuisance.plot & is.logical(nuisance.plot)) || is.na(nuisance.plot)){
#         par(def.par)
        par(mar=c(4.5,4.5,1,1))
        plot(normq,sort(rp.eta1),
             xlab='Normal quantiles',ylab=expression(r[(i1)]^'p'),
             cex.axis=0.7,tcl=-0.3)
        lines(normq,env.ri.eta1[,1],lty=2)
        lines(normq,env.ri.eta1[,2],lty=2)
        abline(a=0,b=1,lty=2)
      }

      if((nuisance.plot==2 & !is.logical(nuisance.plot)) || (nuisance.plot & is.logical(nuisance.plot)) || is.na(nuisance.plot)){
#         par(def.par)
        par(mar=c(4.5,4.5,1,1))
        plot(normq,sort(rp.eta2),
             xlab='Normal quantiles',ylab=expression(r[(i2)]^'p'),
             cex.axis=0.7,tcl=-0.3)
        lines(normq,env.ri.eta2[,1],lty=2)
        lines(normq,env.ri.eta2[,2],lty=2)
        abline(a=0,b=1,lty=2)
      }
    }

    if(rstar.plot){
      if((nuisance.plot==1 & !is.logical(nuisance.plot)) || (!nuisance.plot & is.logical(nuisance.plot)) || is.na(nuisance.plot)){
#         par(def.par)
        par(mar=c(4.5,4.5,1,1))
        plot(normq,sort(rs.eta1),
             xlab='Normal quantiles',ylab=expression(r[(i1)]^'p*'),
             cex.axis=0.7,tcl=-0.3)
        lines(normq,env.rs.eta1[,1],lty=2)
        lines(normq,env.rs.eta1[,2],lty=2)
        abline(a=0,b=1,lty=2)
      }

      if((nuisance.plot==2 & !is.logical(nuisance.plot)) || (nuisance.plot & is.logical(nuisance.plot)) || is.na(nuisance.plot)){
#         par(def.par)
        par(mar=c(4.5,4.5,1,1))
        plot(normq,sort(rs.eta2),
             xlab='Normal quantiles',ylab=expression(r[(i2)]^'p*'),
             cex.axis=0.7,tcl=-0.3)
        lines(normq,env.rs.eta2[,1],lty=2)
        lines(normq,env.rs.eta2[,2],lty=2)
        abline(a=0,b=1,lty=2)
      }
    }
    return(invisible(env.st))
  }
}
