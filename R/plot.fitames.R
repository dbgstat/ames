#' Residual plots and model diagnostics and for fitames objects
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


plot.fitames <- function(obj,
                         smooth=T,
                         deviance.plots=F,
                         std.deviance.plots=F,
                         pearson.plots=F,
                         std.pearson.plots=F,
                         qqplot=F,
                         lag.plot=F,lag=1,
                         diag.measures=F,
                         envelope=F,nsim=99,probs=c(0.005,0.995),
                         all=F,
                         quiet=T){
  if(class(obj)!="fitames") stop("Requires a fitames object.")

  if(all==T){
    ret.now <- plot.all(obj,nsim=nsim)
    return(ret.now)
  }

  rets  <- NULL
  x     <- obj$data$dose
  y     <- obj$data$count
  def.par <- par(no.readonly = TRUE)

  if(smooth==T & obj$fit.method!="vglm"){
    x.smooth  <- seq(min(x),max(x),length.out=10000)
    predictor <- as.character( obj$model$Predictor )
    link      <- as.character( obj$model$Link )
    theta.hat <- obj$coef$Estimate

    switch(predictor,
           krewski = {
             alpha <- obj$controls$alpha
             Eta <- function(x,theta) (theta[1]+theta[2]*x)*exp(-theta[3]*(x^alpha))
           },
           myers = {
             Eta <- function(x,theta) (theta[1]+theta[2]*x)*exp(-theta[3]*x)
           },
           bernstein = {
             Eta <- function(x,theta) theta[1]+theta[2]*x
           },
           breslow = {
             delta <- obj$controls$delta
             Eta <- function(x,theta) exp(theta[1])*((x+delta)^theta[2])*exp(-theta[3]*x)
           },
           stead = {
             Eta <- function(x,theta) (theta[1]+theta[2]*(x^theta[3]))*exp(-theta[4]*x)
           },
           svetliza = {
             Eta <- function(x,theta) exp(theta[1]-exp(theta[2]-theta[3]*x))
           },
           margolin = {
             m <- eval(call(link, obj$controls$cfu))
             Eta <- function(x,theta) m*(1-exp(-(theta[1]+theta[2]*x)))*exp(-theta[3]*x)
           }
    )
    y.smooth <- exp( Eta(x.smooth,theta.hat)  )

    layout(rbind(matrix(1,12,12,byrow=T),2))
    par(mar=c(4.5,4.5,1,1))
    plot(x,y,pch=19,cex=0.8,ylim=c(min(c(y,y.smooth)),max(c(y,y.smooth))),
         xlab=expression(dose[i]),ylab=expression(count[i]))
    lines(x.smooth,y.smooth,col=2)
    par(mar=c(0,0,0,0))
    plot(1:5,1:5,type='n',bty='n',xaxt='n',yaxt='n')
    legend(x=2.25,y=5,pch=c(19,NA,NA),lty=c(0,0,1),col=c(1,1,2),legend=c('observed','','predicted'),bty='n',ncol=3,cex=0.9)
  }else{
    layout(rbind(matrix(1,12,12,byrow=T),2))
    par(mar=c(4.5,4.5,1,1))
    plot(x,y,pch=19,cex=0.8,ylim=c(min(c(y,obj$fitted.values)),max(c(y,obj$fitted.values))),
         xlab=expression(dose[i]),ylab=expression(count[i]))
    points(x,obj$fitted.values,col=2)
    par(mar=c(0,0,0,0))
    plot(1:5,1:5,type='n',bty='n',xaxt='n',yaxt='n')
    legend(x=2.4,y=5,pch=c(19,NA,1),col=c(1,NA,2),legend=c('observed','','fitted'),bty='n',ncol=3,cex=0.9)
  }

  if(obj$fit.method!="gnlm"){
    x    <- obj$data$dose
    X    <- obj$matrices$X
    H    <- obj$hat
    Wi   <- obj$matrices$Wi
    nobs <- dim(Wi)[3]
    nij  <- sum(diag(H))

    ris.eta1 <- obj$residuals$pearson.1
    ris.eta2 <- obj$residuals$pearson.2
    rys      <- obj$residuals$response

    hi1  <- diag(H)[seq(1,2*nobs,2)]
    hi2  <- diag(H)[seq(2,2*nobs,2)]
    rstar.eta1 <- ris.eta1/sqrt(1-hi1)
    rstar.eta2 <- ris.eta2/sqrt(1-hi2)

    ress.min <- min( c(ris.eta1,ris.eta2,rstar.eta1,rstar.eta2,rys) )
    ress.max <- max( c(ris.eta1,ris.eta2,rstar.eta1,rstar.eta2,rys) )
    ylimits <- c( min(c(-2.5,ress.min)),max(c(2.5,ress.max)))

    if(deviance.plots==T | std.deviance.plots==T){
      warning("Deviance residuals are not defined for VG(N)LMs")
    }
    if(pearson.plots==T){
      par(def.par)
      par(mar=c(4.5,4.5,1,1))
      plot(x,rys,
           xlab=expression(dose[i]),ylab=expression(r[i]^'z'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(def.par)
      par(mar=c(4.5,4.5,1,1))
      plot(obj$fitted.values,rys,
           xlab=expression(hat(mu)[i]),ylab=expression(r[i]^'z'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(def.par)
      par(mar=c(4.5,4.5,1,1))
      plot(obj$predictors[,1],rys,
           xlab=expression(hat(eta)[i1]),ylab=expression(r[i]^'z'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(def.par)
      par(mar=c(4.5,4.5,1,1))
      plot(obj$predictors[,2],rys,
           xlab=expression(hat(eta)[i2]),ylab=expression(r[i]^'z'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(def.par)
      par(mar=c(4.5,4.5,1,1))
      plot(x,ris.eta1,
           xlab=expression(dose[i]),ylab=expression(r[i1]^'p'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(def.par)
      par(mar=c(4.5,4.5,1,1))
      plot(x,ris.eta2,
           xlab=expression(dose[i]),ylab=expression(r[i2]^'p'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(def.par)
      par(mar=c(4.5,4.5,1,1))
      plot(obj$fitted.values,ris.eta1,
           xlab=expression(hat(mu)[i]),ylab=expression(r[i1]^'p'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(def.par)
      par(mar=c(4.5,4.5,1,1))
      plot(obj$fitted.values,ris.eta2,
           xlab=expression(hat(mu)[i]),ylab=expression(r[i2]^'p'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(def.par)
      par(mar=c(4.5,4.5,1,1))
      plot(obj$predictors[,1],ris.eta1,
           xlab=expression(hat(eta)[i1]),ylab=expression(r[i1]^'p'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(def.par)
      par(mar=c(4.5,4.5,1,1))
      plot(obj$predictors[,2],ris.eta2,
           xlab=expression(hat(eta)[i2]),ylab=expression(r[i2]^'p'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(def.par)
      par(mar=c(4.5,4.5,1,1))
      plot(obj$predictors[,2],ris.eta1,
           xlab=expression(hat(eta)[i2]),ylab=expression(r[i1]^'p'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(def.par)
      par(mar=c(4.5,4.5,1,1))
      plot(obj$predictors[,1],ris.eta2,
           xlab=expression(hat(eta)[i1]),ylab=expression(r[i2]^'p'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)
    }
    if(std.pearson.plots==T){
      warning("Studentized pearson residuals are experimental for VG(N)LMs. Use carefully.")
      if(pearson.plots==F){
        par(def.par)
        par(mar=c(4.5,4.5,1,1))
        plot(x,rys,
             xlab=expression(dose[i]),ylab=expression(r[i]^'z'),
             ylim=ylimits,cex.axis=0.7,tcl=-0.3)

        par(def.par)
        par(mar=c(4.5,4.5,1,1))
        plot(obj$fitted.values,rys,
             xlab=expression(hat(mu)[i]),ylab=expression(r[i]^'z'),
             ylim=ylimits,cex.axis=0.7,tcl=-0.3)

        par(def.par)
        par(mar=c(4.5,4.5,1,1))
        plot(obj$predictors[,1],rys,
             xlab=expression(hat(eta)[i1]),ylab=expression(r[i]^'z'),
             ylim=ylimits,cex.axis=0.7,tcl=-0.3)
      }
      par(def.par)
      par(mar=c(4.5,4.5,1,1))
      plot(x,rstar.eta1,
           xlab=expression(dose[i]),ylab=expression(r[i1]^'p*'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(def.par)
      par(mar=c(4.5,4.5,1,1))
      plot(x,rstar.eta2,
           xlab=expression(dose[i]),ylab=expression(r[i2]^'p*'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(def.par)
      par(mar=c(4.5,4.5,1,1))
      plot(obj$fitted.values,rstar.eta1,
           xlab=expression(hat(mu)[i]),ylab=expression(r[i1]^'p*'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(def.par)
      par(mar=c(4.5,4.5,1,1))
      plot(obj$fitted.values,rstar.eta2,
           xlab=expression(hat(mu)[i]),ylab=expression(r[i2]^'p*'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(def.par)
      par(mar=c(4.5,4.5,1,1))
      plot(obj$predictors[,1],rstar.eta1,
           xlab=expression(hat(eta)[i1]),ylab=expression(r[i1]^'p*'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(def.par)
      par(mar=c(4.5,4.5,1,1))
      plot(obj$predictors[,2],rstar.eta2,
           xlab=expression(hat(eta)[i2]),ylab=expression(r[i2]^'p*'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(def.par)
      par(mar=c(4.5,4.5,1,1))
      plot(obj$predictors[,2],rstar.eta1,
           xlab=expression(hat(eta)[i2]),ylab=expression(r[i1]^'p*'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(def.par)
      par(mar=c(4.5,4.5,1,1))
      plot(obj$predictors[,1],rstar.eta2,
           xlab=expression(hat(eta)[i1]),ylab=expression(r[i2]^'p*'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)
    }

    if(envelope==F){
      if(qqplot==T){
        normq <- qnorm(seq(probs[1],probs[2],length.out=length(obj$fitted.values)))
        par(def.par)
        par(mar=c(4.5,4.5,1,1))
        plot(normq,sort(rys),
             xlab='Normal quantiles',
             ylab=expression(r[i]^'z'),
             ylim=ylimits,cex.axis=0.7,tcl=-0.3)
        abline(a=0,b=1,lty=2)
      }
    }else{
      if(!pearson.plots & !deviance.plots & !std.deviance.plots & !std.pearson.plots){
        if(obj$fit.method=="gnlm"){
          std.deviance.plots <- T
        }
      }
      par(def.par)
      env.obj <- envelope.fitames(obj=obj,nsim=nsim,
                                  ris.plot=pearson.plots,
                                  rstar.plot=std.pearson.plots,
                                  dis.plot=deviance.plots,
                                  dstar.plot=std.deviance.plots,
                                  eta.plot=NA)
      rets    <- c(rets,list(envelopes=env.obj))
    }

    if(diag.measures==T){
      Mi.1 <- max( c(max(hi1)+0.05, (2*nij/nobs)) )
      mi.1 <- min( c(max(hi1)+0.05, (2*nij/nobs)) )

      Mi.2 <- max( c(max(hi2)+0.05, (2*nij/nobs)) )
      mi.2 <- min( c(max(hi2)+0.05, (2*nij/nobs)) )

      par(def.par)
      par(mar=c(4.5,4.5,1,1))
      plot(1:nobs,hi1,ylim=c(0,Mi.2),
           xlab="i",ylab=expression(h[i1]),
           cex.axis=0.7,tcl=-0.3)
      abline(h=nij/nobs,lty=2)
      abline(h=2*nij/nobs,lty=2)

      par(def.par)
      par(mar=c(4.5,4.5,1,1))
      plot(1:nobs,hi2,ylim=c(0,Mi.2),
           xlab="i",ylab=expression(h[i2]),
           cex.axis=0.7,tcl=-0.3)
      abline(h=nij/nobs,lty=2)
      abline(h=2*nij/nobs,lty=2)

      par(def.par)
      diags.obj <- diagnostics.fitames(obj,eta.plot=NA)
      rets      <- c(rets,diags.obj)
    }
    par(def.par)
    if(quiet==T){
      return(invisible(rets))
    }
    else return(rets)

  }else{
    dstar <- obj$residuals$deviance.std
    dis   <- obj$residuals$deviance
    ris   <- obj$residuals$pearson
    rstar <- obj$residuals$pearson.std
    rys   <- obj$residuals$response
    ress.min <- min( c(dstar,dis,ris,rstar) )
    ress.max <- max( c(dstar,dis,ris,rstar) )
    ylimits <- c( min(c(-2.5,ress.min)),max(c(2.5,ress.max)))

    par(def.par)
    if(deviance.plots==T){
      par(mar=c(4.5,4.5,1,1))
      plot(x,dis,
           xlab=expression(dose[i]),ylab=expression(r[i]^'d'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(mar=c(4.5,4.5,1,1))
      plot(obj$fitted.values,dis,
           xlab=expression(hat(mu)[i]),ylab=expression(r[i]^'d'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(mar=c(4.5,4.5,1,1))
      plot(obj$predictors,dis,
           xlab=expression(hat(eta)[i]),ylab=expression(r[i]^'d'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)
    }
    if(std.deviance.plots==T){
      par(mar=c(4.5,4.5,1,1))
      plot(x,dstar,
           xlab=expression(dose[i]),ylab=expression(r[i]^'d*'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(mar=c(4.5,4.5,1,1))
      plot(obj$fitted.values,dstar,
           xlab=expression(hat(mu)[i]),ylab=expression(r[i]^'d*'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(mar=c(4.5,4.5,1,1))
      plot(obj$predictors,dstar,
           xlab=expression(hat(eta)[i]),ylab=expression(r[i]^'d*'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)
    }
    if(pearson.plots==T){
      par(mar=c(4.5,4.5,1,1))
      plot(x,ris,
           xlab=expression(dose[i]),ylab=expression(r[i]^'p'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(mar=c(4.5,4.5,1,1))
      plot(obj$fitted.values,ris,
           xlab=expression(hat(mu)[i]),ylab=expression(r[i]^'p'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(mar=c(4.5,4.5,1,1))
      plot(obj$predictors,ris,
           xlab=expression(hat(eta)[i]),ylab=expression(r[i]^'p'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)
    }
    if(std.pearson.plots==T){
      par(mar=c(4.5,4.5,1,1))
      plot(x,rstar,
           xlab=expression(dose[i]),ylab=expression(r[i]^'p*'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(mar=c(4.5,4.5,1,1))
      plot(obj$fitted.values,rstar,
           xlab=expression(hat(mu)[i]),ylab=expression(r[i]^'p*'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)

      par(mar=c(4.5,4.5,1,1))
      plot(obj$predictors,rstar,
           xlab=expression(hat(eta)[i]),ylab=expression(r[i]^'p*'),
           ylim=ylimits,cex.axis=0.7,tcl=-0.3)
    }

    if(lag.plot==T){
      if(deviance.plots==T){
        plot(dis[-c( (length(dis)-lag+1):length(dis) )],
             dis[-c(1:lag)],
             xlab=substitute(d[i-l], list(l=lag)),
             ylab=expression(r[i]^'d'),cex.axis=0.7,tcl=-0.3)
      }
      if(std.deviance.plots==T){
        plot(dstar[-c( (length(dstar)-lag+1):length(dstar) )],
             dstar[-c(1:lag)],
             xlab=substitute(d[i-l]^'*', list(l=lag)),
             ylab=expression(r[i]^'d*'),cex.axis=0.7,tcl=-0.3)
      }
      if(pearson.plots==T){
        plot(ris[-c( (length(ris)-lag+1):length(ris) )],
             ris[-c(1:lag)],
             xlab=substitute(r[i-l]^'p', list(l=lag)),
             ylab=expression(r[i]^'p'),cex.axis=0.7,tcl=-0.3)
      }
      if(std.pearson.plots==T){
        plot(rstar[-c( (length(rstar)-lag+1):length(rstar) )],
             rstar[-c(1:lag)],
             xlab=substitute(r[i-l]^'p*', list(l=lag)),
             ylab=expression(r[i]^'p*'),cex.axis=0.7,tcl=-0.3)
      }
    }

    if(envelope==F){
      if(qqplot==T){
        normq <- qnorm(seq(probs[1],probs[2],length.out=length(obj$fitted.values)))

        par(mar=c(4.5,4.5,1,1))
        plot(normq,sort(dis),
             xlab='Normal quantiles',
             ylab=expression(d^'(i)'),
             ylim=ylimits,cex.axis=0.7,tcl=-0.3)
        abline(a=0,b=1,lty=2)
      }
    }else{
      if(!pearson.plots & !deviance.plots & !std.deviance.plots & !std.pearson.plots){
        if(obj$fit.method=="gnlm"){
          std.deviance.plots <- T
        }
      }
      env.obj <- envelope.fitames(obj=obj,nsim=nsim,
                                  ris.plot=pearson.plots,
                                  rstar.plot=std.pearson.plots,
                                  dis.plot=deviance.plots,
                                  dstar.plot=std.deviance.plots)
      rets    <- c(rets,list(envelopes=env.obj))
    }

    if(diag.measures==T){
      hii <- diag(obj$hat)
      nii <- length(hii)
      pii <- length(obj$coef[,2])
      mii <- min(hii)
      Mii <- max(hii)

      plot(1:nii,hii,ylab=expression(h[i]),xlab='i',
           ylim=c(0,max( c(max(hii)+max(hii)/5,(2*pii)/nii +max(hii)/5) )))
      abline(h=(2*pii)/nii,lty=2)
      abline(h=(pii)/nii,lty=2)


      diags.obj <- diagnostics.fitames(obj)
      rets      <- c(rets,diags.obj)
    }

    par(def.par)
    if(quiet==T){
      return(invisible(rets))
    }
    else return(rets)
  }
}
