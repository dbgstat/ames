#' Akaike measures for fitames objects
#'
#' Extracts the Akaike information criterion (AIC) of fitted models,
#' its second-order, bias-corrected version for small samples (AICc),
#' or extracts the QAIC/QAICc, which are based on the quasi-likelihood approach.
#' If several models are supplied, then Delta values are also
#' computed and plots the Akaike weights.
#'
#' @param ... One or more "fitames" class objects,
#' usually fitted obtained by fitting models using \code{\link{fit.fitames}}.
#' @param model.names Character (optional), specifying names for each model.
#' @param C Logical. If FALSE, computes AIC or QAIC
#' depening on the fitted model. If TRUE (default) computes the AICc or QAICc.
#' @param wis.plot Logical. If TRUE, creates an index plot using the
#' Akaike weights, and highlights those satisfying Delta < SEE. Defaults to FALSE.
#' @param emql.all Logical. If TRUE, computes the AIC or AICc based on the
#' extended quasi-likelihood function for each model. Defaults to FALSE.
#' @param SEE Numeric. Positive value indicating the strength of empyrical evidence (SEE).
#' Default value is 2.
#'
#' @keywords aic fitames
#'
#' @author Davi Butturi-Gomes
#'
#' Silvio S. Zocchi
#'
#' @examples
#' aic.fitames()
#'
#' @seealso \code{\link{fit.ames}}
#'
#' @export



aic.fitames <- function(...,model.names,C=T,wis.plot=F,emql.all=F,SEE=2){
  if(missing(model.names)) l.names <- as.list(substitute(list(...)))[-1L]
  fits      <- list(...)

  nfit <- length(fits)
  methods  <- NULL
  ndata    <- NULL
  family   <- NULL
  nparam   <- NULL
  aic      <- NULL
  fit.name <- NULL
  for(i in 1:nfit){
    if(class(fits[[i]])=="fitames"){
      methods  <- c(methods, as.character(fits[[i]]$model$Dispersion))
      ndata    <- c(ndata, nrow(fits[[i]]$data))
      family   <- c(family,as.character(fits[[i]]$model$Likelihood ) )
      nparam   <- c(nparam,ifelse(family[i]=="poisson",nrow(fits[[i]]$coef)-1,nrow(fits[[i]]$coef) ))
      if(missing(model.names)) fit.name <- c(fit.name, as.character( l.names[[i]] ) )
      if(!C){
        if(!emql.all){
          aic <- c(aic, -2*fits[[i]]$convergence$LogLikelihood + 2*nparam[i] )
        }else{

        }
      }else{
        b2.cor <- (2*nparam[i]*(nparam[i]+1))/(ndata[i] - nparam[i] - 1)
        if(!emql.all){
          aic <- c(aic, -2*fits[[i]]$convergence$LogLikelihood + 2*nparam[i] + b2.cor)
        }else{
          switch(family[i],
                 nb2={
                   mu <- fits[[i]]$fitted.values
                   k  <- fits[[i]]$coef$Estimate[length(fits[[i]]$coef$Estimate)]
                   y  <- fits[[i]]$data$count
                   emql <- sum( (-1/2)*log(2*pi*y*(1+(y/k))) + (log(mu)*y + log(mu+k)*(-y-k) -
                                                               log(y)*y - log(y+k)*(-y-k)   ) )
                 },
                 poisson={
                   mu <- fits[[i]]$fitted.values
                   y  <- fits[[i]]$data$count
                   emql <- sum( (-1/2)*log(2*pi*y) + (log(mu)*y-mu - log(y)*y + y ) )
                 },
                 quasipoisson={
                   emql <- fits[[i]]$convergence$LogLikelihood
                 },
                 quasipower={
                   emql <- fits[[i]]$convergence$LogLikelihood
                 }
          )
          aic <- c(aic, -2*emql + 2*nparam[i] + b2.cor)
        }
      }
    }else{
      stop('only fitames objects are allowed')
    }
  }

  if(!missing(model.names)) fit.name <- model.names

  if(emql.all){
    if(length(which(is.na(methods)==T))>0) stop('EMQL cannot be applied to VG(N)LMs')
  }else{
    cnd1.1 <- methods=='ml' | is.na(methods) | methods=='pearson'
    cnd1.2 <- methods=='profile' | methods=='mpl'
    if(length(is.na(cnd1.2))>0) cnd1.2[which(is.na(cnd1.2))] <- FALSE
    if(!(all(cnd1.1) | all(cnd1.2))) stop('Such models are not comparable with AIC or AICc.')
  }


  cnd2 <- sd(ndata)==0 || is.na(sd(ndata))
  if(!cnd2) stop('Models must be fitted to the same data.')

  aics       <- data.frame(model=fit.name,aic=aic,p=nparam)

  if(nfit>1){
    delta.aics <- aic - min(aic)
    akaike.wis <- exp(-delta.aics/2)/sum(exp(-delta.aics/2))
    ret0     <- data.frame(aics=aics,deltas=delta.aics,wis=akaike.wis)
    ret0     <- ret0[order(ret0$deltas) , ]
    ret0$SEE <- ifelse(ret0$deltas<=SEE,T,F)

    if(wis.plot){
      par(mar=c(4.5,4.5,1,1))
      plot(1:nfit,akaike.wis,type='h',xaxt='n',
           xlab="Models",ylab=expression(omega^(AIC)) )
      axis(1,at=1:nfit,labels=as.roman(1:nfit))
      lines(as.numeric(rownames(ret0[which(ret0$SEE==T), ])),ret0[which(ret0$SEE==T),]$wis,type='h',col=2 )
    }

    ret <- list(summary=ret0,aics=aics,delta.aics=delta.aics,akaike.wis=akaike.wis)
  }else{
    ret <- list(aics=aics,delta.aics=NULL,akaike.wis=NULL)
  }
  return(ret)
}
