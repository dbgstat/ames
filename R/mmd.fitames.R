#' Maximum Mutagenic Dose (MMD)
#'
#' Computes the MMD values, given one or more previously fitted models
#' @param ... One or more "fitames" class objects, usually fitted obtained by fitting models using fit.ames().
#' @param model.names Character (optional), specifying names for each model.
#' @param akaike.weights Logical. If TRUE, internally calls aic.fitames(), recovering the Akaike weights and uses these values to compute the multimodel MMD. Defaults to FALSE.
#' @param emql.all Logical. Argument passed to aic.fitames(), if akaike.weights is TRUE, otherwise ignored. Defaults to FALSE.
#' @keywords mmd, ames, mi
#' @export
#' @examples
#' mmd.fitames()

mmd.fitames <- function(...,model.names,akaike.weights=F,emql.all=F){
  l.names   <- as.list(substitute(list(...)))[-1L]
  fits      <- list(...)
  nfit      <- length(fits)
  
  preds    <- NULL
  mmds     <- NULL
  fit.name <- NULL
  
  for(i in 1:nfit){
    if(class(fits[[i]])=="fitames"){
      preds    <- c(preds, as.character(fits[[i]]$model$Predictor) )
      switch(preds[i],
             bernstein={
               b1 <- fits[[i]]$coef$Estimate[2]
               if(b1 > 0){
                 mmds <- c(mmds, max(fits[[i]]$data$dose) )
               }else{
                 mmds <- c(mmds, min(fits[[i]]$data$dose) )
               }
             },
             breslow={
               delta <- fits[[i]]$controls$delta
               b0    <- fits[[i]]$coef$Estimate[1]
               b1    <- fits[[i]]$coef$Estimate[2]
               g     <- fits[[i]]$coef$Estimate[3]
               mmds  <- c(mmds, (b1-g*delta)/g)
             },
             krewski={
               alpha <- fits[[i]]$controls$alpha
               b0    <- fits[[i]]$coef$Estimate[1]
               b1    <- fits[[i]]$coef$Estimate[2]
               g     <- fits[[i]]$coef$Estimate[3]
               if(alpha==2){
                 mmds <- c(mmds,(-g*b0 + sqrt( (g*b0)^2 + 2*g*(b1)^2  ))/(2*g*b1))
               }else stop('Pending implementation')
             },
             margolin={
               link <- as.character( fits[[i]]$model$Link )
               m     <- fits[[i]]$controls$cfu
               b0    <- fits[[i]]$coef$Estimate[1]
               b1    <- fits[[i]]$coef$Estimate[2]
               g     <- fits[[i]]$coef$Estimate[3]
               mmds  <- c(mmds,(-1)*((b0 + log(g/(b1+g)))/b1))
             },
             myers={
               b0    <- fits[[i]]$coef$Estimate[1]
               b1    <- fits[[i]]$coef$Estimate[2]
               g     <- fits[[i]]$coef$Estimate[3]
               mmds <- c(mmds, (-1)*( (g*b0 - b1)/(g*b1) ) )
             },
             stead={
               b0    <- fits[[i]]$coef$Estimate[1]
               b1    <- fits[[i]]$coef$Estimate[2]
               b2    <- fits[[i]]$coef$Estimate[3]
               g     <- fits[[i]]$coef$Estimate[4]
               stead.der <- function(x){
                 ((b1*b2*(x^b2)*exp(-g*x))/x) - g*(b0+b1*(x^b2))*exp(-g*x)
               }
               xl <- c( min( fits[[i]]$data$dose ), max( fits[[i]]$data$dose ) )
               xstar <- try(uniroot(stead.der, lower = xl[1], upper = xl[2], tol = .Machine$double.eps^0.25),silent=T)
               if(class(xstar)=="try-error"){
                 xstar <- max(uniroot.all(stead.der, lower = xl[1], upper = xl[2], tol = .Machine$double.eps^0.25))
                 
               }else{
                 xstar <- xstar$root
               }
               mmds <- c(mmds,xstar)
             }
      )
      if(missing(model.names)) fit.name <- c(fit.name, as.character(l.names[[i]]))
    }else{
      stop('only fitames objects are allowed')
    }
  }
  if(!missing(model.names)) fit.name <- model.names
  if(akaike.weights){
    wis <- aic.fitames(...,emql.all=emql.all)$akaike.wis
    mmd <- sum(wis*mmds)
    return(mmd)
  }else{
    ret <- data.frame(model=fit.name,mmd=mmds)
    return(ret)
  }
}