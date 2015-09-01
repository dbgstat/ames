#' Model diagnostics for fitames objects
#'
#'
#' @param ...
#'
#'
#' @author Davi Butturi-Gomes
#'
#' Silvio S. Zocchi
#'
#' @seealso \code{\link{plot.fitames}}
#'
#' @examples
#' iv.ames()
#'
#'
#' @export



diagnostics.fitames <- function(obj,lmaxes.plot=T,curvs.plot=T,nuisance.plot=T,eta.plot=NULL){
  if(class(obj)!="fitames") stop("Requires a fitames object.")
  def.par <- par(no.readonly = TRUE)
  if(obj$fit.method=="gnlm"){
    lh        <- as.character( obj$model$Likelihood )
    predictor <- as.character(obj$model$Predictor)
    link      <- as.character( obj$model$Link )
    theta <- obj$coef$Estimate
    x     <- obj$data$dose
    y     <- obj$data$count
    mu    <- obj$fitted.values
    np    <- length(obj$coef$Estimate)
    ny    <- length(y)
    W  <- obj$matrices$W
    A  <- obj$matrices$A
    X  <- obj$matrices$X
    coefnam <- as.character(obj$coef$Param)

    switch(predictor,
           krewski={
             alpha <- obj$controls$alpha
             e11 <- 0;e12 <- 0;e13 <- -(x^alpha)*exp(-theta[3]*(x^alpha))
             e21 <- 0;e22 <- 0;e23 <- -(x^(alpha+1))*exp(-theta[3]*(x^alpha))
             e31 <- e13;e32 <- e23; e33 <- (x^(2*alpha))*(theta[1]+theta[2]*x)*exp(-theta[3]*(x^alpha))
             Z <- array( 0 ,dim=c(3,3,length(x)))
             Z[1,1,] <- e11;Z[1,2,] <- e21;Z[1,3,] <- e31
             Z[2,1,] <- e12;Z[2,2,] <- e22;Z[2,3,] <- e32
             Z[3,1,] <- e13;Z[3,2,] <- e23;Z[3,3,] <- e33
           },
           myers={
             e11 <- 0;e12 <- 0;e13 <- -x*exp(-theta[3]*x)
             e21 <- 0;e22 <- 0;e23 <- -(x^2)*exp(-theta[3]*x)
             e31 <- e13;e32 <- e23; e33 <- (x^2)*(theta[1]+theta[2]*x)*exp(-theta[3]*x)
             Z <- array( 0 ,dim=c(3,3,length(x)))
             Z[1,1,] <- e11;Z[1,2,] <- e21;Z[1,3,] <- e31
             Z[2,1,] <- e12;Z[2,2,] <- e22;Z[2,3,] <- e32
             Z[3,1,] <- e13;Z[3,2,] <- e23;Z[3,3,] <- e33
           },
           bernstein={
             Z <- array( 0 ,dim=c(2,2,length(x)))
           },
           breslow={
             delta <- obj$controls$delta
             Z <- array( 0 ,dim=c(3,3,length(x)))
           },
           stead={
             e11 <- 0;e12 <- 0;e13 <- 0;e14 <- -x*exp(-theta[4]*x)
             e21 <- 0;e22 <- 0;e23 <- (x^theta[3])*exp(-theta[4]*x)*log(x);e24 <- -(x^(theta[3]+1))*exp(-theta[4]*x)
             e31 <- 0;e32 <- e23; e33 <- 2*theta[2]*e23; e34 <- theta[2]*log(x)*e24
             e41 <- e14; e42 <- e24; e43 <- e34; e44 <- (x^2)*(theta[1]+theta[2]*(x^theta[3]))*exp(-theta[4]*x)
             Z <- array( 0 ,dim=c(4,4,length(x)))
             Z[1,1,] <- e11;Z[1,2,] <- e21;Z[1,3,] <- e31;Z[1,4,] <- e41
             Z[2,1,] <- e12;Z[2,2,] <- e22;Z[2,3,] <- e32;Z[2,4,] <- e42
             Z[3,1,] <- e13;Z[3,2,] <- e23;Z[3,3,] <- e33;Z[3,4,] <- e43
             Z[4,1,] <- e13;Z[4,2,] <- e23;Z[4,3,] <- e33;Z[4,4,] <- e44
           },
           margolin={
             m <- eval(call(link, obj$controls$cfu))
             c1 <- -m*exp(-theta[1]-theta[2]*x-theta[3]*x)
             e11 <- c1;e12 <- x*c1; e13 <- x*c1
             e21 <- e12;e22 <- (x^2)*c1;e23 <- (x^2)*c1
             e31 <- e13;e32 <- e23; e33 <- m*(x^2)*exp(-theta[3]*x)*(1-exp(-theta[1]-theta[2]*x))
             Z <- array( 0 ,dim=c(3,3,length(x)))
             Z[1,1,] <- e11;Z[1,2,] <- e21;Z[1,3,] <- e31
             Z[2,1,] <- e12;Z[2,2,] <- e22;Z[2,3,] <- e32
             Z[3,1,] <- e13;Z[3,2,] <- e23;Z[3,3,] <- e33
           }
    )

    switch(lh,
           poisson={
             V  <- diag(mu)
             Set <- (y-mu)/mu^2
             Net <- y-mu
             coefnam <- coefnam[-np]
             np  <- np-1
             ofi <- np+1
             labs  <- c(expression(theta),coefnam)
             labs2 <- expression(paste(C[i],"(",bold(theta),")"))
           },
           nb2={
             k  <- theta[np]
             V  <- diag(mu*(1+mu/k ))
             Set <- (k*(2*mu+k)*(y-mu)/((mu^2)*((mu+k)^2)))
             Net <- (y-mu)/( 1 + (mu/k) )
             ofi  <- np+2
             coefnam[np] <- "k"
             labs  <- c(expression(paste(theta,",k")),coefnam,expression(bold(theta)))
             labs2 <- c(expression(paste(C[i],"(",bold(theta),",k)")),expression(paste(C[i],"(",bold(theta),")")))
           },
           quasipoisson={
             V  <- theta[np]*diag(mu)
             Set <- (y-mu)/mu^2
             Net <- y-mu
             coefnam <- coefnam[-np]
             np  <- np-1
             ofi <- np+1
             labs  <- c(expression(theta),coefnam)
             labs2 <- expression(paste(C[i],"(",bold(theta),")"))
           },
           quasipower={
             V  <- theta[(np-1)]*diag((mu^theta[np]))
             Set <- (y-mu)/mu^2
             Net <- y-mu
             coefnam <- coefnam[-c((np-1),np)]
             np  <- np-2
             ofi <- np+1
             labs  <- c(expression(theta),coefnam)
             labs2 <- expression(paste(C[i],"(",bold(theta),")"))
           }
    )

    rA <- t(((y-mu)/diag(V)) %*% solve(A))
    SS <- diag(Set)
    NN <- diag(Net)

    Lbb <- t(X)%*%( -(W + (solve(A))^2 %*% SS - NN) )%*%X + do.call(cbind,lapply(seq_len(dim(Z)[1]),function(i) Z[i,,]%*%rA))
    Delta.theta <- t(X)%*%diag( (y-mu)/diag(V) )%*% solve(A)

    switch(lh,
           poisson={
             L     <- Lbb
             Delta <- Delta.theta
           },
           nb2={
             Lkk <- sum(trigamma(k+y)+((y-2*mu-k)/((k+mu)^2)))+ny*((1-k*trigamma(k))/k )
             Lbk <- t(X)%*%diag( 1/((mu+k)^2) )%*%solve(A)%*%(y-mu);Lkb <- t(Lbk)
             Delta.k <- digamma(k+y)-digamma(k)-((k+y)/(k+mu))+log(k/(k+mu))+1
             L     <- rbind(cbind(Lbb,Lbk),cbind(Lkb,Lkk))
             Delta <- rbind(Delta.theta,Delta.k)
             B11    <- matrix(0,ncol=ncol(L),nrow=nrow(L));B11[nrow(B11),ncol(B11)] <- solve(Lkk)
             B1     <- t(Delta) %*% (solve(L)-B11) %*% Delta
             lmax1  <- eigen(B1)$vectors[,which(abs(eigen(B1)$values)==max(abs(eigen(B1)$values)))]
           },
           quasipoisson={
             L     <- Lbb
             Delta <- Delta.theta
           },
           quasipower={
             L     <- Lbb
             Delta <- Delta.theta
           }
    )

    B    <- t(Delta) %*% solve(L) %*% Delta
    lmax <- eigen(B)$vectors[,which(abs(eigen(B)$values)==max(abs(eigen(B)$values)))]

    npB <- NULL
    for(i in 1:nrow(L)) npB <- cbind(npB, c(i,1:nrow(L))[ (-i-1) ] )
    Bs <- NULL
    for(i in 1:nrow(L)){
      Li <- L[npB[,i],npB[,i]]
      Bi  <- rbind( rep(0,np), cbind( rep(0,(np-1)), solve( Li[2:nrow(L),2:nrow(L)] ) ) )
      Bs[[i]]  <- t(Delta[npB[,i],]) %*% (solve(Li)-Bi) %*% Delta[npB[,i],]
    }

    lmaxes <- NULL
    for(i in 1:length(Bs)){
      lmaxes[[i]]  <- eigen(Bs[[i]])$vectors[,which(abs(eigen(Bs[[i]])$values)==max(abs(eigen(Bs[[i]])$values)))]
    }

    if(lmaxes.plot){
      par(mar=c(4.5,4.5,1,1))
      plot(rep(1,ny),abs(lmax),xlim=c(1,ofi),ylim=c(0,1),pch=19,cex=0.3,xaxt='n',
           ylab=expression(paste("|",l[max],"|",sep="")),cex.axis=0.7,tcl=-0.3,
           xlab='Parameters')
      if( length(which(abs(lmax)>0.4))>0 ){
        text(x=1.1,y=abs(lmax[which(abs(lmax)>0.4)] ),
             labels=which(abs(lmax)>0.4),cex=0.6)
      }
      for(i in 1:length(lmaxes)){
        points(rep((i+1),ny),abs(lmaxes[[i]]),xlim=c(0,4),ylim=c(0,1),pch=19,cex=0.3)
        if( length(which(abs(lmaxes[[i]]) > 0.4))>0){
          text(x=(1.1+i),y=abs(lmaxes[[i]][which(abs(lmaxes[[i]])>0.4)]),
               labels=which(abs(lmaxes[[i]]) > 0.4),cex=0.6)
        }
      }
      if(exists("lmax1")==T){
        points(rep((length(lmaxes)+2),ny),abs(lmax1),xlim=c(0,4),ylim=c(0,1),pch=19,cex=0.3)
        if(length(which(abs(lmax1) > 0.4)) > 0){
          text(x=(1.1+length(lmaxes)+1),y=abs( lmax1[which(abs(lmax1)>0.4)] ),
               labels=which(abs(lmax1) > 0.4),cex=0.6)
        }
      }
      axis(1,at=1:ofi,labels=labs,cex.axis=0.7,tcl=-0.3,)
    }

    Curv <- NULL
    for(i in 1:length(y)){
      li <- rep(0,ny); li[i] <- 1
      Curv <- rbind(Curv, 2*abs( t(li) %*% t(Delta) %*% solve(L) %*% Delta %*% li) )
    }
    cutoff <- sum(Curv)/(2*np)

    if(exists("B1")==T){
      Curv.b1 <- NULL
      for(i in 1:length(y)){
        li.b1 <- rep(0,ny); li.b1[i] <- 1
        Curv.b1 <- rbind(Curv.b1, 2*abs( t(li.b1) %*% t(Delta) %*% (solve(L)-B11) %*% Delta %*% li.b1) )
      }
      cutoff.1 <- sum(Curv.b1)/(2*(np-1))

      curvmaxes <- max(c(Curv,Curv.b1))
      cutmaxes  <- max(c(cutoff,cutoff.1))
      curvlim   <- max(c( (curvmaxes + curvmaxes/5), (cutmaxes+0.1)))

      if(curvs.plot){
        if(nuisance.plot){
          par(mar=c(4.5,4.5,1,1))
          plot(1:ny,Curv,ylim=c(0,curvlim),
               xlab='i',pch=19,cex=0.7,
               cex.axis=0.7,tcl=-0.3,
               ylab=labs2[1])
          abline(h=cutoff,lty=2)
        }

        par(mar=c(4.5,4.5,1,1))
        plot(1:ny,Curv.b1,ylim=c(0,curvlim),
             xlab='i',pch=19,cex=0.7,
             cex.axis=0.7,tcl=-0.3,
             ylab=labs2[2])
        abline(h=cutoff.1,lty=2)
      }
    }else{
      curvlim <- max(c((max(Curv)+(max(Curv)/5)), (cutoff + 0.1) ))
      if(curvs.plot){
        par(mar=c(4.5,4.5,1,1))
        plot(1:ny,Curv,ylim=c(0,curvlim),
             xlab='i',pch=19,cex=0.7,
             cex.axis=0.7,tcl=-0.3,
             ylab=labs2[1])
        abline(h=cutoff,lty=2)
      }
    }

    rownames(Delta) <- coefnam
    rownames(L) <- coefnam; colnames(L) <- coefnam

    if(exists("B11")==T){
      return(invisible(list(Delta=Delta,L=L,B1=B11,Ci=data.frame(Ci.all=Curv,Ci.theta=Curv.b1) )))
    }else return(invisible(list(Delta=Delta,L=L,Ci=Curv)))
  }else if(obj$fit.method=="vglm" | obj$fit.method=="vgnlm"){
    if(is.null(eta.plot)){
      eta.plot <- nuisance.plot
    }else{
      nuisance.plot <- eta.plot
    }

    H  <- obj$hat
    X  <- obj$matrices$X
    W  <- obj$matrices$W
    Xi <- obj$matrices$Xi
    Wi <- obj$matrices$Wi
    di <- obj$matrices$di

    thetaI <- obj$coef$Estimate

    X1 <- obj$matrices$X[seq(1,nrow(obj$matrices$X),2),]
    X1 <- as.matrix( X1[,-which( apply(X1,2,sum) == 0)] )
    X2 <- obj$matrices$X[seq(2,nrow(obj$matrices$X),2),]
    X2 <- as.matrix( X2[,-which( apply(X2,2,sum) == 0)] )

    np <- ncol(X)

    L <- matrix(0,ncol=length(thetaI),nrow=length(thetaI))
    for(i in 1:length(di)){
      L <- L + t(Xi[[i]]) %*% (di[[i]] %*% t(di[[i]])) %*% Xi[[i]]
    }
    Delta1 <- NULL
    Delta2 <- NULL
    for(i in 1:length(di)){
      Delta1 <- rbind(Delta1, X1[i,] * di[[i]][1] )
      Delta2 <- rbind(Delta2, X2[i,] * di[[i]][2] )
    }

    Delta <- t(cbind(Delta1,Delta2))

    B <- t(Delta) %*% solve(L) %*% Delta
    lmax <- eigen(B)$vectors[,which(abs(eigen(B)$values)==max(abs(eigen(B)$values)))]
    #
    B1.eta1 <- rbind(
      cbind( matrix(0,ncol(X1),ncol(X1)),matrix(0,ncol(X1),ncol(X2))),
      cbind( t(matrix(0,ncol(X1),ncol(X2))),
             solve(L[(ncol(X1)+1):(ncol(X1)+ncol(X2)),(ncol(X1)+1):(ncol(X1)+ncol(X2))])
      )
    )
    BB1 <- B1.eta1
    B.eta1 <- t(Delta) %*% (solve(L)-B1.eta1) %*% Delta
    lmax.eta1 <- eigen(B.eta1)$vectors[,which(abs(eigen(B.eta1)$values)==max(abs(eigen(B.eta1)$values)))]

    B1.eta2 <- rbind(
      cbind(solve( L[1:ncol(X1),1:ncol(X1) ] ),
            matrix(0,ncol(X1),ncol(X2))
      ),
      cbind(t(matrix(0,ncol(X1),ncol(X2))),
            matrix(0,ncol(X2),ncol(X2))
      )
    )
    BB2 <- B1.eta2
    B.eta2 <- t(Delta) %*% (solve(L)-B1.eta2) %*% Delta
    lmax.eta2 <- eigen(B.eta2)$vectors[,which(abs(eigen(B.eta2)$values)==max(abs(eigen(B.eta2)$values)))]

    ny  <- length(lmax)
    ofi <- ncol(X) + length(di[[1]]) + 1

    npB <- NULL
    for(i in 1:nrow(L)) npB <- cbind(npB, c(i,1:nrow(L))[ (-i-1) ] )

    Bs <- NULL
    for(i in 1:nrow(L)){
      Li <- L[npB[,i],npB[,i]]
      Bi  <- rbind( rep(0,np), cbind( rep(0,(np-1)), solve( Li[2:nrow(L),2:nrow(L)] ) ) )
      Bs[[i]]  <- t(Delta[npB[,i],]) %*% (solve(Li)-Bi) %*% Delta[npB[,i],]
    }

    lmaxes <- NULL
    for(i in 1:length(Bs)){
      lmaxes[[i]]  <- eigen(Bs[[i]])$vectors[,which(abs(eigen(Bs[[i]])$values)==max(abs(eigen(Bs[[i]])$values)))]
    }

    if(lmaxes.plot){
      par(mar=c(4.5,4.5,1,1))
      plot(rep(1,ny),abs(lmax),xlim=c(1,ofi),ylim=c(0,1),pch=19,cex=0.3,xaxt='n',
           ylab=expression(paste("|",l[max],"|",sep="")),cex.axis=0.7,tcl=-0.3,
           xlab='Parameters',cex.axis=0.7,tcl=-0.3)
      axis(1, at=1:ofi,labels=c(expression(bold(theta)),
                                as.character(obj$coef$Param),
                                expression(bold(theta[1])),
                                expression(bold(theta[2]))
      ),cex.axis=0.7,tcl=-0.3
      )

      if( length( which(abs(lmax)>0.4) )>0){
        text(x=1.1,y=abs(lmax[which(abs(lmax)>0.4)] ),
             labels=which(abs(lmax)>0.4),cex=0.6)
      }
      for(i in 1:length(lmaxes)){
        points(rep((i+1),ny),abs(lmaxes[[i]]),xlim=c(0,4),ylim=c(0,1),pch=19,cex=0.3)
        if( length(which(abs(lmaxes[[i]]) > 0.4))>0){
          text(x=(1.1+i),y=abs(lmaxes[[i]][which(abs(lmaxes[[i]])>0.4)]),
               labels=which(abs(lmaxes[[i]]) > 0.4),cex=0.6)
        }
      }

      points(rep((length(lmaxes)+2),length(lmax.eta1)),abs(lmax.eta1),xlim=c(0,4),ylim=c(0,1),pch=19,cex=0.3)
      if( length( which(abs(lmax.eta1)>0.4) )>0){
        text(x=(1.1+length(lmaxes)+1),y=abs( lmax.eta1[which(abs(lmax.eta1)>0.4)] ),
             labels=which(abs(lmax.eta1) > 0.4),cex=0.6)
      }
      points(rep((length(lmaxes)+3),length(lmax.eta2)),abs(lmax.eta2),xlim=c(0,4),ylim=c(0,1),pch=19,cex=0.3)
      if( length( which(abs(lmax.eta2)>0.4) )>0){
        text(x=(1.1+length(lmaxes)+2),y=abs( lmax.eta2[which(abs(lmax.eta2)>0.4)] ),
             labels=which(abs(lmax.eta2) > 0.4),cex=0.6)
      }
    }
    Curv <- NULL
    for(i in 1:length(lmax)){
      li <- rep(0,length(lmax)); li[i] <- 1
      Curv <- rbind(Curv, 2*abs( t(li) %*% t(Delta) %*% solve(L) %*% Delta %*% li) )
    }
    cutoff1 <- sum(Curv)/(2*np)

    Curv.eta1 <- NULL
    for(i in 1:length(lmax)){
      li <- rep(0,length(lmax)); li[i] <- 1
      Curv.eta1 <- rbind(Curv.eta1, 2*abs( t(li) %*% t(Delta) %*% (solve(L) - B1.eta1) %*% Delta %*% li) )
    }
    cutoff2 <- sum(Curv.eta1)/(2*(np-ncol(X2)))

    par(mar=c(4.5,4.5,1,1))

    Curv.eta2 <- NULL
    for(i in 1:length(lmax)){
      li <- rep(0,length(lmax)); li[i] <- 1
      Curv.eta2 <- rbind(Curv.eta2, 2*abs( t(li) %*% t(Delta) %*% (solve(L) - B1.eta2) %*% Delta %*% li) )
    }
    cutoff3 <- sum(Curv.eta2)/(2*(np - ncol(X1)))

    curvmaxes <- max(c(Curv,Curv.eta1,Curv.eta2))
    cutmaxes  <- max(c(cutoff1,cutoff2,cutoff3))
    curvlim   <- max(c( (curvmaxes + curvmaxes/5), (cutmaxes+0.1)))

    if(curvs.plot){
      par(mar=c(4.5,4.5,1,1))
      if( (nuisance.plot==0 & !is.logical(nuisance.plot)) || is.na(nuisance.plot) ){
        plot(1:nrow(X1),Curv,pch=19,cex=0.7,ylim=c(0,curvlim),
             cex.axis=0.7,tcl=-0.3,
             xlab=expression(i),
             ylab=expression(paste(C[i],"(",theta[1],',',theta[2],")")))
        abline(h=cutoff1,lty=2)
      }

      if( is.na(nuisance.plot) || (!nuisance.plot & is.logical(nuisance.plot)) || (nuisance.plot==1 & !is.logical(nuisance.plot))){
      plot(1:nrow(X1),Curv.eta1,pch=19,cex=0.7,ylim=c(0,curvlim),
           cex.axis=0.7,tcl=-0.3,
           xlab=expression(i),
           ylab=expression(paste(C[i],"(",theta[1],")")))
      abline(h=cutoff2,lty=2)
      }

      if( is.na(nuisance.plot) || (nuisance.plot & is.logical(nuisance.plot)) || (nuisance.plot==2 & !is.logical(nuisance.plot))){
        par(mar=c(4.5,4.5,1,1))
        plot(1:nrow(X1),Curv.eta2,pch=19,cex=0.7,ylim=c(0,curvlim),
             cex.axis=0.7,tcl=-0.3,
             xlab=expression(i),
             ylab=expression(paste(C[i],"(",theta[2],")")))
        abline(h=cutoff3,lty=2)
      }

    }

    return(invisible(
      list(Delta=t(Delta),L=L,B1=list(B.theta1=BB1,B.theta2=BB2),Ci=data.frame(Ci.all=Curv,Ci.theta1=Curv.eta1,Curv.theta2=Curv.eta2)
      ))
    )
  }
}
