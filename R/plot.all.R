#' (internal)
#'
#' @export

plot.all <- function(obj,nsim=1000){
  if(obj$fit.method=="gnlm"){
    plot.fitames(obj,
                 smooth=T,
                 deviance.plots=T,
                 std.deviance.plots=T,
                 pearson.plots=T,
                 std.pearson.plots=T,
                 qqplot=T,
                 lag.plot=T,lag=1,
                 diag.measures=T,
                 envelope=T,nsim=nsim,probs=c(0.005,0.995),
                 quiet=T
    )
  }else{
    plot.fitames(obj,
                 smooth=T,
                 deviance.plots=F,
                 std.deviance.plots=F,
                 pearson.plots=T,
                 std.pearson.plots=T,
                 qqplot=T,
                 lag.plot=T,lag=1,
                 diag.measures=T,
                 envelope=T,nsim=nsim,probs=c(0.005,0.995),
                 quiet=T
    )
  }
}
