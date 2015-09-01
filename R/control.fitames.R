#' Controls of fitting procedures in fit.ames
#'
#'
#' @param ...
#'
#' @author Davi Butturi-Gomes
#'
#' Silvio S. Zocchi
#'
#' @examples
#' fit.ames()
#'
#' @export

control.fitames <- function(convergence.tol,
                            max.it,
                            zerocounts,
                            strain,
                            cfu,
                            alpha,
                            delta,
                            allow.jitter,
                            max.jitter,
                            dec,
                            tau.profile,
                            conv.crit,
                            halfstepping,
                            stepsize,
                            diagW.min,
                            detas.method,
                            nsim
                            ){

  if(missing(convergence.tol)) convergence.tol <- 1e-6
  if(missing(max.it)) max.it <- 30
  if(missing(zerocounts)) zerocounts <- 1e-3
  if(missing(strain)) strain <- "ta100"
  if(missing(cfu)) cfu <- 1e8
  if(missing(alpha)) alpha <- 2
  if(missing(delta)) delta <- 1
  if(missing(allow.jitter)) allow.jitter <- T
  if(missing(max.jitter)) max.jitter <- 20
  if(missing(dec)) dec <- 4
  if(missing(tau.profile)) tau.profile <- seq(0,4,length.out=1000)
  if(missing(conv.crit)) conv.crit <- "likelh"
  if(missing(halfstepping)) halfstepping <- T
  if(missing(stepsize)) stepsize <- 2
  if(missing(diagW.min)) diagW.min <- .Machine$double.eps^0.75
  if(missing(detas.method)) detas.method <- 1
  if(missing(nsim))  nsim <- 1000

  return(invisible(list(convergence.tol=convergence.tol,
       max.it=max.it,
       zerocounts=zerocounts,
       strain=strain,
       cfu=cfu,
       alpha=alpha,
       delta=delta,
       allow.jitter=allow.jitter,
       max.jitter=max.jitter,
       dec=dec,
       tau.profile=tau.profile,
       conv.crit=conv.crit,
       halfstepping=halfstepping,
       stepsize=stepsize,
       diagW.min=diagW.min,
       detas.method=detas.method,
       nsim=nsim
       )
  ))
}
