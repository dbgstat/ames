#' Log-likelihood on convergence of a fitames object
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

logLik.fitames <- function(obj){
  if(class(obj)!="fitames") stop("Requires a \'fit.ames()\' output object")
  return(obj$convergence$LogLikelihood )
}

