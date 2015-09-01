#' (internal)
#'
#' @export

print.fitames <- function(obj){
  if(class(obj)=="fitames"){
    print(obj$model)
    cat("\n\n")
    print(data.frame(Param=obj$coef$Param,Estimate=round(obj$coef$Estimate,4),Sd=round(obj$coef$Sd,4)))
    cat("\n\n")
    print(obj$convergence)
    if(is.null(obj$profiles)==F){
      cat("\n\n")
      print(obj$profiles)
    }
  }else stop("Requires a fitames object.")
}
