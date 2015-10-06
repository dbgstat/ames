#' @rdname ames-deprecated
#' @export

dgenpois <- function(x,mu,phi){
  .Deprecated("dgp1",package="ames")
  return( dgp1(x,mu,phi,log=F)  )
}
