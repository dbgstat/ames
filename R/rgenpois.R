#' @rdname ames-deprecated
#' @export

rgenpois<-function(n,mu,phi){
  .Deprecated("rgp1",package="ames")
  return(rgp1(n,mu,phi))
}
