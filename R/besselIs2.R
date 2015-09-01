#' Computes the modified Bessel function of first kind
#'
#' Computes the modified Bessel I function, using one of its basic definitions as an infinite series.
#' This function is a simple modification of \code{besselIs} \{Bessel\} to allow for vectors.
#'
#' @param ...
#'
#' @author Davi Butturi-Gomes
#'
#' Silvio S. Zocchi
#'
#' @seealso \code{\link[Bessel]{besselIs}}
#'
#' @export

besselIs2<-function (x, nu, nterm = 800, expon.scaled = FALSE, log = FALSE,
                    Ceps = if (isNum) 8e-16 else 2^(-x@.Data[[1]]@prec))
{
  if(length(nu) > 1 && length(nu)==length(x) ){
    a <- c()
    for(i in 1:length(nu)){
      a <- c(a,besselIs(x=x[i],nu=nu[i],nterm=nterm,expon.scaled=expon.scaled,log=log) )
    }
    a
  }else{
  j <- (nterm - 1):0
  n <- length(x)
  if (n == 0)
    return(x)
  if (is(nu, "mpfr"))
    x <- mpfr(x, precBits = max(64, .getPrec(nu)))
  l.s.j <- outer(j, (x/2), function(X, Y) X * 2 * log(Y))
  isNum <- is.numeric(x) || is.complex(x)
  if (is(l.s.j, "mpfr"))
    j <- mpfr(j, precBits = max(sapply(l.s.j@.Data, slot,
                                       "prec")))
  else if (!isNum)
    j <- as(j, class(x))
  log.s.j <- if (expon.scaled)
    l.s.j - rep(x, each = nterm) - lgamma(j + 1) - lgamma(nu +
                                                            1 + j)
  else l.s.j - lgamma(j + 1) - lgamma(nu + 1 + j)
  s.j <- if (log)
    lsum(log.s.j)
  else exp(log.s.j)
  if (log) {
    if (any(i0 <- x == 0))
      s.j[i0] <- 0
    if (any(lrgS <- log.s.j[1, ] > log(Ceps) + s.j))
      lapply(x[lrgS], function(x) warning(sprintf("bI(x=%g): 'nterm' may be too small",
                                                  x), call. = FALSE))
    nu * log(x/2) + s.j
  }
  else {
    s <- colSums(s.j)
    if (any(i0 <- x == 0))
      s[i0] <- 0
    if (!all(iFin <- is.finite(s)))
      stop(sprintf("infinite s for x=%g", x[!iFin][1]))
    if (any(lrgS <- s.j[1, ] > Ceps * s))
      lapply(x[lrgS], function(x) warning(sprintf("bI(x=%g): 'nterm' may be too small",
                                                  x), call. = FALSE))
    (x/2)^nu * s
  }
  }
}
environment(besselIs2) <- environment(besselIs)
