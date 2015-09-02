#' Skellam distribution
#'
#' Density of the Skellam distribution, slightly modified from code{dskellam} \{VGAM\}.
#'
#' @param x vector of quantiles
#' @param mu1 First parameter for the mean of the distribution
#' @param mu2 Second parameter for the mean of the distribution
#' @param log Logical. If TRUE, returns the density in the natural logarithm scale.
#' Defaults to FALSE.
#'
#' @author Davi Butturi-Gomes
#'
#' Silvio S. Zocchi
#'
#' @seealso \code{\link[VGAM]{dskellam}}
#'
#' @export

dskellam2 <- function (x, mu1, mu2, log = FALSE)
{
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)
  L <- max(length(x), length(mu1), length(mu2))
  if (length(x) != L)
    x <- rep(x, len = L)
  if (length(mu1) != L)
    mu1 <- rep(mu1, len = L)
  if (length(mu2) != L)
    mu2 <- rep(mu2, len = L)
  ok2 <- is.finite(mu1) & is.finite(mu2) & (mu1 >= 0) & (mu2 >=
                                                           0)
  ok3 <- (mu1 == 0) & (mu2 > 0)
  ok4 <- (mu1 > 0) & (mu2 == 0)
  ok5 <- (mu1 == 0) & (mu2 == 0)
  if (log.arg) {
    ans <- -mu1 - mu2 + 2 * sqrt(mu1 * mu2) + 0.5 * x * log(mu1) -
      0.5 * x * log(mu2) + besselIs2(2 * sqrt(mu1 * mu2),
                                       nu = abs(x), expon.scaled = TRUE,log=T)
    ans[ok3] <- dpois(x = -x[ok3], lambda = mu2[ok3], log = TRUE)
    ans[ok4] <- dpois(x = -x[ok4], lambda = mu1[ok4], log = TRUE)
    ans[ok5] <- dpois(x = x[ok5], lambda = 0, log = TRUE)
    ans[x != round(x)] = log(0)
  }
  else {
    ans <- (mu1/mu2)^(x/2) * exp(-mu1 - mu2 + 2 * sqrt(mu1 *
                                                         mu2)) * exp( besselIs2(2 * sqrt(mu1 * mu2), nu = abs(x),
                                                                         expon.scaled = TRUE,log=T))
    ans[ok3] <- dpois(x = -x[ok3], lambda = mu2[ok3])
    ans[ok4] <- dpois(x = -x[ok4], lambda = mu1[ok4])
    ans[ok5] <- dpois(x = x[ok5], lambda = 0)
    ans[x != round(x)] <- 0
  }
  #ans[!ok2] <- NaN
  ans
}
environment(dskellam2) <- environment(VGAM::dskellam)
