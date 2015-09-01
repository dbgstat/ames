#' Skellam Distribution Family Function
#'
#' Estimates the two parameters of a Skellam distribution by maximum likelihood estimation.
#' This is to be used along with \code{vglm} function, available within \{VGAM\} package.
#'
#'
#' @author Davi Butturi-Gomes
#'
#' Silvio S. Zocchi
#'
#' @seealso \code{\link[VGAM]{skellam}}, \code{\link[VGAM]{vglm}}
#'
#' @export

skellam2 <- function (lmu1 = "loge", lmu2 = "loge", imu1 = NULL, imu2 = NULL,
          nsimEIM = 100, parallel = FALSE, zero = NULL,...)
{
  lmu1 <- as.list(substitute(lmu1))
  emu1 <- link2list(lmu1)
  lmu1 <- attr(emu1, "function.name")
  lmu2 <- as.list(substitute(lmu2))
  emu2 <- link2list(lmu2)
  lmu2 <- attr(emu2, "function.name")
  if (length(imu1) && !is.Numeric(imu1, positive = TRUE))
    stop("bad input for argument 'imu1'")
  if (length(imu2) && !is.Numeric(imu2, positive = TRUE))
    stop("bad input for argument 'imu2'")
  if (!is.Numeric(nsimEIM, length.arg = 1, integer.valued = TRUE) ||
        nsimEIM <= 50)
    stop("argument 'nsimEIM' should be an integer greater than 50")
  new("vglmff", blurb = c("Skellam distribution\n\n", "Links:    ",
                          namesof("mu1", lmu1, earg = emu1, tag = FALSE), ", ",
                          namesof("mu2", lmu2, earg = emu2, tag = FALSE), "\n",
                          "Mean:     mu1-mu2", "\n", "Variance: mu1+mu2"), constraints = eval(substitute(expression({
                            constraints <- cm.VGAM(matrix(1, M, 1), x = x, bool = .parallel,
                                                   constraints = constraints, apply.int = TRUE)
                            constraints <- cm.zero.VGAM(constraints, x, .zero, M)
                          }), list(.parallel = parallel, .zero = zero))), initialize = eval(substitute(expression({
                            temp5 <- w.y.check(w = w, y = y, ncol.w.max = 1, ncol.y.max = 1,
                                               Is.integer.y = TRUE, out.wy = TRUE, maximize = TRUE)
                            w <- temp5$w
                            y <- temp5$y
                            predictors.names <- c(namesof("mu1", .lmu1, earg = .emu1,
                                                          tag = FALSE), namesof("mu2", .lmu2, earg = .emu2,
                                                                                tag = FALSE))
                            if (!length(etastart)) {
                              junk <- lm.wfit(x = x, y = c(y), w = c(w))
                              var.y.est <- sum(c(w) * junk$resid^2)/junk$df.residual
                              mean.init <- weighted.mean(y, w)
                              mu1.init <- max((var.y.est + mean.init)/2, 0.01)
                              mu2.init <- max((var.y.est - mean.init)/2, 0.01)
                              mu1.init <- rep(if (length(.imu1)) .imu1 else mu1.init,
                                              length = n)
                              mu2.init <- rep(if (length(.imu2)) .imu2 else mu2.init,
                                              length = n)
                              etastart <- cbind(theta2eta(mu1.init, .lmu1, earg = .emu1),
                                                theta2eta(mu2.init, .lmu2, earg = .emu2))
                            }
                          }), list(.lmu1 = lmu1, .lmu2 = lmu2, .imu1 = imu1, .imu2 = imu2,
                                   .emu1 = emu1, .emu2 = emu2))), linkinv = eval(substitute(function(eta,
                                                                                                     extra = NULL) {
                                     mu1 <- eta2theta(eta[, 1], link = .lmu1, earg = .emu1)
                                     mu2 <- eta2theta(eta[, 2], link = .lmu2, earg = .emu2)
                                     mu1 - mu2
                                   }, list(.lmu1 = lmu1, .lmu2 = lmu2, .emu1 = emu1, .emu2 = emu2))),
      last = eval(substitute(expression({
        misc$link <- c(mu1 = .lmu1, mu2 = .lmu2)
        misc$earg <- list(mu1 = .emu1, mu2 = .emu2)
        misc$expected <- TRUE
        misc$nsimEIM <- .nsimEIM
      }), list(.lmu1 = lmu1, .lmu2 = lmu2, .emu1 = emu1, .emu2 = emu2,
               .nsimEIM = nsimEIM))), loglikelihood = eval(substitute(function(mu,
                                                                               y, w, residuals = FALSE, eta, extra = NULL, summation = TRUE) {
                 mu1 <- eta2theta(eta[, 1], link = .lmu1, earg = .emu1)
                 mu2 <- eta2theta(eta[, 2], link = .lmu2, earg = .emu2)
                 if (residuals) {
                   stop("loglikelihood residuals not implemented yet")
                 } else {
                   ll.elts <- if (is.logical(.parallel) && length(.parallel) ==
                                    1 && .parallel) c(w) * besselIs2(2 * mu1,
                                                                       nu = abs(y), expon = TRUE, log=T, ...) else c(w) * (-mu1 -
                                                                                                             mu2 + 0.5 * y * log(mu1) - 0.5 * y * log(mu2) +
                                                                                                             2 * sqrt(mu1 * mu2) + besselIs2(2 * sqrt(mu1 * mu2), nu = abs(y), expon = TRUE,log=T,...))
                   if (summation) {
                     sum(ll.elts)
                   } else {
                     ll.elts
                   }
                 }
               }, list(.lmu1 = lmu1, .lmu2 = lmu2, .emu1 = emu1, .emu2 = emu2,
                       .parallel = parallel))), vfamily = c("skellam"),
      simslot = eval(substitute(function(object, nsim) {
        pwts <- if (length(pwts <- object@prior.weights) >
                      0) pwts else weights(object, type = "prior")
        if (any(pwts != 1)) warning("ignoring prior weights")
        eta <- predict(object)
        mu1 <- eta2theta(eta[, 1], link = .lmu1, earg = .emu1)
        mu2 <- eta2theta(eta[, 2], link = .lmu2, earg = .emu2)
        rskellam(nsim * length(mu1), mu1, mu2)
      }, list(.lmu1 = lmu1, .lmu2 = lmu2, .emu1 = emu1, .emu2 = emu2,
              .parallel = parallel))), deriv = eval(substitute(expression({
                mu1 <- eta2theta(eta[, 1], link = .lmu1, earg = .emu1)
                mu2 <- eta2theta(eta[, 2], link = .lmu2, earg = .emu2)
                dmu1.deta <- dtheta.deta(mu1, link = .lmu1, earg = .emu1)
                dmu2.deta <- dtheta.deta(mu2, link = .lmu2, earg = .emu2)
                temp8 <- 2 * sqrt(mu1 * mu2)
                temp9 <- besselIs2(temp8, nu = abs(y), expon = TRUE,log=T,...)
                temp7 <- besselIs2(temp8, nu = abs(y - 1), expon = TRUE,log=T,...) - log(2)
                temp10 <- besselIs2(temp8, nu = abs(y + 1), expon = TRUE,log=T,...) - log(2)

                temp6 <- exp(temp7-temp9)+exp(temp10-temp9)

                #print(cbind(temp6,temp7,temp9,temp10,temp8))

                dl.dmu1 <- -1 + 0.5 * y/mu1 + sqrt(mu2/mu1) * temp6
                dl.dmu2 <- -1 - 0.5 * y/mu2 + sqrt(mu1/mu2) * temp6
                c(w) * cbind(dl.dmu1 * dmu1.deta, dl.dmu2 * dmu2.deta)
              }), list(.lmu1 = lmu1, .lmu2 = lmu2, .emu1 = emu1, .emu2 = emu2,
                       .nsimEIM = nsimEIM))), weight = eval(substitute(expression({
                         run.var <- run.cov <- 0
                         for (ii in 1:(.nsimEIM)) {
                           ysim <- rskellam(n, mu1 = mu1, mu2 = mu2)
                           temp9 <- besselIs2(temp8, nu = abs(ysim), expon = TRUE,log=T,...)
                           temp7 <- besselIs2(temp8, nu = abs(ysim - 1), expon = TRUE,log=T,...) - log(2)
                           temp10 <- besselIs2(temp8, nu = abs(ysim + 1), expon = TRUE,log=T,...) - log(2)

                           temp6 <- exp(temp7-temp9)+exp(temp10-temp9)
                           dl.dmu1 <- -1 + 0.5 * ysim/mu1 + sqrt(mu2/mu1) *
                             temp6
                           dl.dmu2 <- -1 - 0.5 * ysim/mu2 + sqrt(mu1/mu2) *
                             temp6
                           rm(ysim)
                           temp3 <- cbind(dl.dmu1, dl.dmu2)
                           run.var <- ((ii - 1) * run.var + temp3^2)/ii
                           run.cov <- ((ii - 1) * run.cov + temp3[, 1] *
                                         temp3[, 2])/ii
                         }
                         wz <- if (intercept.only) matrix(colMeans(cbind(run.var,
                                                                         run.cov)), n, dimm(M), byrow = TRUE) else cbind(run.var,
                                                                                                                         run.cov)
                         dtheta.detas <- cbind(dmu1.deta, dmu2.deta)
                         index0 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
                         wz <- wz * dtheta.detas[, index0$row] * dtheta.detas[,
                                                                              index0$col]
                         c(w) * wz
                       }), list(.lmu1 = lmu1, .lmu2 = lmu2, .emu1 = emu1, .emu2 = emu2,
                                .nsimEIM = nsimEIM))))
}

environment(skellam2) <- environment(skellam)
