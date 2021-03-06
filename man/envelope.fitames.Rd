% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/envelope.fitames.R
\name{envelope.fitames}
\alias{envelope.fitames}
\title{Simulated envelope in normal plots for fitames objects}
\usage{
envelope.fitames(obj, nsim, ris.plot = F, rstar.plot = F, dis.plot = F,
  dstar.plot = F, nuisance.plot = F, eta.plot = NULL)
}
\arguments{
\item{obj}{A model fitted by \code{\link{fit.ames}}}

\item{nsim}{Integer, specifying the number of simulations}

\item{ris.plot}{See \code{pearson.plots} of \code{\link{plot.fitames}}}

\item{rstar.plot}{See \code{std.pearson.plots} of \code{\link{plot.fitames}}}

\item{dis.plot}{See \code{deviance.plots} of \code{\link{plot.fitames}}}

\item{dstar.plot}{See \code{std.deviance.plots} of \code{\link{plot.fitames}}}

\item{nuisance.plot}{See \code{\link{diagnostics.fitames}}}

\item{eta.plot}{See \code{\link{diagnostics.fitames}}}
}
\description{
Used for computing simulated envelopes in ordered residual plots \emph{versus}
normal quantiles (normal plots). Although internally called within \code{\link{plot.fitames}}
given the right arguments, may be used in an independent fashion.
}
\examples{
fit.ames()
}
\author{
Davi Butturi-Gomes

Silvio S. Zocchi
}
\seealso{
\code{\link{diagnostics.fitames}} \code{\link{plot.fitames}}
}

