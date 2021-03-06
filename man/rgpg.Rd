% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/gpg.R
\name{dgpg}
\alias{dgpg}
\alias{rgpg}
\title{Generalized Poisson Type-G distribution (GP-G)}
\usage{
dgpg(x, mu, phi, g, log = F, use.gp12 = T)

rgpg(n, mu, phi, g, use.gp12 = T)
}
\arguments{
\item{x}{Integer, vector of integers or numeric that can be coerced into integers.}

\item{mu}{Positive numeric. Mean parameter of the distribution}

\item{phi}{Positive numeric. Dispersion parameter.}

\item{g}{Positive numeric <= 1. Functional parameter of the distribution.}

\item{log}{Logical. If TRUE, return the probabilities in log scale. The default is FALSE.}

\item{use.gp12}{Logical. If TRUE (default) recovers the GP-1 parametrization for g=1 or,
for g=2, recovers the GP-2 parametrization. This is recommended due otimization purposes.}

\item{n}{Integer or numeric that can be coerced into integer. Number of random deviations
to be generated}
}
\description{
\code{dgpg} - Generalized Poisson type-G distribution probability mass function

\code{rgpg} - Generalized Poisson type-G distribution random number generator
}
\examples{
x <- rgpg(1000,mu=mu<-20,phi=phi<-0.01,g=g<-2) # Same as using rgp2.
mu; mu*((1+phi*(mu^(g-1)))^2) # nominal mean and variance
mean(x);var(x) # "observed" mean and variance
#
dx <- (dgpg(10:30,mu=mu,phi=phi,g=g))
plot(10:30,dx,type='h')
}
\author{
Davi Butturi-Gomes

Silvio S. Zocchi
}
\seealso{
\code{\link{dgp1}}, \code{\link{rgp1}}, \code{\link{dgp2}}, \code{\link{rgp2}}
}

