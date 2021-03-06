% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/besselIs2.R
\name{besselIs2}
\alias{besselIs2}
\title{Computes the modified Bessel function of first kind}
\usage{
besselIs2(x, nu, nterm = 800, expon.scaled = FALSE, log = FALSE,
  Ceps = if (isNum) 8e-16 else 2^(-x@.Data[[1]]@prec))
}
\arguments{
\item{...}{}
}
\description{
Computes the modified Bessel I function, using one of its basic definitions as an infinite series.
This function is a simple modification of \code{besselIs} \{Bessel\} to allow for vectors.
}
\author{
Davi Butturi-Gomes

Silvio S. Zocchi
}
\seealso{
\code{\link[Bessel]{besselIs}}
}

