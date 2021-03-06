% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/dtest.ames.R
\name{dtest.ames}
\alias{dtest.ames}
\title{Dean score test for overdispersion}
\usage{
dtest.ames(data, count, dose, alternative.var = c("gp1", "nb2"),
  alternative = c("greater", "less", "two.sided"), predictor = c("krewski",
  "myers", "bernstein", "breslow", "stead", "margolin", "svetliza"),
  link = log, theta.start = "smart", controls = control.fitames(),
  quiet = F)
}
\arguments{
\item{...}{}
}
\description{
Dean score test for overdispersion
}
\examples{
iv.ames()
}
\author{
Davi Butturi-Gomes

Silvio S. Zocchi
}

