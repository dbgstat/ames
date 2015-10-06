rm(list=ls())
require(devtools)
require(roxygen2)
# require(distr)


######
install_github('dbgstat/ames')
require(ames)
?rgp1
citation("ames")


#####
document()
setwd("..")
install("ames")
setwd("./ames")

x <- rgenpois(1000,mu=mu<-20,phi=phi<-2)
?dgp1
?dgenpois


x <- rgp1(1000,mu=mu<-20,phi=phi<-2)
mu; mu*phi # Nominal mean and variance
mean(x);var(x) # "Observed" mean and variance
#
dx <- (dgp1(15:25,mu=mu,phi=phi))
plot(15:25,dx,type='h')
