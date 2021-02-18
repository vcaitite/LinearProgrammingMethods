rm(list = ls())
source("~/Documents/UFMG/9/Otimização/Codes_R/InteriorPointsFunction.R")

A <- cbind(c(1, 0, 3, -1, 0),c(0,1,2,0,-1))
print(A)
c <- matrix(c(3,5), ncol = 1)
print(c)
b <- matrix(c(4,6,18,0,0), ncol = 1)
print(b)
x0 <- matrix(c(1,2), ncol = 1) # É um ponto viável
gama <- 0.5
maxIter <- 100
epslon <- 10^-3
resultList <- InteriorPoint(A, b, c, gama, x0, maxIter, epslon, plot=TRUE)