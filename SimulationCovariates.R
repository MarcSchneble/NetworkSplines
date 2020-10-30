# adjust path
setwd("~/NetworkSplines")
rm(list = ls())

# load required packages
library(spatstat)
library(igraph)
library(Matrix)
library(MASS)
library(VCA) 
library(splines) 
library(parallel) 
library(dplyr)
library(tidyr)

# load functions
source(file = "Functions.R")

delta = 0.05
h = 0.01
r = 1
L = augment.linnet(simplenet, delta, h)
beta <- c(1, 2)

t <- runif(1)
mu <- exp(beta[1] + t*beta[2])
n <- rpois(1, mu)  
covariates <- rep(t, n)
L.lpp <- runiflpp(n = n, L)

for (i in 2:100) {
  t <- rnorm(1)
  mu <- exp(beta[1] + t*beta[2])
  n <- rpois(1, mu)  
  covariates <- c(covariates, rep(t, n))
  L.lpp <- superimpose.lpp(L.lpp, runiflpp(n = n, L))
}


L.lpp$data$t <- round(covariates, 2)

lins <- NULL
smooths <- NULL
fit <- fit.lpp(L.lpp, r, smooths = smooths, lins = lins)
print(fit)

design <- get.design(L.lpp, r, smooths = smooths, lins = lins)
y <- exp(as.vector(design$Z[, design$ind.smooths[[1]]]%*%fit$theta.hat[design$ind.smooths[[1]]]))

#x <- design$data$t
#y <- as.vector(design$Z[, design$ind.smooths[[2]]]%*%fit$theta.hat[design$ind.smooths[[2]]]) 
#plot(x, y) 

#x <- design$data$dist2V 
#y <- as.vector(design$Z[, design$ind.smooths[[3]]]%*%fit$theta.hat[design$ind.smooths[[3]]]) 
#plot(x, y