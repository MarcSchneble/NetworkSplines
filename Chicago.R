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

# load functions
source(file = "Functions.R")

# augmented Chicago network
delta <- 10
h <- 2
r <- 2
L.lpp <- unmark(chicago)
X <- as.ppp(L.lpp)
L <- augment.linnet(as.linnet(L.lpp), delta, h)
L.lpp <- lpp(X, L)
K.ginv.net <- ginv(as.matrix(Matrix(K.linnet(L, r))))

# intensity estimate with penalized spline method
intens.pspline <- intensity.pspline.lpp(L.lpp, r, K.ginv.net = K.ginv.net)

# intensity estimate with kernel based approach
sigma <- as.numeric(bw.lppl(L.lpp))
intens.kernel <- density.lpp(L.lpp, sigma = sigma)

# intensity with Voronoi estimate
f <- as.numeric(bw.voronoi(chicago))
intens.voronoi <- densityVoronoi(chicago, f = f, nrep = 100)

# simulate
varphi <- as.linfun(intens.pspline)
a <- simulation.chicago(1, delta, h, r, 200, varphi, kernel = TRUE, K.ginv.net = K.ginv.net)
a

