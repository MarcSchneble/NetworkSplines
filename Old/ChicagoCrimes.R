# adjust path
setwd("~/NetworkSplines")

# load required packages
library(spatstat)
library(igraph) 
library(Matrix)
library(MASS)
library(VCA) 
library(splines) 

# load functions
source(file = "Old/FunctionsOld.R")

# augmented Chicago network
delta = 5
h = 1
r = 1
L.lpp = unmark(chicago)
X = as.ppp(L.lpp)
L = augment.linnet(as.linnet(L.lpp), delta, h)
L.lpp = lpp(X, L)

# intensity estimate with penalized spline method
intens.pspline = intensity.pspline.lpp(L.lpp, r)

# intensity estimate with kernel based approach
sigma = bw.lppl(L.lpp)
intens.kernel = density.lpp(L.lpp, sigma = sigma)

# common color scale for both plots
min.intens = min(min(intens.pspline$v[-which(is.na(intens.pspline$v))]), min(intens.kernel$v[-which(is.na(intens.kernel$v))]))
max.intens = max(max(intens.pspline$v[-which(is.na(intens.pspline$v))]), max(intens.kernel$v[-which(is.na(intens.kernel$v))]))

# plots
pdf(file = "Plots/ChicagoNetwork.pdf", width = 10, height = 8)
par(mar=c(0, 0, 2, 0), cex = 2.4)
plot(L.lpp, main = "Chicago Crime Network", lwd = 3,
     cols = "red", chars = 16, size = 1)
dev.off()

pdf(file = "Plots/ChicagoIntensityPSpline.pdf", width = 10, height = 8)
par(mar=c(0, 0, 3, 0.5), cex = 2.4)
plot.linim(intens.pspline, main = "Penalized Spline Estimate" , zlim = c(min.intens, max.intens), box = TRUE)
dev.off()

pdf(file = "Plots/ChicagoIntensityKernel.pdf", width = 10, height = 8)
par(mar=c(0, 0, 3, 0.5), cex = 2.4)
plot.linim(intens.kernel, main = "Kernel Based Estimate", zlim = c(min.intens, max.intens), box = TRUE) 
dev.off()