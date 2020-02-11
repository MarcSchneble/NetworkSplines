# adjust path
setwd("X:/Projekte/Network Splines/R Code/Git")

# load required packages
library(spatstat)
library(igraph)
library(Matrix)
library(MASS)
library(VCA) #Moore Penrose inverse for sparse matrices
library(colorRamps) # color palettes for plotting
library(splines) 

# load functions
source(file = "Functions.R")



######################################################################################################################################
######################################## Create Plots in Figure 3 ####################################################################
######################################################################################################################################


# specify the geometric network and augment it by several attributes that we need for out estimation routine
delta = 0.05
h = 0.01
r = 1
L = augment.linnet(simplenet, delta, h)

# plots for different n
n = c(5, 20, 100, 1000)
for (i in 1:length(n)) {
  set.seed(7)
  L.lpp = runiflpp(n = n[i], L)
  
  # intensity with penalized splines
  intens.psplines = intensity.psplines.lpp(L.lpp, r, rho = 10)
  
  # intensity with kernel based method
  sigma = bw.lppl(L.lpp)
  intens.kernel = density.lpp(L.lpp, sigma = sigma)
  
  # common color scale for both plots
  min.intens = min(min(intens.psplines$v[-which(is.na(intens.psplines$v))]), min(intens.kernel$v[-which(is.na(intens.kernel$v))]))
  max.intens = max(max(intens.psplines$v[-which(is.na(intens.psplines$v))]), max(intens.kernel$v[-which(is.na(intens.kernel$v))]))
  co.map = colourmap(matlab.like2(1000), range=c(min.intens, max.intens))
  
  # plots
  pdf(file = paste("Plots/SimpleNetwork/SimpleNetworkn", n[i], ".pdf", sep = ""), width = 10, height = 8)
  par(mar=c(0, 0, 2, 0), cex = 2.5)
  plot(L.lpp, main = paste("Simple Network with \n n = ", n[i], " Random Points", sep = ""), lwd = 3,
       cols = "red", chars = 16, size = 0.6)
  dev.off()
  
  pdf(file = paste("Plots/SimpleNetwork/IntensityPSplinesn", n[i], ".pdf", sep = ""), width = 10, height = 8)
  par(mar=c(0, 0, 1, 0.5), cex = 2.5)
  plot.linim(intens.psplines, main = paste("Kernel Estimate n = ", n[i], sep = ""), col = co.map)
  dev.off()
  
  pdf(file = paste("Plots/SimpleNetwork/IntensityKerneln", n[i], ".pdf", sep = ""), width = 10, height = 8)
  par(mar=c(0, 0, 1, 0.5), cex = 2.5)
  plot.linim(intens.kernel, main = paste("Penalized Splines Estimate n = ", n[i], sep = ""), col = co.map)
  dev.off()
}



######################################################################################################################################
######################################## Simulation Study with varying n #############################################################
######################################################################################################################################


# specify the geometric network and augment it by several attributes that we need for out estimation routine
delta = 0.05
h = 0.01
r = 1
L = augment.linnet(simplenet, delta, h)

# setting for simulation
S = 100
n = c(5, 10, 20, 50, 100, 200, 500, 1000)
ISE.n = matrix(0, S, length(n))
sd.n = rep(0, length(n))
set.seed(1)
for (i in 1:length(n)) {
  for (s in 1:S) {
    message(paste("Simulation number ", s, " for n = ", n[i], ".", sep = ""))
    L.lpp = runiflpp(n = n[i], L)
    ISE.n[s, i] = sum(ise.intensity.uniform(L.lpp, r))/n[i]^2
  }
  sd.n[i] = sd(ISE.n[, i])
  save(ISE.n, file = "ISEn.RData")
}


pdf(file = "Plots/MeanSampleSize.pdf", width = 5, height = 5)
par(mar=c(4.1, 4.5, 1.5, 0.5))
plot(log(sample.size), log(colMeans(ISE.n)), xlab = "log(n)", ylab = expression(log(bar(ISE)[n])), col = "red", pch = 16,
     main = "Mean of ISE for Different Sample Sizes n")
dev.off()

pdf(file = "Plots/SDSampleSize.pdf", width = 5, height = 5)
par(mar=c(4.1, 4.5, 1.5, 0.5))
plot(log(sample.size), log(sd.sample.size), xlab = "log(n)", ylab = expression(log(hat(sigma)(ISE[n]))), col = "red", pch = 16,
     main = "SD of ISE for Different Sample Sizes n")
dev.off()



######################################################################################################################################
######################################## Simulation Study with varying delta and h ###################################################
######################################################################################################################################


# setting for simulation
S = 100
n = 100
r = 1
delta = c(0.1, 0.05, 0.01)
h = c(0.1, 0.05, 0.01, 0.005) 
ISE.delta.h = array(0, c(S, length(delta), length(h)))
mean.ISE.delta.h = sd.ISE.delta.h = matrix(0, length(delta), length(h))
set.seed(1)
for (i in 1:length(delta)) {
  for (j in 1:length(h)) {
    if (h[j] <= delta[i]){
      for (s in 1:S) {
        message(paste("Simulation number ", s, " for delta = ", delta[i], " and h = ", h[j], ".", sep = ""))
        L = augment.linnet(simplenet, delta[i], h[j])
        L.lpp = runiflpp(n = n, L)
        ISE.delta.h[s, i, j] = sum(ise.intensity.uniform(L.lpp, r))/n^2
      }
    }
    mean.ISE.delta.h[i, j] = mean(ISE.delta.h[, i, j])
    sd.ISE.delta.h[i, j] = sd(ISE.delta.h[, i, j])
    save(mean.ISE.delta.h, file = "meanISEdeltah.RData")
    save(sd.ISE.delta.h, file = "meanISEdeltah.RData")
  }
}