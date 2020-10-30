# adjust path
setwd("~/NetworkSplines")

# load required packages
library(spatstat)
library(igraph)
library(Matrix)
library(MASS)
library(VCA) 
library(splines) 
library(parallel) 

# load functions
source(file = "Old/FunctionsOld.R")

# choose kind of simulation study
create.plots = FALSE
uniform.n = FALSE
uniform.delta.h = TRUE
exp.n = TRUE
exp.delta.h = TRUE


##### Create Plots in Figure 3 ####


if (create.plots){
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
    intens.psplines = intensity.pspline.lpp(L.lpp, r, rho = 10)
    
    # intensity with kernel based method
    sigma = bw.lppl(L.lpp)
    intens.kernel = density(L.lpp, sigma = sigma)
    
    # common color scale for both plots
    min.intens = min(min(intens.psplines$v[-which(is.na(intens.psplines$v))]), 
                     min(intens.kernel$v[-which(is.na(intens.kernel$v))]))
    max.intens = max(max(intens.psplines$v[-which(is.na(intens.psplines$v))]), 
                     max(intens.kernel$v[-which(is.na(intens.kernel$v))]))

    # plots
    pdf(file = paste("Plots/SimpleNetwork/SimpleNetworkn", n[i], ".pdf", sep = ""), width = 10, height = 8)
    par(mar=c(0, 0, 2, 0), cex = 2.5)
    plot(L.lpp, main = paste("Simple Network with \n n = ", n[i], " Random Points", sep = ""), lwd = 3,
         cols = "red", chars = 16, size = 1)
    dev.off()
    
    pdf(file = paste("Plots/SimpleNetwork/IntensityPSplinen", n[i], ".pdf", sep = ""), width = 10, height = 8)
    par(mar=c(0, 0, 1, 0.5), cex = 2.5)
    plot.linim(intens.psplines, main = paste("Penalized Spline Estimate, n = ", n[i], sep = ""), zlim = c(min.intens, max.intens),
               box= TRUE)
    dev.off()
    
    pdf(file = paste("Plots/SimpleNetwork/IntensityKerneln", n[i], ".pdf", sep = ""), width = 10, height = 8)
    par(mar=c(0, 0, 1, 0.5), cex = 2.5)
    plot.linim(intens.kernel, main = paste("Kernel Based Estimate, n = ", n[i], sep = ""), zlim = c(min.intens, max.intens),
               box = TRUE)
    dev.off()
  }
}


#### Uniform intensity with varying n ####


if (uniform.n){
  # specify the geometric network and augment it by several attributes that we need for out estimation routine
  delta = 0.05
  h = 0.01
  r = 1
  L = augment.linnet(simplenet, delta, h)
  
  S = 1000
  n = c(5, 10, 20, 50, 100, 200, 500, 1000)
  ISE.uniform.n = data.frame(matrix(0, length(n), 5))
  colnames(ISE.uniform.n) = c("n", "Mean.PSpline", "SD.Spline", "Mean.Kernel", "SD.Kernel")
  ISE.uniform.n$n = n
  
  no_cores = detectCores() - 2
  set.seed(1)
  for (i in 1:length(n)) {
    message(paste("S = ", S, " simulations with n = ", n[i], ".", sep = ""))
    # specification of uniform intensity
    H = function(x, y, seg, tp) {n[i]/L$d.L}
    varphi = linfun(H, L)
    
    cl = makeCluster(no_cores)
    clusterExport(cl, c("delta", "h", "r", "n", "L", "varphi", "i",
                        "augment.linnet", "B.linnet", "gamma.hat.lpp",
                        "intensity.edge", "intensity.pspline.lpp",
                        "ise.intensity.splines", "ise.linfun", 
                        "IWLS", "K.linnet", "offset.bin", "simulation.study",
                        "squared.diff.linfun", "y.lpp"))
    
    clusterEvalQ(cl, library(spatstat)) 
    clusterEvalQ(cl, library(igraph))
    clusterEvalQ(cl, library(Matrix))   
    clusterEvalQ(cl, library(MASS))   
    clusterEvalQ(cl, library(VCA))   
    clusterEvalQ(cl, library(splines)) 
    
    ISE = parLapply(cl, 1:S, simulation.study, delta = delta, h = h, r = r, n = n[i], varphi = varphi, L = L)
    
    stopCluster(cl) 
    
    ISE.uniform.n$Mean.PSpline[i] = mean(unlist(ISE)[which(names(unlist(ISE)) == "ISE.pspline")])
    ISE.uniform.n$SD.Spline[i] = sd(unlist(ISE)[which(names(unlist(ISE)) == "ISE.pspline")])
    ISE.uniform.n$Mean.Kernel[i] = mean(unlist(ISE)[which(names(unlist(ISE)) == "ISE.kernel")])
    ISE.uniform.n$SD.Kernel[i] = sd(unlist(ISE)[which(names(unlist(ISE)) == "ISE.kernel")])
    
    save(ISE.uniform.n, file = "RData/ISEuniformn.RData")
  }
  
  pdf(file = "Plots/MeanUniformn.pdf", width = 10, height = 10)
  par(mar=c(4.1, 4.5, 1.5, 0.5), cex = 2.4)
  ylim = c(min(log(ISE.uniform.n$Mean.PSpline), log(ISE.uniform.n$Mean.Kernel)), 
           max(log(ISE.uniform.n$Mean.PSpline), log(ISE.uniform.n$Mean.Kernel))) 
  plot(log(n), log(ISE.uniform.n$Mean.PSpline), xlab = "log(n)", ylab = expression(log(bar(ISE)[n])), col = "red", pch = 16,
       main = expression(bold("Mean ISE")), type = "b", ylim = ylim)
  points(log(n), log(ISE.uniform.n$Mean.Kernel), col = "blue", pch = 17, type = "b")
  legend("bottomleft", legend = c("Penalized Spline Estimate", "Kernel Based Estimate"), col = c("red", "blue"), pch = c(16, 17))
  dev.off()
  
  pdf(file = "Plots/SDUniformn.pdf", width = 10, height = 10)
  par(mar=c(4.1, 4.5, 1.5, 0.5), cex = 2.4)
  ylim = c(min(log(ISE.uniform.n$SD.Spline), log(ISE.uniform.n$SD.Kernel)), 
           max(log(ISE.uniform.n$SD.Spline), log(ISE.uniform.n$SD.Kernel)))
  plot(log(n), log(ISE.uniform.n$SD.Spline), xlab = "log(n)", ylab = expression(log(hat(sigma)(ISE[n]))), col = "red", pch = 16,
       main = expression(bold("Standard Deviation of ISE")), type = "b", ylim = ylim)
  points(log(n), log(ISE.uniform.n$SD.Kernel), col = "blue", pch = 17, type = "b")
  legend("bottomleft", legend = c("Penalized Spline Estimate", "Kernel Based Estimate"), col = c("red", "blue"), pch = c(16, 17))
  dev.off()
}


#### Uniform intensity with varying delta and h ####


if (uniform.delta.h){
  # setting for simulation
  S = 1000
  n = 100
  r = 1
  delta = c(0.1, 0.05, 0.01)
  h = c(0.1, 0.05, 0.01, 0.005) 
  mean.uniform.delta.h = sd.uniform.delta.h = matrix(0, length(delta), length(h))
  
  no_cores = detectCores() - 2
  set.seed(1)
  for (i in 1:length(delta)) {
    for (j in 1:length(h)) {
      if (h[j] <= delta[i]){
        message(paste("S = ", S, " simulations with delta = ", delta[i], " and h = ", h[j], ".", sep = ""))
        # specification of uniform intensity
        L = augment.linnet(simplenet, delta[i], h[j])
        H = function(x, y, seg, tp) {n/L$d.L}
        varphi = linfun(H, L)
        
        cl = makeCluster(no_cores)
        clusterExport(cl, c("delta", "h", "r", "n", "L", "varphi", "i", "j",
                            "augment.linnet", "B.linnet", "gamma.hat.lpp",
                            "intensity.edge", "intensity.pspline.lpp",
                            "ise.intensity.splines", "ise.linfun", 
                            "IWLS", "K.linnet", "offset.bin", "simulation.study",
                            "squared.diff.linfun", "y.lpp"))
        
        clusterEvalQ(cl, library(spatstat)) 
        clusterEvalQ(cl, library(igraph))
        clusterEvalQ(cl, library(Matrix))   
        clusterEvalQ(cl, library(MASS))   
        clusterEvalQ(cl, library(VCA))   
        clusterEvalQ(cl, library(splines)) 
        
        ISE = parLapply(cl, 1:S, simulation.study, delta = delta[i], h = h[j], 
                        r = r, n = n, varphi = varphi, L = L, kernel = FALSE)
        
        stopCluster(cl) 
        
        mean.uniform.delta.h[i, j] = mean(unlist(ISE)[which(names(unlist(ISE)) == "ISE.pspline")])
        sd.uniform.delta.h[i, j] = sd(unlist(ISE)[which(names(unlist(ISE)) == "ISE.pspline")])
        save(mean.uniform.delta.h, file = "RData/meanISEuniformdeltah.RData")
        save(sd.uniform.delta.h, file = "RData/sdISEuniformdeltah.RData")
      }
    }
  }
}


#### Non-uniform intensity with varying n ####


if (exp.n){
  # specify the geometric network and augment it by several attributes that we need for out estimation routine
  delta = 0.05
  h = 0.01
  r = 1
  L = augment.linnet(simplenet, delta, h)
  
  S = 1000
  n = c(5, 10, 20, 50, 100, 200, 500, 1000)
  ISE.exp.n = data.frame(matrix(0, length(n), 5))
  colnames(ISE.exp.n) = c("n", "Mean.PSpline", "SD.Spline", "Mean.Kernel", "SD.Kernel")
  ISE.exp.n$n = n
  
  no_cores = detectCores() - 2
  set.seed(1)
  for (i in 1:length(n)) {
    message(paste("S = ", S, " simulations with n = ", n[i], ".", sep = ""))
    # specification of true intensity
    H = function(x, y, seg, tp) {sqrt(y)*exp(-x*y)*n[i]/integral.linfun(linfun(function(x, y, seg, tp) {sqrt(y)*exp(-x*y)}, L))}
    varphi = linfun(H, L)
    
    cl = makeCluster(no_cores)
    clusterExport(cl, c("delta", "h", "r", "n", "L", "varphi", "i",
                        "augment.linnet", "B.linnet", "gamma.hat.lpp",
                        "intensity.edge", "intensity.pspline.lpp",
                        "ise.intensity.splines", "ise.linfun", 
                        "IWLS", "K.linnet", "offset.bin", "simulation.study",
                        "squared.diff.linfun", "y.lpp"))
    
    clusterEvalQ(cl, library(spatstat)) 
    clusterEvalQ(cl, library(igraph))
    clusterEvalQ(cl, library(Matrix))   
    clusterEvalQ(cl, library(MASS))   
    clusterEvalQ(cl, library(VCA))   
    clusterEvalQ(cl, library(splines)) 
    
    ISE = parLapply(cl, 1:S, simulation.study, delta = delta, h = h, r = r, n = n[i], varphi = varphi, L = L)
    
    stopCluster(cl) 
    
    ISE.exp.n$Mean.PSpline[i] = mean(unlist(ISE)[which(names(unlist(ISE)) == "ISE.pspline")])
    ISE.exp.n$SD.Spline[i] = sd(unlist(ISE)[which(names(unlist(ISE)) == "ISE.pspline")])
    ISE.exp.n$Mean.Kernel[i] = mean(unlist(ISE)[which(names(unlist(ISE)) == "ISE.kernel")])
    ISE.exp.n$SD.Kernel[i] = sd(unlist(ISE)[which(names(unlist(ISE)) == "ISE.kernel")])
    
    save(ISE.exp.n, file = "RData/ISEexpn.RData")
  }
  
  pdf(file = "Plots/IntensityExp.pdf", width = 10, height = 8)
  par(mar=c(0, 0, 1.0, 1.5), cex = 2.4)
  H = function(x, y, seg, tp) {sqrt(y)*exp(-x*y)/integral.linfun(linfun(function(x, y, seg, tp) {sqrt(y)*exp(-x*y)}, L))}
  varphi = linfun(H, L)
  plot(varphi, box = TRUE, main ="")
  title("Non-Uniform Intensity (n=1)")
  dev.off()
  
  pdf(file = "Plots/MeanExpn.pdf", width = 10, height = 10)
  par(mar=c(4.1, 4.5, 1.5, 0.5), cex = 2.4)
  ylim = c(min(log(ISE.exp.n$Mean.PSpline), log(ISE.exp.n$Mean.Kernel)), 
           max(log(ISE.exp.n$Mean.PSpline), log(ISE.exp.n$Mean.Kernel))) 
  plot(log(n), log(ISE.exp.n$Mean.PSpline), xlab = "log(n)", ylab = expression(log(bar(ISE)[n])), col = "red", pch = 16,
       main = expression(bold("Mean  ISE")), type = "b", ylim = ylim)
  points(log(n), log(ISE.exp.n$Mean.Kernel), col = "blue", pch = 17, type = "b")
  legend("bottomleft", legend = c("Penalized Spline Estimate", "Kernel Based Estimate"), col = c("red", "blue"), pch = c(16, 17))
  dev.off()
  
  pdf(file = "Plots/SDExpn.pdf", width = 10, height = 10)
  par(mar=c(4.1, 4.5, 1.5, 0.5), cex = 2.4)
  ylim = c(min(log(ISE.exp.n$SD.Spline), log(ISE.exp.n$SD.Kernel)), 
           max(log(ISE.exp.n$SD.Spline), log(ISE.exp.n$SD.Kernel)))
  plot(log(n), log(ISE.exp.n$SD.Spline), xlab = "log(n)", ylab = expression(log(hat(sigma)(ISE[n]))), col = "red", pch = 16,
       main = expression(bold("Standard Deviation of ISE")), type = "b", ylim = ylim)
  points(log(n), log(ISE.exp.n$SD.Kernel), col = "blue", pch = 17, type = "b")
  legend("bottomleft", legend = c("Penalized Spline Estimate", "Kernel Based Estimate"), col = c("red", "blue"), pch = c(16, 17))
  dev.off()
}


#### Non-uniform intensity with varying delta and h ####


if (exp.delta.h){
  # setting for simulation
  S = 1000
  n = 100
  r = 1
  delta = c(0.1, 0.05, 0.01)
  h = c(0.1, 0.05, 0.01, 0.005) 
  mean.exp.delta.h = sd.exp.delta.h = matrix(0, length(delta), length(h))
  
  no_cores = detectCores() - 2
  set.seed(1)
  for (i in 1:length(delta)) {
    for (j in 1:length(h)) {
      if (h[j] <= delta[i]){
        message(paste("S = ", S, " simulations with delta = ", delta[i], " and h = ", h[j], ".", sep = ""))
        # specification of non-uniform intensity
        L = augment.linnet(simplenet, delta[i], h[j])
        H = function(x, y, seg, tp) {sqrt(y)*exp(-x*y)*n/integral.linfun(linfun(function(x, y, seg, tp) {sqrt(y)*exp(-x*y)}, L))}
        varphi = linfun(H, L)
        
        cl = makeCluster(no_cores)
        clusterExport(cl, c("delta", "h", "r", "n", "L", "varphi", "i", "j",
                            "augment.linnet", "B.linnet", "gamma.hat.lpp",
                            "intensity.edge", "intensity.pspline.lpp",
                            "ise.intensity.splines", "ise.linfun", 
                            "IWLS", "K.linnet", "offset.bin", "simulation.study",
                            "squared.diff.linfun", "y.lpp"))
        
        clusterEvalQ(cl, library(spatstat)) 
        clusterEvalQ(cl, library(igraph))
        clusterEvalQ(cl, library(Matrix))   
        clusterEvalQ(cl, library(MASS))   
        clusterEvalQ(cl, library(VCA))   
        clusterEvalQ(cl, library(splines)) 
        
        ISE = parLapply(cl, 1:S, simulation.study, delta = delta[i], h = h[j], 
                        r = r, n = n, varphi = varphi, L = L, kernel = FALSE)
        
        stopCluster(cl) 
        
        mean.exp.delta.h[i, j] = mean(unlist(ISE)[which(names(unlist(ISE)) == "ISE.pspline")])
        sd.exp.delta.h[i, j] = sd(unlist(ISE)[which(names(unlist(ISE)) == "ISE.pspline")])
        save(mean.exp.delta.h, file = "RData/meanISEexpdeltah.RData")
        save(sd.exp.delta.h, file = "RData/sdISEexpdeltah.RData")
      }
    }
  }
}