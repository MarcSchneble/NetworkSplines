# adjust path
setwd("~/NetworkSplines")
rm(list = ls())

# load required packages
library(spatstat)
library(igraph)
library(Matrix)
library(splines) 
library(dplyr)
library(tidyr)
library(parallel)
library(ggplot2)
library(scales)


# load functions
source("Functions.R")

simulation.intensity.kernel <- TRUE
simulation.intensity.edges <- TRUE
simulation.internal <- FALSE
simulation.external <- FALSE

if (simulation.intensity.kernel | simulation.external){
  # augmented Chicago network
  delta <- 10
  h <- 2
  r <- 2
  L <- augment.linnet(as.linnet(chicago), delta, h, r)
  X <- as.ppp(chicago)
  L.lpp <- lpp(X, L)
  
  # intensity estimate with penalized spline method
  intens.pspline <- intensity.pspline.lpp(L.lpp)
  
  # intensity estimate with kernel based approach
  sigma <- bw.lppl(L.lpp, sigma = seq(50, 200, 0.1))
  intens.kernel <- density.lpp(L.lpp, sigma = as.numeric(sigma), dimyx = c(256, 256))
  
  # intensity with Voronoi estimate
  set.seed(1)
  f <- bw.voronoi(chicago)
  intens.voronoi <- densityVoronoi(chicago, f = as.numeric(f), nrep = 100, dimyx = c(256, 256))
  
  # plots
  pdf(file = "Plots/ChicagoNetwork.pdf", width = 10, height = 8)
  par(mar=c(0, 0, 2, 0), cex = 1.6)
  plot(L.lpp, main = "Chicago Crime Network", lwd = 3, leg.side = "right",
       cols = hue_pal()(7), chars = 0:6, size = 1)
  dev.off()
  
  max.intens <- max(intens.pspline$v, intens.kernel$v, na.rm = TRUE)
  
  pdf(file = "Plots/ChicagoIntensityPSpline.pdf", width = 10, height = 8)
  par(mar=c(0, 0, 3, 0.5), cex = 1.6)
  plot.linim(intens.pspline, main = "Penalized Spline Estimate" , zlim = c(0, max.intens))
  points(X)
  dev.off()
  
  pdf(file = "Plots/ChicagoIntensityKernel.pdf", width = 10, height = 8)
  par(mar=c(0, 0, 3, 0.5), cex = 1.6)
  plot.linim(intens.kernel, main = "Kernel Based Estimate", zlim = c(0, max.intens)) 
  points(X)
  dev.off()
  
  pdf(file = "Plots/ChicagoIntensityVoronoi.pdf", width = 10, height = 8)
  par(mar=c(0, 0, 3, 0.5), cex = 1.6)
  plot.linim(intens.voronoi, main = "Smoothed Voronoi Estimate") 
  points(X)
  dev.off()
  
  # intensity with marks
  
  fit.marks.pspline <- fit.lpp(L.lpp, lins = "marks")
  fit.marks.poisson <- lppm(L.lpp ~ marks, data = list(images = as.im(intens.kernel)))
  phi <- summary(fit.marks.poisson$fit$internal$glmfit)$dispersion
  
  df <- tibble(nr = c(seq(1, 11, 2), seq(2, 12, 2)),
               beta = c(tail(fit.marks.pspline$theta.hat, 6), fit.marks.poisson$fit$coef[2:7]),
               se = c(tail(fit.marks.pspline$se.hat, 6), sqrt(diag(vcov(fit.marks.poisson$fit$internal$glmfit)))[2:7]/sqrt(phi)),
               lower = beta - 1.96*se, upper = beta + 1.96*se,
               kind = factor(c(rep("pspline", 6), rep("poisson", 6)), levels = c("pspline", "poisson")),
               name = rep(c("burglary", "cartheft", "damage", "robbery", "theft", "trespass"), 2)) 
  
  g <- ggplot(df) + geom_point(aes(x = nr, y = beta, color = kind, shape = kind)) + 
    geom_errorbar(aes(x = nr, ymin = lower, ymax = upper, color = kind), width = 0.7, alpha = 0.7, show.legend = FALSE) + 
    theme_bw() + 
    scale_color_hue(name = "Model", labels = c("Penalized spline", "Poisson process")) +
    scale_shape_manual(name = "Model", labels = c("Penalized spline", "Poisson process"), values = c(16, 17)) +
    theme(legend.justification = c(0.01, 0.99), legend.position = c(0.01, 0.995), 
          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + 
    labs(x = "Kind of crime", y = "Effect size on the log-scale") + 
    scale_x_continuous(breaks = seq(1.5, 11.5, 2), labels = c("burglary", "cartheft", "damage", "robbery", "theft", "trespass")) + 
    scale_y_continuous(breaks = seq(-2.5, 1, 0.5)) 
  
  pdf(file = "Plots/ChicagoEffects.pdf", height = 4, width = 6)
  print(g)
  dev.off()
}



# simulation of network ----

if (simulation.intensity.kernel){
  
  R <- 100
  delta <- 10
  h <- 2
  r <- 2
  varphi <- as.linfun(intens.kernel)
  n <- c(50, 75, 100, 150, 200, 500, 1000)
  
  df.pspline <- tibble(n = as.factor(rep(n, each = R)), ISE = NA, kind = "pspline")
  df.kernel <- tibble(n = as.factor(rep(n, each = R)), ISE = NA, kind = "kernel")
  df.voronoi <- tibble(n = as.factor(rep(n, each = R)), ISE = NA, kind = "voronoi")
  
  for (i in 1:length(n)) {
    
    no_cores <- detectCores() - 2
    set.seed(1)
    
    
    cl <- makeCluster(no_cores)
    clusterExport(cl, ls())
    
    clusterEvalQ(cl, library(spatstat)) 
    clusterEvalQ(cl, library(igraph))
    clusterEvalQ(cl, library(Matrix))   
    clusterEvalQ(cl, library(MASS))   
    clusterEvalQ(cl, library(splines)) 
    clusterEvalQ(cl, library(dplyr))
    clusterEvalQ(cl, library(tidyr))
    
    message(n[i])
    ISE <- parLapply(cl, 1:R, simulation.chicago, delta = delta, h = h, r = r, n = n[i], varphi = varphi, kernel = TRUE, voronoi = TRUE)
    stopCluster(cl) 
    
    df.pspline$ISE[which(df.pspline$n == n[i])] <- unlist(ISE)[which(names(unlist(ISE)) == "ISE.pspline")]
    df.kernel$ISE[which(df.kernel$n == n[i])] <- unlist(ISE)[which(names(unlist(ISE)) == "ISE.kernel")] 
    df.voronoi$ISE[which(df.voronoi$n == n[i])] <- unlist(ISE)[which(names(unlist(ISE)) == "ISE.voronoi")]
  }
  
  df <- bind_rows(df.pspline, df.kernel, df.voronoi) %>% mutate(kind = factor(kind, levels = c("pspline", "kernel", "voronoi")))
  g <- ggplot(df, aes(x = n, y = ISE)) + 
    geom_boxplot(aes(color = kind)) + 
    theme_bw() + 
    theme(legend.justification = c(0.99, 0.99), legend.position = c(0.99, 0.995)) + 
    labs(color = "Estimate") + 
    scale_color_hue(labels = c("Penalized Spline", "Kernel Based", "Smoothed Voronoi"))
  
  pdf(file = "Plots/simulation_n.pdf", height = 6, width = 6)
  print(g)
  dev.off()
}

# simulation of network edges ----

if (simulation.intensity.edges){
  
  R <- 100
  delta <- 10
  h <- 2
  r <- 2
  H <- function(x, y, seg, tp){
    1*(seg %% 10 == 0)
  }
  varphi <- linfun(H, L)
  n <- c(50, 75, 100, 150, 200, 500, 1000)
  
  df.pspline <- tibble(n = as.factor(rep(n, each = R)), ISE = NA, kind = "pspline")
  df.kernel <- tibble(n = as.factor(rep(n, each = R)), ISE = NA, kind = "kernel")
  df.voronoi <- tibble(n = as.factor(rep(n, each = R)), ISE = NA, kind = "voronoi")
  
  for (i in 1:length(n)) {
    
    no_cores <- detectCores() - 2
    set.seed(1)
    
    
    cl <- makeCluster(no_cores)
    clusterExport(cl, ls())
    
    clusterEvalQ(cl, library(spatstat)) 
    clusterEvalQ(cl, library(igraph))
    clusterEvalQ(cl, library(Matrix))   
    clusterEvalQ(cl, library(MASS))   
    clusterEvalQ(cl, library(splines)) 
    clusterEvalQ(cl, library(dplyr))
    clusterEvalQ(cl, library(tidyr))
    
    message(n[i])
    ISE <- parLapply(cl, 1:R, simulation.chicago, delta = delta, h = h, r = r, n = n[i], varphi = varphi, kernel = TRUE, voronoi = TRUE)
    stopCluster(cl) 
    
    df.pspline$ISE[which(df.pspline$n == n[i])] <- unlist(ISE)[which(names(unlist(ISE)) == "ISE.pspline")]
    df.kernel$ISE[which(df.kernel$n == n[i])] <- unlist(ISE)[which(names(unlist(ISE)) == "ISE.kernel")] 
    df.voronoi$ISE[which(df.voronoi$n == n[i])] <- unlist(ISE)[which(names(unlist(ISE)) == "ISE.voronoi")]
  }
  
  df <- bind_rows(df.pspline, df.kernel, df.voronoi) %>% mutate(kind = factor(kind, levels = c("pspline", "kernel", "voronoi")))
  g <- ggplot(df, aes(x = n, y = ISE)) + 
    geom_boxplot(aes(color = kind)) + 
    theme_bw() + 
    theme(legend.justification = c(0.99, 0.99), legend.position = c(0.99, 0.995)) + 
    labs(color = "Estimate") + 
    scale_color_hue(labels = c("Penalized Spline", "Kernel Based", "Smoothed Voronoi"))
  
  pdf(file = "Plots/simulation_n.pdf", height = 6, width = 6)
  print(g)
  dev.off()
}


# simulation with network dependend covariates ----

if (simulation.internal){
  
  R <- 100
  delta <- 10
  h <- 1
  r <- 2
  L <- augment.linnet(as.linnet(chicago), delta, h, r)
  H <- function(x, y, seg, tp){
    exp(2*tp + x/1000)
  }
  varphi <- linfun(H, L)
  n <- c(50, 75, 100, 150, 200, 500, 1000)
  
  df.tp <- tibble(n = as.factor(rep(n, each = R)), beta = NA, kind = "tp")
  df.x <- tibble(n = as.factor(rep(n, each = R)), beta = NA, kind = "x")
  
  for (i in 1:length(n)) {
    
    no_cores <- min(detectCores() - 2, 10)
    set.seed(1)
    
    
    cl <- makeCluster(no_cores)
    clusterExport(cl, ls())
    
    clusterEvalQ(cl, library(spatstat)) 
    clusterEvalQ(cl, library(igraph))
    clusterEvalQ(cl, library(Matrix))   
    clusterEvalQ(cl, library(MASS))   
    clusterEvalQ(cl, library(splines)) 
    clusterEvalQ(cl, library(dplyr))
    clusterEvalQ(cl, library(tidyr))
    
    message(n[i])
    fit <- parLapply(cl, 1:R, simulation.chicago.covariates.internal, delta = delta, h = h, r = r, n = n[i], varphi = varphi)
    stopCluster(cl) 
    
    df.tp$beta[which(df.tp$n == n[i])] <- unlist(fit)[which(names(unlist(fit)) == "beta1")]
    df.x$beta[which(df.x$n == n[i])] <- unlist(fit)[which(names(unlist(fit)) == "beta2")] 
  }
  
  df <- bind_rows(df.tp, df.x) %>% mutate(kind = factor(kind, levels = c("tp", "x")))
  saveRDS(df, file = "df_simulation_internal.rds")
  g <- ggplot(df, aes(x = n, y = beta)) + 
    geom_boxplot(aes(color = kind)) + 
    theme_bw() + 
    theme(legend.justification = c(0.99, 0.99), legend.position = c(0.99, 0.995)) + 
    labs(color = "Model") + 
    scale_color_hue(labels = c("Relative position on the edge", "x-coordinate"))
  
  pdf(file = "Plots/simulation_internal.pdf", height = 6, width = 6)
  print(g)
  dev.off()
}






# simulation with external covariates ----

if (simulation.external){
  
  delta <- 10
  h <- 5
  r <- 2
  L <- augment.linnet(as.linnet(chicago), delta, h, r)
  beta <- c(2, 1, -1)
  varphi <- as.linfun(intens.pspline)
  
  a <- simulation.chicago.covariates.external(1, delta, h, r, beta, varphi)
  
}





