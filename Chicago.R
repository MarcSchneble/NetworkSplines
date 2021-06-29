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
library(gridExtra)
library(gtable)

# load functions
source("Functions.R")

simulation.intensity.kernel <- TRUE
simulation.intensity.edges <- FALSE
simulation.intensity.delta.h <- FALSE
simulation.internal <- FALSE
simulation.external <- FALSE

if (simulation.intensity.kernel | simulation.external | simulation.intensity.delta.h){
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
  sigma <- bw.lppl(L.lpp, distance = "path")
  intens.kernel <- density.lpp(L.lpp, sigma = as.numeric(sigma), dimyx = c(256, 256))
  
  # intensity estimate with two-dimensional kernel
  sigma <- bw.scott.iso(L.lpp)
  intens.kernel2d <- density.lpp(L.lpp, sigma = sigma, distance = "euclidean", dimyx = c(256, 256))
  
  # plots
  pdf(file = "Plots/ChicagoNetwork.pdf", width = 8, height = 5.5)
  par(mar=c(0, 0, 0, 0), cex = 1.6)
  plot(L.lpp, main = "", lwd = 3, leg.side = "right",
       cols = scales::hue_pal()(7), chars = 0:6, size = 1, leg.args = list(sep = 0.0001))
  dev.off()
  
  max.intens <- max(intens.pspline$v, intens.kernel$v, intens.kernel2d$v, na.rm = TRUE)
  
  pdf(file = "Plots/ChicagoIntensityPSpline.pdf", width = 8, height = 6.5)
  par(mar=c(0, 0, 0, 1), cex = 1.6)
  plot.linim(intens.pspline, main = "" , zlim = c(0, max.intens), ribsep = -0.05, 
             box = TRUE, ribwid = 0.08)
  points(X)
  dev.off()
  
  pdf(file = "Plots/ChicagoIntensityKernel.pdf", width = 8, height = 6.5)
  par(mar=c(0, 0, 0, 1), cex = 1.6)
  plot.linim(intens.kernel, main = "" , zlim = c(0, max.intens), ribsep = -0.05, 
             box = TRUE, ribwid = 0.08) 
  points(X)
  dev.off()
  
  pdf(file = "Plots/ChicagoIntensityKernel2d.pdf", width = 8, height = 6.5)
  par(mar=c(0, 0, 0, 1), cex = 1.6)
  plot.linim(intens.kernel2d, main = "" , zlim = c(0, max.intens), ribsep = -0.05, 
             box = TRUE, ribwid = 0.08) 
  points(X)
  dev.off()
  
  # intensity with marks
  
  fit.marks.pspline <- fit.lpp(L.lpp, lins = "marks")
  fit.marks.poisson <- lppm(L.lpp ~ marks, data = list(images = as.im(intens.kernel)))
  phi <- summary(fit.marks.poisson$fit$internal$glmfit)$dispersion
  
  df <- tibble(nr = c(seq(1, 11, 2), seq(2, 12, 2)),
               beta = c(tail(fit.marks.pspline$theta.hat, 6), fit.marks.poisson$fit$coef[2:7]),
               se = c(tail(fit.marks.pspline$effects$linear$se, 6), sqrt(diag(vcov(fit.marks.poisson$fit$internal$glmfit)))[2:7]/sqrt(phi)),
               lower = beta - 1.96*se, upper = beta + 1.96*se,
               kind = factor(c(rep("pspline", 6), rep("poisson", 6)), levels = c("pspline", "poisson"))) 
  
  g <- ggplot(df) + geom_point(aes(x = nr, y = beta, color = kind, shape = kind), size = 2.5) + 
    geom_errorbar(aes(x = nr, ymin = lower, ymax = upper, color = kind), width = 0.7, alpha = 0.7) + 
    theme_bw() + 
    scale_color_hue(labels = c("penalized spline based", "Poisson process")) +
    scale_shape(labels = c("penalized spline based", "Poisson process")) +
    theme(legend.justification = c(0.01, 0.99), legend.position = c(0.01, 0.995), 
          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + 
    labs(color = "Estimate", shape = "Estimate", 
         x = "Kind of crime", y = "Effect size on the log-scale\nand 95% confidence interval") + 
    scale_x_continuous(breaks = seq(1.5, 11.5, 2), labels = c("Burglary", "Cartheft", "Damage", "Robbery", "Theft", "Trespass")) + 
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
  varphi <- intens.kernel
  N <- c(100, 200, 500, 1000)
  
  #(test <- lapply(1:2, simulation.chicago.2, delta = delta, h = h, r = r, n = 200, varphi = varphi, kernel = TRUE, kernel2d = TRUE))
  
  intens <- vector("list", length(N))
  for (i in 1:length(N)) {
    
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
    
    message(N[i])
    intens[[i]] <- parLapply(cl, 1:R, simulation.chicago.2, delta = delta, h = h, r = r, n = N[i], varphi = varphi, kernel = TRUE, kernel2d = TRUE)
    stopCluster(cl) 
    
    saveRDS(intens, file = "simulation_intensity_n.rds")
  }
  #intens <- readRDS("simulation_intensity_n.rds")
  
  df.pspline <- tibble(r = rep(1:R, length(N)), n = as.factor(rep(N, each = R)), ISE = NA, kind = "pspline")
  df.kernel <- tibble(r = rep(1:R, length(N)), n = as.factor(rep(N, each = R)), ISE = NA, kind = "kernel")
  df.kernel2d <- tibble(r = rep(1:R, length(N)), n = as.factor(rep(N, each = R)), ISE = NA, kind = "kernel2d")
  
  mise <- tibble(n = rep(N, 3), kind = rep(c("pspline", "kernel", "kernel2d"), each = length(N)), mise = 0, mise.test = 0, var = 0, bias = 0)
  
  for (i in 1:length(N)) {
    varphi.n <- varphi/integral.linim(varphi)*N[i]
    A <- B <- C <- matrix(0, 256, 256)
    for (r in 1:R) {
      A <- A + intens[[i]][[r]]$intens.pspline$v
      B <- B + intens[[i]][[r]]$intens.kernel$v
      C <- C + intens[[i]][[r]]$intens.kernel2d$v
      df.pspline$ISE[which(df.pspline$r == r & df.pspline$n == N[i])] <- integral.linfun((intens[[i]][[r]]$intens.pspline - varphi.n)^2)/N[i]^2 
      df.kernel$ISE[which(df.kernel$r == r & df.kernel$n == N[i])] <- integral.linfun((intens[[i]][[r]]$intens.kernel - varphi.n)^2)/N[i]^2 
      df.kernel2d$ISE[which(df.kernel2d$r == r & df.kernel2d$n == N[i])] <- integral.linfun((intens[[i]][[r]]$intens.kernel2d - varphi.n)^2)/N[i]^2 
    }
    mise$mise[which(mise$n == N[i] & mise$kind == "pspline")] <- mean(df.pspline %>% filter(n == N[i], kind == "pspline") %>% pull(ISE))
    mise$mise[which(mise$n == N[i] & mise$kind == "kernel")] <- mean(df.kernel %>% filter(n == N[i], kind == "kernel") %>% pull(ISE))
    mise$mise[which(mise$n == N[i] & mise$kind == "kernel2d")] <- mean(df.kernel2d %>% filter(n == N[i], kind == "kernel2d") %>% pull(ISE))
    A <- A/R
    B <- B/R
    C <- C/R
    mean.intens.pspline <- mean.intens.kernel <- mean.intens.kernel2d <- varphi.n
    mean.intens.pspline$v <- A
    mean.intens.kernel$v <- B
    mean.intens.kernel2d$v <- C
    mise$bias[which(mise$n == N[i] & mise$kind == "pspline")] <- integral.linfun((mean.intens.pspline - varphi.n)^2)/N[i]^2
    mise$bias[which(mise$n == N[i] & mise$kind == "kernel")] <- integral.linfun((mean.intens.kernel - varphi.n)^2)/N[i]^2
    mise$bias[which(mise$n == N[i] & mise$kind == "kernel2d")] <- integral.linfun((mean.intens.kernel2d - varphi.n)^2)/N[i]^2
    
    for (r in 1:R) {
      mise$var[which(mise$n == N[i] & mise$kind == "pspline")] <- mise$var[which(mise$n == N[i] & mise$kind == "pspline")] + 
        1/R*integral.linfun((intens[[i]][[r]]$intens.pspline - mean.intens.pspline)^2)/N[i]^2
      mise$var[which(mise$n == N[i] & mise$kind == "kernel")] <- mise$var[which(mise$n == N[i] & mise$kind == "kernel")] + 
        1/R*integral.linfun((intens[[i]][[r]]$intens.kernel - mean.intens.kernel)^2)/N[i]^2
      mise$var[which(mise$n == N[i] & mise$kind == "kernel2d")] <- mise$var[which(mise$n == N[i] & mise$kind == "kernel2d")] + 
        1/R*integral.linfun((intens[[i]][[r]]$intens.kernel2d - mean.intens.kernel2d)^2)/N[i]^2
    }
  }
  mise$mise.test <- mise$var + mise$bias
  mise[, 3:6] <- mise[, 3:6]*10000
  
  df <- bind_rows(df.pspline, df.kernel, df.kernel2d) %>% 
    mutate(kind = factor(kind, levels = c("pspline", "kernel", "kernel2d")))
  
  g <- ggplot(df %>% mutate(ISE = 10000*ISE), aes(x = n, y = ISE)) + 
    geom_boxplot(aes(color = kind)) + 
    theme_bw() + 
    theme(legend.justification = c(0.01, 0.99), legend.position = c(0.01, 0.995)) + 
    labs(color = "Estimate") + 
    scale_y_continuous(limits = c(0, 0.5)) + 
    scale_color_hue(labels = c("penalized spline based estimate", "kernel estimate based on shortest path distance", "kernel estimate based on Euclidean distance"))
  
  pdf(file = "Plots/simulation_n.pdf", height = 6, width = 6)
  print(g)
  dev.off()
  
  df.mise <- tibble(n = factor(rep(mise$n, 2), labels = c("n = 100", "n = 200", "n = 500", "n = 1000")), kind = factor(rep(mise$kind, 2)),
                    error = c(mise$var, mise$bias),
                    kind2 = factor(c(rep("IVar", nrow(mise)), rep("ISBias", nrow(mise))), levels = c("IVar", "ISBias"))) %>% 
    mutate(kind = factor(kind, levels = c("pspline", "kernel", "kernel2d"))) 
  
  g.kernel <- ggplot(df.mise, aes(x = kind, y = error)) + 
    geom_col(aes(fill = kind2, linetype = kind), position = "stack", color = "black", size = 1, width = 0.8) + 
    scale_y_continuous(breaks = seq(0, 0.3, 0.05), limits = c(0, 0.3)) +
    theme_bw() + 
    facet_wrap( ~ n, nrow = 1) + 
    theme(panel.spacing = unit(0, "lines"), axis.text.x = element_blank(), 
          legend.position = "bottom", legend.direction = "vertical") + 
    guides(linetype = guide_legend(keyheight = 0.8, keywidth = 1.6, override.aes = list(fill = NA, col = "black")),
           fill = guide_legend(keyheight = 0.8, keywidth = 1.6)) +
    labs(y = "MISE", x = "", fill = "", linetype = "Method") + 
    scale_linetype(labels = c("penalized spline based estimate", "kernel estimate based on shortest path distance", "kernel estimate based on Euclidean distance"))
}

# simulation of network ----

if (simulation.intensity.edges){
  
  R <- 100
  delta <- 10
  h <- 2
  r <- 2
  L <- augment.linnet(as.linnet(chicago), delta, h, r)
  N <- c(100, 200, 500, 1000)
  
  H <- function(x, y, seg, tp){
    1*(seg %% 10 == 0)
  }
  varphi <- as.linim(linfun(H, L), dimyx = c(256, 256))
  varphi <- varphi/integral.linim(varphi)
  
  intens.edges <- vector("list", length(N))
  for (i in 1:length(N)) {
    
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
    
    message(N[i])
    intens.edges[[i]] <- parLapply(cl, 1:R, simulation.chicago.2, delta = delta, h = h, r = r, n = N[i], 
                                   varphi = varphi, kernel = TRUE, kernel2d = TRUE)
    stopCluster(cl) 
    
    saveRDS(intens.edges, file = "simulation_intensity_edges_n.rds")
  }
  #intens.edges <- readRDS("simulation_intensity_edges_n.rds")
  
  df.pspline <- tibble(r = rep(1:R, length(N)), n = as.factor(rep(N, each = R)), ISE = NA, kind = "pspline")
  df.kernel <- tibble(r = rep(1:R, length(N)), n = as.factor(rep(N, each = R)), ISE = NA, kind = "kernel")
  df.kernel2d <- tibble(r = rep(1:R, length(N)), n = as.factor(rep(N, each = R)), ISE = NA, kind = "kernel2d")
  
  mise <- tibble(n = rep(N, 3), kind = rep(c("pspline", "kernel", "kernel2d"), each = length(N)), mise = 0, mise.test = 0, var = 0, bias = 0)
  
  for (i in 1:length(N)) {
    varphi.n <- varphi/integral.linim(varphi)*N[i]
    A <- B <- C <- matrix(0, 256, 256)
    for (r in 1:R) {
      A <- A + intens.edges[[i]][[r]]$intens.pspline$v
      B <- B + intens.edges[[i]][[r]]$intens.kernel$v
      C <- C + intens.edges[[i]][[r]]$intens.kernel2d$v
      df.pspline$ISE[which(df.pspline$r == r & df.pspline$n == N[i])] <- integral.linfun((intens.edges[[i]][[r]]$intens.pspline - varphi.n)^2)/N[i]^2 
      df.kernel$ISE[which(df.kernel$r == r & df.kernel$n == N[i])] <- integral.linfun((intens.edges[[i]][[r]]$intens.kernel - varphi.n)^2)/N[i]^2 
      df.kernel2d$ISE[which(df.kernel2d$r == r & df.kernel2d$n == N[i])] <- integral.linfun((intens.edges[[i]][[r]]$intens.kernel2d - varphi.n)^2)/N[i]^2 
    }
    mise$mise[which(mise$n == N[i] & mise$kind == "pspline")] <- mean(df.pspline %>% filter(n == N[i], kind == "pspline") %>% pull(ISE))
    mise$mise[which(mise$n == N[i] & mise$kind == "kernel")] <- mean(df.kernel %>% filter(n == N[i], kind == "kernel") %>% pull(ISE))
    mise$mise[which(mise$n == N[i] & mise$kind == "kernel2d")] <- mean(df.kernel2d %>% filter(n == N[i], kind == "kernel2d") %>% pull(ISE))
    A <- A/R
    B <- B/R
    C <- C/R
    mean.intens.pspline <- mean.intens.kernel <- mean.intens.kernel2d <- varphi.n
    mean.intens.pspline$v <- A
    mean.intens.kernel$v <- B
    mean.intens.kernel2d$v <- C
    mise$bias[which(mise$n == N[i] & mise$kind == "pspline")] <- integral.linfun((mean.intens.pspline - varphi.n)^2)/N[i]^2
    mise$bias[which(mise$n == N[i] & mise$kind == "kernel")] <- integral.linfun((mean.intens.kernel - varphi.n)^2)/N[i]^2
    mise$bias[which(mise$n == N[i] & mise$kind == "kernel2d")] <- integral.linfun((mean.intens.kernel2d - varphi.n)^2)/N[i]^2
    
    for (r in 1:R) {
      mise$var[which(mise$n == N[i] & mise$kind == "pspline")] <- mise$var[which(mise$n == N[i] & mise$kind == "pspline")] + 
        1/R*integral.linfun((intens.edges[[i]][[r]]$intens.pspline - mean.intens.pspline)^2)/N[i]^2
      mise$var[which(mise$n == N[i] & mise$kind == "kernel")] <- mise$var[which(mise$n == N[i] & mise$kind == "kernel")] + 
        1/R*integral.linfun((intens.edges[[i]][[r]]$intens.kernel - mean.intens.kernel)^2)/N[i]^2
      mise$var[which(mise$n == N[i] & mise$kind == "kernel2d")] <- mise$var[which(mise$n == N[i] & mise$kind == "kernel2d")] + 
        1/R*integral.linfun((intens.edges[[i]][[r]]$intens.kernel2d - mean.intens.kernel2d)^2)/N[i]^2
    }
  }
  mise$mise.test <- mise$var + mise$bias
  mise[, 3:6] <- mise[, 3:6]*10000
  
  df <- bind_rows(df.pspline, df.kernel, df.kernel2d) %>% 
    mutate(kind = factor(kind, levels = c("pspline", "kernel", "kernel2d")))
  
  g <- ggplot(df %>% mutate(ISE = 10000*ISE), aes(x = n, y = ISE)) + 
    geom_boxplot(aes(color = kind)) + 
    theme_bw() + 
    theme(legend.justification = c(0.01, 0.99), legend.position = c(0.01, 0.995)) + 
    labs(color = "Estimate") + 
    scale_y_continuous(limits = c(0, 5)) + 
    scale_color_hue(labels = c("penalized spline based estimate", "kernel estimate based on shortest path distance", "kernel estimate based on Euclidean distance"))
  
  pdf(file = "Plots/simulation_n.pdf", height = 6, width = 6)
  print(g)
  dev.off()
  
  df.mise <- tibble(n = factor(rep(mise$n, 2), labels = c("n = 100", "n = 200", "n = 500", "n = 1000")), kind = factor(rep(mise$kind, 2)),
                    error = c(mise$var, mise$bias),
                    kind2 = factor(c(rep("IVar", nrow(mise)), rep("ISBias", nrow(mise))), levels = c("IVar", "ISBias"))) %>% 
    mutate(kind = factor(kind, levels = c("pspline", "kernel", "kernel2d"))) 
  
  g.edges <- ggplot(df.mise, aes(x = kind, y = error)) + 
    geom_col(aes(fill = kind2, linetype = kind), position = "stack", color = "black", size = 1, width = 0.8) + 
    scale_y_continuous(breaks = seq(0, 3, 0.5), limits = c(0, 3)) +
    theme_bw() + 
    facet_wrap( ~ n, nrow = 1) + 
    theme(panel.spacing = unit(0, "lines"), axis.text.x = element_blank(),
          legend.position = "bottom", legend.direction = "vertical") + 
    guides(linetype = guide_legend(keyheight = 0.8, keywidth = 1.6, override.aes = list(fill = NA, col = "black")),
           fill = guide_legend(keyheight = 0.8, keywidth = 1.6)) +
    labs(y = "MISE", x = "", fill = "", linetype = "Method") + 
    scale_linetype(labels = c("penalized spline based estimate", "kernel estimate based on shortest path distance", "kernel estimate based on Euclidean distance"))
  
  legend <- gtable_filter(ggplot_gtable(ggplot_build(g.kernel + theme(legend.position = "bottom"))), "guide-box") 
  
  pdf(file = "Plots/mise.pdf", height = 6, width = 12)
  grid.arrange(arrangeGrob(g.kernel + theme(legend.position="none"), g.edges + theme(legend.position="none"), nrow = 1), 
               legend, nrow = 2, heights = c(9, 2))
  dev.off() 
}


# simulation of network edges ----


if (simulation.intensity.delta.h) {
  R <- 100
  delta <- c(5, 10, 20, 50)
  r <- 2
  varphi <- intens.kernel
  n <- 200
  
  df <- tibble(delta = rep(delta, each = 3*R), h = delta/rep(rep(c(2, 5, 10), each = R), 3), ISE = NA)
  
  for (i in 1:length(delta)) {
    for (j in c(2, 5, 10)) {
      h <<- delta[i]/j
      
      no_cores <- min(detectCores() - 2, 20)
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
      
      message(c(delta[i], h))
      ISE <- parLapply(cl, 1:R, simulation.chicago, delta = delta[i], h = h, r = r, n = n, varphi = varphi)
      stopCluster(cl) 
      
      df$ISE[which(df$delta == delta[i] & df$h == h)] <- unlist(ISE)
    }
  }
  saveRDS(df, "df_simulation_delta_h.rds")
  df <- mutate(df, h = factor(h), delta = factor(delta))
  g <- list()
  d <- delta
  for (i in 1:length(delta)) {
    g[[i]] <- ggplot(df %>% mutate(ISE = ISE*10000) %>% filter(delta == d[i]), aes(x = h, y = ISE)) + 
      geom_boxplot(aes(color = h)) + 
      theme_bw() + 
      theme(legend.direction = "horizontal",
            legend.justification = c(0.01, 0.99), legend.position = c(0.01, 0.995)) +
      scale_y_continuous(limits = c(0, 0.5))
    labs(color = "Global bin width h")
    pdf(file = paste0("Plots/simulation_delta_", d[i], ".pdf"), height = 4, width = 4)
    print(g[[i]])
    dev.off()
  }
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
  n <- c(100, 200, 500, 1000)
  start <- Sys.time()
  test <- simulation.chicago.covariates.internal(1, delta, h, r, 200, varphi)
  print(Sys.time() - start)
  
  df.tp <- tibble(n = as.factor(rep(n, each = R)), beta = NA, kind = "tp")
  df.x <- tibble(n = as.factor(rep(n, each = R)), beta = NA, kind = "x")
  
  for (i in 1:length(n)) {
    
    no_cores <- min(detectCores() - 2, 20)
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
    print(fit)
    stopCluster(cl) 
    
    df.tp$beta[which(df.tp$n == n[i])] <- unlist(fit)[which(names(unlist(fit)) == "beta1")]
    df.x$beta[which(df.x$n == n[i])] <- unlist(fit)[which(names(unlist(fit)) == "beta2")] 
  }
  
  df <- bind_rows(df.tp, df.x) %>% mutate(kind = factor(kind, levels = c("tp", "x")))
  saveRDS(df, file = "df_simulation_internal.rds")
  df <- df[-which(df$beta < 0), ]
  g <- ggplot(df %>% filter(n %in% c(100, 200, 500, 1000)), aes(x = n, y = beta)) + 
    geom_boxplot(aes(color = kind)) + 
    theme_bw() + 
    theme(legend.justification = c(0.99, 0.99), legend.position = c(0.99, 0.995),
          legend.direction = "horizontal") + 
    labs(color = "Effect", y = expression(paste(hat(beta)))) + 
    scale_color_hue(labels = c("Relative position on the edge", "x-coordinate")) + 
    scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.5))
  
  pdf(file = "Plots/simulation_internal.pdf", height = 4, width = 6)
  print(g)
  dev.off()
}


# simulation with external covariates ----

if (simulation.external){
  
  R <- 100
  delta <- 10
  h <- 5
  r <- 2
  L <- augment.linnet(as.linnet(chicago), delta, h, r)
  beta <- list(
    c(2, 1, -1),
    c(1, -1, 1),
    c(1, 0, -1),
    c(4, 1, 1)
  )
  varphi <- as.linfun(intens.kernel)
  
  df.t <- tibble(i = rep(as.factor(1:4), each = R), beta = NA, kind = "t")
  df.type <- tibble(i = rep(as.factor(1:4), each = R), beta = NA, kind = "type")
  
  for (i in 1:length(beta)) {
    
    no_cores <- min(detectCores() - 2, 20)
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
    
    message(beta[[i]])
    fit <- parLapply(cl, 1:R, simulation.chicago.covariates.external, delta = delta, h = h, r = r, beta = beta[[i]], varphi = varphi)
    print(fit)
    stopCluster(cl) 
    
    df.t$beta[which(df.t$i == i)] <- unlist(fit)[which(names(unlist(fit)) == "beta1")]
    df.type$beta[which(df.type$i == i)] <- unlist(fit)[which(names(unlist(fit)) == "beta2")] 
  }
  
  df <- bind_rows(df.t, df.type) %>% mutate(kind = factor(kind, levels = c("t", "type")))
  saveRDS(df, file = "df_simulation_external.rds")
  g <- ggplot(df, aes(x = i, y = beta)) + 
    geom_boxplot(aes(color = kind)) + 
    theme_bw() + 
    theme(legend.direction = "horizontal") + 
    theme(legend.justification = c(0.99, 0.99), legend.position = c(0.99, 0.995)) + 
    labs(color = "Effect", y = expression(paste(hat(beta))),
         x = expression(beta)) + 
    scale_color_hue(labels = c("Time", "Type")) + 
    scale_x_discrete(labels = c(expression(beta^(1)), expression(beta^(2)), expression(beta^(3)), expression(beta^(4)))) + 
    scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 0.5))
  
  pdf(file = "Plots/simulation_external.pdf", height = 4, width = 6)
  print(g)
  dev.off()
}