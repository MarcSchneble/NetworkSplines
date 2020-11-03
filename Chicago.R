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


# load functions
source("Functions.R")

# augmented Chicago network
delta <- 10
h <- 2
r <- 2
L <- augment.linnet(as.linnet(chicago), delta, h, r)
X <- as.ppp(unmark(chicago))
L.lpp <- lpp(X, L)

# intensity estimate with penalized spline method
intens.pspline <- intensity.pspline.lpp(L.lpp)

# intensity estimate with kernel based approach
sigma <- bw.lppl(L.lpp, sigma = seq(10, 500, 10))
intens.kernel <- density.lpp(L.lpp, sigma = as.numeric(sigma), dimyx = c(256, 256))

# intensity with Voronoi estimate
f <- as.numeric(bw.voronoi(chicago))
intens.voronoi <- densityVoronoi(chicago, f = f, nrep = 100, dimyx = c(256, 256))

# simulate
varphi <- as.linfun(intens.pspline)
a <- simulation.chicago(1, delta, h, r, 200, varphi, kernel = TRUE)

# simulation of network ----

R <- 50
delta <- 10
h <- 2
r <- 2
varphi <- as.linfun(intens.pspline)
n <- c(50, 100, 200, 500)

ISE.chicago <- data.frame(matrix(0, length(n), 5)) %>% as_tibble()
colnames(ISE.chicago) <- c("n", "Mean.PSpline", "SD.Spline", "Mean.Kernel", "SD.Kernel")
ISE.chicago$n <- n

for (i in 1:length(n)) {
  
  no_cores <- detectCores() - 2
  set.seed(1)
  
  
  cl <- makeCluster(no_cores)
  clusterExport(cl, ls())
  
  clusterEvalQ(cl, library(spatstat)) 
  clusterEvalQ(cl, library(igraph))
  clusterEvalQ(cl, library(Matrix))   
  clusterEvalQ(cl, library(MASS))   
  clusterEvalQ(cl, library(VCA))   
  clusterEvalQ(cl, library(splines)) 
  clusterEvalQ(cl, library(dplyr))
  clusterEvalQ(cl, library(tidyr))
  
  message(n[i])
  ISE <- parLapply(cl, 1:R, simulation.chicago, delta = delta, h = h, r = r, n = n[i], varphi = varphi, kernel = TRUE)
  
  stopCluster(cl) 
  
  ISE.chicago$Mean.PSpline[i] = mean(unlist(ISE)[which(names(unlist(ISE)) == "ISE.pspline")])
  ISE.chicago$SD.Spline[i] = sd(unlist(ISE)[which(names(unlist(ISE)) == "ISE.pspline")])
  ISE.chicago$Mean.Kernel[i] = mean(unlist(ISE)[which(names(unlist(ISE)) == "ISE.kernel")])
  ISE.chicago$SD.Kernel[i] = sd(unlist(ISE)[which(names(unlist(ISE)) == "ISE.kernel")])
}

# simulation with network dependend covariates ----

delta <- 10
h <- 1
r <- 2
L <- augment.linnet(as.linnet(chicago), delta, h, r)
H <- function(x, y, seg, tp){
  exp(2*tp + x/1000)
}
varphi <- linfun(H, L)
n <- 1000

a <- simulation.chicago.covariates.internal(1, delta, h, r, n, varphi)


# simulation with external covariates ----

delta <- 10
h <- 5
r <- 2
L <- augment.linnet(as.linnet(chicago), delta, h, r)
beta <- c(2, 1, -1)
varphi <- as.linfun(intens.pspline)

a <- simulation.chicago.covariates.external(1, delta, h, r, beta, varphi)



