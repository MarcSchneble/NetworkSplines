augment.linnet = function(L, delta, h, r, geninv = FALSE){
  # this function augments an object of class "linnet" by several attributes depending on delta, h and r
  
  # count of line segments
  M = L$lines$n
  
  # count of vertices
  W = L$vertices$n
  
  # adjacency matrix of vertices
  A = L$m*1
  
  # ends of line segments (x0 y0 x1 y1)
  ends = as.matrix(L$lines$ends)
  
  # length of line segments
  d = diag(L$dpath[L$from, L$to])
  # total length of the network
  d.L = sum(d)
  
  # ensure that J.m > 0
  #delta = min(delta, min(L$dpath[which(L$dpath > 0)])/2)
  # line specific knot distances
  delta <- d*(delta > d) + delta*(delta <= d)
  delta.m = pmin(d/floor(d/delta)*(d/delta - floor(d/delta) < 0.5) + d/ceiling(d/delta)*(d/delta - floor(d/delta) >= 0.5), d/2)
  
  # ensure that h <= delta
  h = min(h, min(delta.m))
  # line specific bin widths
  h.m = d/floor(d/h)*(d/h - floor(d/h) < 0.5) + d/ceiling(d/h)*(d/h - floor(d/h) >= 0.5)
  
  
  # initializing...
  tau = b = z = vector("list", M) 
  V.left = V.right = vector("list", W)
  N.m = J.m = rep(0, M)
  
  # do for every line segment
  for (m in 1:M) {
    
    # knot sequences tau
    tau[[m]] = seq(0, d[m], delta.m[m])
    
    # bin boundaries b
    b[[m]] = seq(0, d[m], h.m[m])
    
    # characterization of bins by midpoints z
    z[[m]] = (b[[m]][1:(length(b[[m]])-1)] + b[[m]][2:length(b[[m]])])/2
    
    # total count of bins in the geometric network
    N.m[m] = length(z[[m]])
    
    # count of linear B-splines on line segment
    J.m[m] = length(tau[[m]]) - 2
  }
  
  # degrees of vertices in the network
  V.deg = degree(graph_from_adjacency_matrix(A, mode = "undirected"))
  
  # incident lines to vertices
  i = 0
  for (v in 1:W) {
    i = i+1
    vv = which(A[v, ] == 1)
    q = 0
    for (j in vv) {
      q = q+1
      V.left[[i]] = c(V.left[[i]], which(L$to == vv[q] & L$from == v))
      V.right[[i]] = c(V.right[[i]], which(L$from == vv[q] & L$to == v))
    }
    V.left[[i]] = sort(unique(V.left[[i]]))
    V.right[[i]] = sort(unique(V.right[[i]]))
  }
  
  # add to linnet object
  L$M = M
  L$W = W
  L$d = d
  L$d.L = d.L
  L$delta = delta.m
  L$h = h.m
  L$tau = tau
  L$b = b
  L$z = z
  L$J.m = J.m
  L$J = sum(J.m)+W
  L$N.m = N.m
  L$A = A
  L$V.deg = V.deg
  L$V.l = V.left
  L$V.r = V.right
  L$ind.v.null <- NULL
  
  # add penalty matrix and generalized inverse (if required)
  K <- K.linnet(L, r)
  L$r <- r
  L$K <- K
  
  return(L)
}

B.linnet = function(L){
  # returns design matrix B with dimension N x J
  
  # design matrix for line segments
  B = Matrix(matrix(0, sum(L$N.m), L$J), sparse = TRUE)
  
  # line specific B-splines
  for (m in 1:L$M) {
    B[((cumsum(L$N.m)-L$N.m)[m]+1):cumsum(L$N.m)[m], ((cumsum(L$J.m)-L$J.m)[m]+1):cumsum(L$J.m)[m]] = 
      splineDesign(knots = L$tau[[m]], x = L$z[[m]], ord = 2, outer.ok = TRUE, sparse = TRUE)
  }
  
  vv <- 0
  # vertex specific B-splines
  for (v in setdiff(1:L$W, L$ind.v.null)) {
    vv <- vv + 1
    # left line ends
    for (m in L$V.l[[v]]) {
      B[((cumsum(L$N.m)-L$N.m)[m]+1):((cumsum(L$N.m)-L$N.m)[m] + length(which(1 - L$z[[m]]/L$delta[m] > 0))), sum(L$J.m)+vv] = 
        (1 - L$z[[m]]/L$delta[m])[which(1 - L$z[[m]]/L$delta[m] > 0)]
    }
    # right line ends
    for (m in L$V.r[[v]]) {
      B[(cumsum(L$N.m)[m]-length(which(1 - (L$d[m]-L$z[[m]])/L$delta[m] > 0))+1):cumsum(L$N.m)[m], sum(L$J.m)+vv] = 
        (1 - (L$d[m]-L$z[[m]])/L$delta[m])[which(1 - (L$d[m]-L$z[[m]])/L$delta[m] > 0)]
    }
  }
  # check if B is a valid design matrix
  if (sum(B) != sum(L$N.m)){
    stop("Error! Rowsums of B are not equal to one!")
  }
  return(B)
}

get.design <- function(L.lpp, smooths = NULL, lins = NULL, offset = NULL, m = 10, l = 2){
  
  L <- as.linnet(L.lpp)
  K <- B <- ind.smooths <- vector("list", length(smooths) + 1)
  K[[1]] <- L.lpp$domain$K
  B[[1]] <- B.linnet(L)
  
  data <- bin.data(L.lpp, L, smooths = smooths, lins = lins)
  
  # baseline intensity of the network
  Z <- B[[1]][data$id, ]
  ind.smooths[[1]] <- 1:ncol(Z)
  names.theta <- paste0("L.", 1:ncol(Z))
  
  # add network covariates if needed
  if (is.element("dist2V", c(smooths, lins))){
    data$dist2V <- dist2V(L)[data$id]
  }
  if (is.element("reldist2V", c(smooths, lins))){
    data$reldist2V <- dist2V(L, rel = TRUE)[data$id]
  }
  if (is.element("reldist2Vfrom", c(smooths, lins))){
    data$reldist2Vfrom <- tp(L)[data$id]
  }
  if (is.element("xcoord", c(smooths, lins))){
    data$xcoord <- get.x(L)[data$id]
  }
  if (is.element("x.km", c(smooths, lins))){
    data$x.km <- get.x.km(L)[data$id]
  }
  if (is.element("y.km", c(smooths, lins))){
    data$y.km <- get.y.km(L)[data$id]
  }
  if ("dist2Vdiscrete" %in% c(smooths, lins)){
    data$dist2Vdiscrete <- dist2Vdiscrete(L)[data$id]
  }
  if ("routetype" %in% lins){
    data$routetype <- L$routetype[data$id]
  } 
  if ("direction" %in% lins){
    data$direction <- L$direction[data$id]
  } 
  
  # smooth effects
  if (length(smooths) > 0){
    for (a in 1:length(smooths)) {
      data.a <- data %>% pull(as.symbol(smooths[a]))
      z <- sort(unique(data.a))
      delta <- (max(z) - min(z))/(m-1)
      k <- seq(min(z) - l*delta, max(z) + l*delta, delta)
      B.uncentered <- splineDesign(knots = k, x = z, ord = l+1, outer.ok = TRUE)
      # remove first column for identification
      B[[a + 1]] <- sweep(B.uncentered, 2, colMeans(B.uncentered))[, -1]
      # compute penalty
      D <- diff(diag(ncol(B[[a + 1]]) + 1), differences = 2)[, -1]
      K[[a + 1]] <- t(D)%*%D
      # add to design matrix
      Z <- cbind(Z, B[[a + 1]][match(data %>% pull(as.symbol(smooths[a])), z), ])
      ind.smooths[[a + 1]] <- (ncol(Z) - ncol(B[[a + 1]]) + 1):ncol(Z) 
      names.theta <- c(names.theta, paste0(smooths[a], ".", 1:ncol(B[[a + 1]])))
    }
  }
  
  # linear effects
  ind.lins <- NULL
  if (length(lins) > 0){
    ff <- as.formula(paste("~", paste(lins, collapse = " + ")))
    model.matrix.lin <- model.matrix(ff, data = data)
    Z <- cbind(Z, model.matrix.lin[, -1])
    ind.lins <- (ncol(Z) - ncol(model.matrix.lin) + 2):ncol(Z)
    names.theta <- c(names.theta, colnames(model.matrix.lin)[-1])
  }
  
  # offset
  if (!is.null(offset)){
    data$offset <- get.offset(offset, L)[data$id]
  } else {
    data$offset <- 1
  }
  
  return(list(data = data, Z = Z, B = B, K = K, ind.smooths = ind.smooths, ind.lins = ind.lins, names.theta = names.theta))
}

K.linnet = function(L, r){
  # returns the first or second penalty matrix K
  
  # adjacency matrix of knots
  A.tau = matrix(0, L$J, L$J)
  
  # on each line
  for (m in 1:L$M) {
    if (L$J.m[m] > 1){
      A.tau[((cumsum(L$J.m)-L$J.m)[m]+2):cumsum(L$J.m)[m], ((cumsum(L$J.m)-L$J.m)[m]+1):(cumsum(L$J.m)[m]-1)] = 
        diag(L$J.m[m]-1)
    }
  }
  
  # around each vertex
  for (v in 1:L$W) {
    # left line ends
    for (m in L$V.l[[v]]) {
      A.tau[sum(L$J.m)+v, (cumsum(L$J.m)-L$J.m)[m]+1] = 1
    }
    # right line ends
    for (m in L$V.r[[v]]) {
      A.tau[sum(L$J.m)+v, cumsum(L$J.m)[m]] = 1
    }
  }
  
  # filling the upper diagonal
  A.tau = A.tau + t(A.tau)
  
  if (r == 1){
    # find for every spline function the adjacent spline functions
    adj = which(A.tau*lower.tri(A.tau) == 1, arr.ind = T)
    
    # initalizing first order difference matrix
    D = Matrix(matrix(0, nrow(adj), L$J), sparse = TRUE)
    
    for(i in 1:nrow(adj)){
      D[i, adj[i, 1]] = 1
      D[i, adj[i, 2]] = -1
    }
    return(t(D)%*%D)
  }
  if (r == 2){
    # find for every spline function the adjacent spline functions
    adj1 = which(shortest.paths(graph_from_adjacency_matrix(A.tau)) == 1, arr.ind = T)
    adj2 = which(tril(shortest.paths(graph_from_adjacency_matrix(A.tau))) == 2, arr.ind = T)
    
    # initalizing second order difference matrix
    D = Matrix(matrix(0, nrow(adj2), L$J), sparse = TRUE)
    
    for (i in 1:nrow(adj2)) {
      D[i, adj2[i, 1]] = D[i, adj2[i, 2]] = 1
      D[i, intersect(adj1[which(adj1[, 1] == adj2[i, 1]), 2], adj1[which(adj1[, 1] == adj2[i, 2]), 2])] = -2
    }
    return(t(D)%*%D)
  }
  stop("Supply either r = 1 or r = 2")
}

bin.data <- function(L.lpp, L, smooths = NULL, lins = NULL){
  
  # get covariates from every point on the network (if applicable)
  name <- intersect(names(L.lpp$data), c(smooths, lins))
  covariates <- as.data.frame(L.lpp$data) %>% select(all_of(name))
  if (ncol(covariates) == 0){
    covariates <- covariates %>% mutate(q = 1)
  }
  
  # get all combinations of covariates and calculate the number of rows of the data matrix
  # get all combinations of covariates and calculate the number of rows of the data matrix
  covariates.comb <- covariates %>% distinct() %>% expand.grid() %>% distinct() %>% as_tibble()
  for (a in 1:length(covariates.comb)) {
    covariates.comb <- arrange(covariates.comb, !!sym(names(covariates.comb[a])))
  }
  
  N <- sum(L$N.m)*nrow(covariates.comb)
  
  # initializing
  data <- setNames(data.frame(matrix(nrow = N, ncol = length(name) + 3)), c("id", "y", "h", name)) %>% as_tibble() %>%
    mutate(id = as.integer(id), y = as.double(y), h = as.double(h))
  
  # set factor variable if applicable
  if (length(name) > 0){
    for (a in 1:length(name)) {
      if (is.factor(covariates.comb %>% pull(sym(name[a])))){
        data <- mutate(data, !!name[a] := factor(NA, levels = levels(covariates.comb %>% pull(sym(name[a])))))
      }
      if (is.double(covariates.comb %>% pull(sym(name[a])))){
        data <- mutate(data, !!name[a] := as.double(NA))
      } 
      if (is.integer(covariates.comb %>% pull(sym(name[a])))){
        data <- mutate(data, !!name[a] := as.integer(NA))
      } 
    }
  }

  ind <- 1
  for (j in 1:nrow(covariates.comb)) {
    ind.cov <- which(do.call(paste, covariates) == do.call(paste, covariates.comb[j, ]))
    data.sub <- L.lpp$data[ind.cov, ]
    for (m in 1:L$M) {
      # positions of data on line m
      ind.m <- which(data.sub$seg == m)
      y.m <- sort(as.numeric(data.sub[ind.m, ]$tp))*L$d[m]
      
      # bin data
      y.b <- rep(0, length(L$z[[m]]))
      for (k in 1:length(L$z[[m]])) {
        y.b[k] <- length(which(y.m < L$b[[m]][k+1] & y.m > L$b[[m]][k]))
      }
      
      # stack into one vector
      data$y[ind:(ind + length(y.b) - 1)] <- y.b
      data$h[ind:(ind + length(y.b) - 1)] <- L$h[m]
      ind <- ind + length(y.b)
    }
    # add bin id for every row
    data$id[((j-1)*sum(L$N.m) + 1):(j*sum(L$N.m))] <- 1:sum(L$N.m)
    
    # add covariates
    if (ncol(data) > 3){
      data[((j-1)*sum(L$N.m) + 1):(j*sum(L$N.m)), 4:ncol(data)] <- covariates.comb[j, ]
    }
  }
  return(data)
}

dist2V <- function(L, rel = FALSE){
  d <- rep(NA, sum(L$N.m))
  ind <- 1
  for (m in 1:L$M) {
    d[ind:(ind + L$N.m[m] - 1)] <- pmin(L$z[[m]], L$d[m] - L$z[[m]])/ifelse(rel, L$d[m], 1)
    ind <- ind + L$N.m[m]
  }
  return(d)
}

tp <- function(L){
  d <- rep(NA, sum(L$N.m))
  ind <- 1
  for (m in 1:L$M) {
    d[ind:(ind + L$N.m[m] - 1)] <- L$z[[m]]/L$d[m]
    ind <- ind + L$N.m[m]
  }
  return(d)
}

get.x <- function(L){
  x <- rep(NA, sum(L$N.m))
  ind <- 1
  for (m in 1:L$M) {
    x[ind:(ind + L$N.m[m] - 1)] <- as.lpp(seg = m, tp = L$z[[m]]/L$d[m], L = L)$data$x/1000
    ind <- ind + L$N.m[m]
  }
  return(x)
}

get.x.km <- function(L){
  x <- rep(NA, sum(L$N.m))
  ind <- 1
  for (m in 1:L$M) {
    x[ind:(ind + L$N.m[m] - 1)] <- as.lpp(seg = m, tp = L$z[[m]]/L$d[m], L = L)$data$x
    ind <- ind + L$N.m[m]
  }
  return(x)
}
get.y.km <- function(L){
  x <- rep(NA, sum(L$N.m))
  ind <- 1
  for (m in 1:L$M) {
    x[ind:(ind + L$N.m[m] - 1)] <- as.lpp(seg = m, tp = L$z[[m]]/L$d[m], L = L)$data$y
    ind <- ind + L$N.m[m]
  }
  return(x)
}

dist2Vdiscrete <- function(L, threshold = 20){
  d <- rep(NA, sum(L$N.m))
  ind <- 1
  for (m in 1:L$M) {
    d[ind:(ind + L$N.m[m] - 1)] <- ifelse(pmin(L$z[[m]], L$d[m] - L$z[[m]]) <= threshold, 1, 0)
    ind <- ind + L$N.m[m]
  }
  return(d)
}

get.routetype <- function(L, data){

}

get.offset <- function(offset, L){
  if ("linim" %in% class(offset)){
    L.linfun <- as.linfun(L.im)
    offset <- rep(NA, sum(L$N.m))
    ind <- 1
    for (m in 1:L$M) {
      offset[ind:(ind + length(L$z[[m]]) - 1)] <- L.linfun(seg = m, tp = L$z[[m]]/L$d[m])
      ind <- ind + length(L$z[[m]])
    }
    return(offset)
  }

  if ("lpp" %in% class(offset)){
    data <- bin.data(L.lpp, L)
    return(data$y + 1)
  }
  
  stop("Offset must be of class ``linim'' or ``lpp''.")

}

fit.lpp = function(L.lpp, smooths = NULL, lins = NULL, offset = NULL,
                   rho = 10, rho.max = 1e5, eps.rho = 0.01, maxit.rho = 100){
  
  # get design of the model
  design <- get.design(L.lpp, smooths = smooths, lins = lins, offset = offset)
  
  # determine optimal smoothing parameter rho with Fellner-Schall method
  rho <- rep(rho, length(design$K))
  Delta.rho <- Inf
  it.rho <- 0
  theta <- rep(0, ncol(design$Z))
  while(Delta.rho > eps.rho){
    it.rho <- it.rho + 1
    
    theta <- scoring(theta, design, rho)

    V <- solve(fisher(theta, design, rho))
    
    # update rho 
    rho.new <- rep(NA, length(design$K))
    for (a in 1:length(design$K)) {
      rho.new[a] <- as.vector(rho[a]*(rankMatrix(design$K[[a]], method = "qr.R")[1]/rho[a] - 
                                        sum(diag(V[design$ind.smooths[[a]], design$ind.smooths[[a]]]%*%design$K[[a]])))/
                                (t(theta[design$ind.smooths[[a]]])%*%design$K[[a]]%*%theta[design$ind.smooths[[a]]]))
    }
    
    if (any(rho.new > rho.max)) break
    if (any(rho.new < 0)){
      warning("rho = 0 has occurred")
    }
    if (it.rho > maxit.rho){
      warning("Stopped estimation of rho because maximum number of iterations has been reached!")
      break
    }
    print(rho.new)
    Delta.rho <- sqrt(sum((rho.new - rho)^2))/sqrt(sum((rho)^2))
    rho <- rho.new
  }
  
  # effects in one table
  effects <- list(linear = NULL, smooth = NULL)
  
  #linear effects
  if (length(lins) > 0) {
    effects$linear <- tibble(name = tail(design$names.theta, length(design$ind.lins)), estimate = NA, se = NA)
    for (i in 1:length(design$ind.lins)) {
      effects$linear$estimate[i] <- round(theta[design$ind.lins[i]], 3)
      effects$linear$se[i] <- round(sqrt(V[design$ind.lins[i], design$ind.lins[i]]), 3)
      effects$linear$rr[i] <- round(exp(effects$linear$estimate[i]), 2)
      effects$linear$rr.lower[i] <- round(exp(effects$linear$estimate[i] - 1.96*effects$linear$se[i]), 2)
      effects$linear$rr.upper[i] <- round(exp(effects$linear$estimate[i] + 1.96*effects$linear$se[i]), 2)
    }
  } 

  # smooth effects
  effects$smooth <- vector("list", length(smooths))
  names(effects$smooth) <- smooths
  if (length(smooths) > 0){
    for (i in 1:length(smooths)) {
      ind <- match(unique(design$data[[smooths[i]]]), design$data[[smooths[i]]])
      confidence.band <- get.confidence.band(theta, V, design, i, ind, smooths)
      effects$smooth[[i]] <- tibble(x = design$data[[smooths[i]]][ind],
                                    y = as.vector(design$Z[, design$ind.smooths[[i + 1]]]%*%theta[design$ind.smooths[[i + 1]]])[ind],
                                    lwr = confidence.band$lower,
                                    upr = confidence.band$upper)
    }
  }
  
  return(list(theta.hat = theta, V = V, ind.smooths = design$ind.smooths, names.theta = design$names.theta, effects = effects))
}

B.new.linnet = function(L, z.new, N.m.new){
  # returns design matrix for new data
  
  # design matrix for line segments
  B = matrix(0, sum(N.m.new), sum(L$J.m)+L$W)
  
  # line specific B-splines
  for (m in 1:L$M) {
    if (nrow(z.new[[m]]) > 0) {
      B[((cumsum(N.m.new)-N.m.new)[m]+1):cumsum(N.m.new)[m], ((cumsum(L$J.m)-L$J.m)[m]+1):cumsum(L$J.m)[m]] = 
        splineDesign(knots = L$tau[[m]], x = z.new[[m]]$pos, ord = 2, outer.ok = TRUE)
    }
  }
  
  # vertex specific B-splines
  for (v in 1:L$W) {
    # left line ends
    for (m in L$V.l[[v]]) {
      if (length(which(1 - z.new[[m]]$pos/L$delta[m] > 0)) > 0){
        B[((cumsum(N.m.new)-N.m.new)[m]+1):((cumsum(N.m.new)-N.m.new)[m] + length(which(1 - z.new[[m]]$pos/L$delta[m] > 0))), sum(L$J.m)+v] = 
          (1 - z.new[[m]]$pos/L$delta[m])[which(1 - z.new[[m]]$pos/L$delta[m] > 0)]
      }
      
    }
    # right line ends
    for (m in L$V.r[[v]]) {
      if (length(which(1 - (L$d[m]-z.new[[m]]$pos)/L$delta[m] > 0)) > 0){
        B[(cumsum(N.m.new)[m]-length(which(1 - (L$d[m]-z.new[[m]]$pos)/L$delta[m] > 0))+1):cumsum(N.m.new)[m], sum(L$J.m)+v] = 
          (1 - (L$d[m]-z.new[[m]]$pos)/L$delta[m])[which(1 - (L$d[m]-z.new[[m]]$pos)/L$delta[m] > 0)]
      }
    }
  }
  # check if B is a valid design matrix
  if (sum(B) != sum(N.m.new)){
    stop("Error! Rowsums of B are not equal to one!")
  }
  return(B)
}

intensity.pspline.lpp = function(L.lpp, smooths = NULL, lins = NULL, offset = NULL, dimyx = c(256, 256),
                                 eps.rho = 1e-3){
  # returns an object of class "linim" according to the function density.lpp in the spatstat package
  
  # maximum likelihood estimate for gamma
  fit <- tryCatch(
    fit.lpp(L.lpp, smooths = smooths, lins = lins, offset = offset, eps.rho = eps.rho),
    error = function(cond) {
      message("An error has occurred. Here's the original error message:")
      message(cond)
      message("\nReduce the value of delta and try again!")
      return(NULL)
    }
  )
  if (is.null(fit)){
    return(NULL)
  }
  
  gamma.hat.net <- fit$theta.hat[fit$ind.smooths[[1]]]
  
  # network represented by different classes
  L = as.linnet(L)
  L.im = pixellate(L, dimyx = dimyx)
  L.psp = as.psp(L.lpp)
  
  # pixels for plotting
  pixels = which(L.im$v > 0, arr.ind = T)
  
  # x and y coordinates of the pixels
  x.coord = y.coord = rep(0, nrow(pixels))
  for (i in 1:nrow(pixels)) {
    x.coord[i] = L.im$xcol[pixels[i, 2]]
    y.coord[i] = L.im$yrow[pixels[i, 1]]
  }
  # points which need to be colored by intensity
  L.ppp = ppp(x.coord, y.coord, window = L$window)
  
  # project points on network
  projection = project2segment(L.ppp, L.psp)
  
  # findest nearest line segment for each image pixel and the location of the point on the line segment
  line.nr = projection$mapXY
  line.pos = projection$tp*L$d[line.nr]
  
  # create data frame
  pixels = data.frame(matrix(c(x.coord, y.coord, line.nr, line.pos), length(line.nr), 4))
  colnames(pixels) = c("x", "y", "m", "pos")
  
  # pixels by line segment from left end to right end
  z.new = vector("list", L$M)
  N.m.new = rep(0, L$M)
  for (m in 1:L$M) {
    z.new[[m]] = pixels[which(pixels$m == m), ]
    z.new[[m]] = z.new[[m]][order(z.new[[m]]$pos), ]
    N.m.new[m] = nrow(z.new[[m]])
  }
  
  # design matrix for pixel locations
  B.new = B.new.linnet(L, z.new, N.m.new)
  
  # matrix with pixel intensities
  X = matrix(NA, L.im$dim[1], L.im$dim[2])
  for (m in 1:L$M) {
    intens.m = exp(B.new[((cumsum(N.m.new)-N.m.new)[m]+1):cumsum(N.m.new)[m], ]%*%gamma.hat.net)
    for (i in 1:nrow(z.new[[m]])) {
      X[which(L.im$yrow == z.new[[m]]$y[i]), which(L.im$xcol == z.new[[m]]$x[i])] = intens.m[i]
    }
  }
  
  # create "linim" object 
  L.linim = as.linim(X, L)
  L.linim$effects <- fit$effects
  
  return(L.linim)
}


simulation.chicago <- function(s, delta, h, r, n, varphi, kernel = FALSE, kernel2d = FALSE) {
  
  # simulate n points with intensity f
  L.lpp <- rlpp(n = n, f = varphi)
  
  # compute estimates and ISE
  intens.pspline <- intensity.pspline.lpp(L.lpp)
  if (kernel) {
    sigma <- as.numeric(bw.lppl(L.lpp, distance = "path"))
    intens.kernel <- density.lpp(L.lpp, sigma = sigma, dimyx = c(256, 256))
  } else {
    intens.kernel <- NULL
  }
  if (kernel2d){
    sigma <- bw.scott.iso(L.lpp)
    intens.kernel2d <- density.lpp(L.lpp, sigma = sigma, distance = "euclidean", dimyx = c(256, 256))
  } else {
    intens.kernel2d <- NULL
  }
  return(list(intens.pspline = intens.pspline, intens.kernel = intens.kernel, intens.kernel2d = intens.kernel2d))
}

simulation.chicago.covariates.internal <- function(s, delta, h, r, n, varphi){
  L.lpp <- rlpp(n, varphi)
  lins <- c("reldist2Vfrom", "xcoord")
  fit <- fit.lpp(L.lpp, lins = lins, rho = 100, rho.max = 1e4)
  return(list(beta = tail(fit$theta.hat, 2), se = fit$effects$linear$se))
}

simulation.chicago.covariates.external <- function(s, delta, h, r, beta, varphi){
  
  cov.t <- cov.type <- NULL
  L.lpp <- rlpp(n = 0, f = varphi)
  for (i in 1:10) {
    t <- round(rnorm(1), 2)
    for (j in c("A", "B")) {
      mu <- exp(beta[1] + t*beta[2] + beta[3]*(j == "B"))
      n <- rpois(1, mu)
      cov.t <- c(cov.t, rep(t, n))
      cov.type <- c(cov.type, rep(j, n))
      L.lpp <- superimpose.lpp(L.lpp, rlpp(n = n, f = varphi))
    }
  }
  L.lpp$data$t <- cov.t
  L.lpp$data$type <- as.factor(cov.type)
  
  lins <- c("t", "type")
  fit <- tryCatch(
    fit.lpp(L.lpp, lins = lins, rho = 10, eps.rho = 1e-2),
    error = function(cond){
      return(NULL)
    }
  ) 
  return(list(beta = tail(fit$theta.hat, 2), se = fit$effects$linear$se))
}

logL <- function(theta, design, rho){
  mu <- exp(as.vector(design$Z%*%theta) + log(design$data$h) + log(design$data$offset))
  logL <- sum(design$data$y*log(mu) - mu)
  for (a in 1:length(design$K)) {
    logL <- logL - 0.5*rho[a]*
      as.vector(theta[design$ind.smooths[[a]]]%*%design$K[[a]]%*%theta[design$ind.smooths[[a]]])
  }
  return(logL)
}

score <- function(theta, design, rho){
  mu <- exp(as.vector(design$Z%*%theta) + log(design$data$h) + log(design$data$offset))
  score <- as.vector(colSums((design$data$y - mu)*design$Z))
  for (a in 1:length(design$K)) {
    score[design$ind.smooths[[a]]] <- 
      score[design$ind.smooths[[a]]] - rho[a]*design$K[[a]]%*%theta[design$ind.smooths[[a]]]  
  }
  
  return(score)
}

fisher <- function(theta, design, rho){
  mu <- exp(as.vector(design$Z%*%theta) + log(design$data$h) + log(design$data$offset))
  Mu <- Matrix(0, nrow(design$Z), nrow(design$Z))
  diag(Mu) <- mu
  fisher <- t(design$Z)%*%Mu%*%design$Z
  for (a in 1:length(design$K)) {
    fisher[design$ind.smooths[[a]], design$ind.smooths[[a]]] <- 
      fisher[design$ind.smooths[[a]], design$ind.smooths[[a]]] + rho[a]*design$K[[a]]
  }
  return(fisher)
}

get.confidence.band <- function(theta, V, design, i, ind, smooths, q = 0.05, R = 10000){
  gamma <- theta[design$ind.smooths[[i + 1]]]
  cov <- V[design$ind.smooths[[i + 1]], design$ind.smooths[[i + 1]]]
  x <- unique(design$data[[smooths[i]]])
  cov <- as.matrix(Matrix(cov))
  
  set.seed(1)
  mu.sim <- matrix(0, R, length(x))
  for (j in 1:R) {
    gamma.sim <- rmvn(1, gamma, cov)
    mu.sim[j, ] <- design$B[[i + 1]]%*%gamma.sim
  }
  lower <- upper <- rep(0, ncol(mu.sim))
  for (j in 1:ncol(mu.sim)) {
    lower[j] <- quantile(mu.sim[, j], probs = q/2)
    upper[j] <- quantile(mu.sim[, j], probs = 1-q/2)
  }
  return(list(lower = lower, upper = upper))
}

adjust.for.vertex.distance <- function(L.linim, L, dist = 20){

  # network represented by different classes
  L.psp <- as.psp(L)
  
  # pixels for plotting
  pixels <- which(L.linim$v > 0, arr.ind = T)
  
  # x and y coordinates of the pixels
  x.coord <- y.coord <- rep(0, nrow(pixels))
  for (i in 1:nrow(pixels)) {
    x.coord[i] <- L.linim$xcol[pixels[i, 2]]
    y.coord[i] <- L.linim$yrow[pixels[i, 1]]
  }
  # points which need to be colored by intensity
  L.ppp <- ppp(x.coord, y.coord, window = L$window)
  
  # project points on network
  projection <- project2segment(L.ppp, L.psp)
  
  # findest nearest line segment for each image pixel and the location of the point on the line segment
  dist.vertex <- pmin(projection$tp*L$d[projection$mapXY], L$d[projection$mapXY] - projection$tp*L$d[projection$mapXY])
  
  ind <- which(dist.vertex < dist)
  
  for (i in ind) {
    L.linim$v[pixels[i, 1], pixels[i, 2]] <- L.linim$v[pixels[i, 1], pixels[i, 2]]*exp(L.linim$effects$linear$estimate[1])
  }
  
  return(L.linim)
}

scoring <- function(theta, design, rho, eps_theta = 1e-5){
  # perform iterative least squares estimation for a Poisson model with offset
  
  Delta_theta <- Inf
  it <- 0
  while (Delta_theta > eps_theta) {
    it <- it + 1
    theta_new <- as.vector(theta + Matrix::solve(fisher(theta, design, rho))%*%score(theta, design, rho))
    Delta_theta <- as.numeric(sqrt(t(theta - theta_new)%*%(theta - theta_new)))/
      as.numeric(sqrt(t(theta)%*%theta))
    theta <- as.vector(theta_new)
  }
  theta
}
