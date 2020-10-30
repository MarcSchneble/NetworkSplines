reduce.linnet <- function(L){
  
  # degrees of vertices in the network
  V.deg <- degree(graph_from_adjacency_matrix(L$A, mode = "undirected"))
  
  # vertex indices with degree 2
  ind.deg.2 <- which(V.deg == 2)
  ind.v.null <- ind.deg.2
  
  # initialize
  ind.m.delete <- NULL
  new.d <- new.from <- new.to <- NULL
  i <- 0
  P <- list()
  
  # remove vertices with degree until none is left
  while(length(ind.deg.2) > 0){
    i <- i+1
    
    # find the two adjacent vertices to vertex with degree 2
    adj <- which(L$A[ind.deg.2[1], ] == 1)
    P[[i]] <- data.frame(from = c(adj[1], ind.deg.2[1]), 
                         to = c(ind.deg.2[1], adj[2]),
                         m = NA,
                         length = NA)

    # go into the direction of the first adjecent vertex and search for more
    # vertices with degree 2
    while(sum(L$A[adj[1], ]) == 2){
      adj.new <- which(L$A[adj[1], ] == 1)
      l <- adj.new[which(!is.element(adj.new, P[[i]]$from))]
      P[[i]] <- rbind(c(l, adj[1], NA, NA), P[[i]])
      adj[1] <- l
    }
    
    # go into the direction of the second adjecent vertex and search for more
    # vertices with degree 2
    while(sum(L$A[adj[2], ]) == 2){
      adj.new <- which(L$A[adj[2], ] == 1)
      r <- adj.new[which(!is.element(adj.new, P[[i]]$to))]
      P[[i]] <- rbind(P[[i]], c(adj[2], r, NA, NA))
      adj[2] <- r
    }
    
    # save the line indices and their corresponding lengths which are removed from the network
    for (k in 1:nrow(P[[i]])) {
      P[[i]]$m[k] <- which(L$from == P[[i]]$from[k] & L$to == P[[i]]$to[k] | L$from == P[[i]]$to[k] & L$to == P[[i]]$from[k])
      P[[i]]$length[k] <- L$d[P[[i]]$m[k]]
    }
    
    # add new connections in the adjacency matrix
    L$A[P[[i]]$from[1], P[[i]]$to[nrow(P[[i]])]] <- L$A[P[[i]]$to[nrow(P[[i]])], P[[i]]$from[1]] <- 1
    # add line segments to the delete vector
    ind.m.delete <- unique(c(ind.m.delete, P[[i]]$m))
    # compute new distances and the to/from vectors
    new.d <- c(new.d, sum(P[[i]]$length))
    new.from <- c(new.from, P[[i]]$from[1])
    new.to <- c(new.to, P[[i]]$to[nrow(P[[i]])])
    
    # remove vertices from the current vector of vertices with degree 2
    ind.deg.2 <- setdiff(ind.deg.2, P[[i]]$to[1:(nrow(P[[i]])-1)])
  }
  # update adjacency matrix
  L$A[ind.v.null, ] <- 0
  L$A[, ind.v.null] <- 0
  
  # delete line segments and add new line segments
  L$d <- c(L$d[- ind.m.delete], new.d)
  L$from <- c(L$from[- ind.m.delete], new.from)
  L$to <- c(L$to[- ind.m.delete], new.to)
  L$M <- length(L$d)
  
  # ensure that J.m > 0
  delta <- min(delta, min(L$d)/2)
  # line specific knot distances
  delta.m <- L$d/floor(L$d/delta)*(L$d/delta - floor(L$d/delta) < 0.5) + L$d/ceiling(L$d/delta)*(L$d/delta - floor(L$d/delta) >= 0.5)
  
  # ensure that h <= delta
  h <- min(h, delta)
  # line specific bin widths
  h.m <- L$d/floor(L$d/h)*(L$d/h - floor(L$d/h) < 0.5) + L$d/ceiling(L$d/h)*(L$d/h - floor(L$d/h) >= 0.5)
  
  # initializing...
  tau <- b <- z <- vector("list", L$M) 
  V.left <- V.right <- vector("list", L$W)
  N.m <- J.m <- rep(0, L$M)
  
  # do for every line segment
  for (m in 1:L$M) {
    
    # knot sequences tau
    tau[[m]] <- seq(0, L$d[m], delta.m[m])
    
    # bin boundaries b
    b[[m]] <- seq(0, L$d[m], h.m[m])
    
    # characterization of bins by midpoints z
    z[[m]] <- (b[[m]][1:(length(b[[m]])-1)] + b[[m]][2:length(b[[m]])])/2
    
    # total count of bins in the geometric network
    N.m[m] <- length(z[[m]])
    
    # count of linear B-splines on line segment
    J.m[m] <- length(tau[[m]]) - 2
  }
  
  # degrees of vertices in the network
  V.deg <- degree(graph_from_adjacency_matrix(L$A, mode = "undirected"))
  
  # incident lines to vertices
  # V.out and V.in instead of V.left and V.right
  for (v in 1:L$W) {
    vv <- which(L$A[v, ] == 1)
    if (length(vv) > 0){
      for (q in vv) {
        V.left[[v]] <- c(V.left[[v]], which(L$to == q & L$from == v))
        V.right[[v]] <- c(V.right[[v]], which(L$from == q & L$to == v))
      }
      V.left[[v]] <- sort(unique(V.left[[v]]))
      V.right[[v]] <- sort(unique(V.right[[v]]))
    } else {
      V.left[[v]] <- which(1 == 0)
      V.right[[v]] <- which(1 == 0)
    }
  }
  
  # replace in linnet object
  L$delta <- delta.m
  L$tau <- tau
  L$b <- b
  L$z <- z
  L$J.m <- J.m
  L$J <- sum(J.m) + L$W - length(ind.v.null)
  L$N.m <- N.m
  L$V.deg <- V.deg
  L$V.l <- V.left
  L$V.r <- V.right
  L$ind.v.null <- ind.v.null

  return(L)
}

augment.linnet = function(L, delta, h){
  # this function augments an object of class "linnet" by several attributes depending on delta and h
  
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

get.design <- function(L.lpp, r, smooths = NULL, lins = NULL, K.net = NULL, m = 10, l = 2){
  
  L <- as.linnet(L.lpp)
  K <- B <- ind.smooths <- vector("list", length(smooths) + 1)
  if (is.null(K.net)){
    K[[1]] <- K.linnet(L, r)
  } else {
    K[[1]] <- K.net
  }

  B[[1]] <- B.linnet(L)
  
  data <- bin.data(L.lpp, L, smooths = smooths, lins = lins)

  Z <- B[[1]][data$id, ]
  ind.smooths[[1]] <- 1:ncol(Z)
  names.Z <- paste0("L.", 1:ncol(Z))

  # add network covariates if needed
  if (is.element("dist2V", c(smooths, lins))){
    data$dist2V <- dist2V(L)[data$id]
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
      names.Z <- c(names.Z, paste0(smooths[a], ".", 1:ncol(B[[a + 1]])))
    }
  }

  # linear effects
  if (length(lins) > 0){
    for (a in 1:length(lins)) {
      Z <- cbind(Z, data %>% pull(as.symbol(lins[a])))
      names.Z <- c(names.Z, lins[a])
    }
  }
  
  return(list(data = data, Z = Z, K = K, ind.smooths = ind.smooths))
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
  covariates.comb <- covariates %>% distinct() %>% expand.grid() %>% as_tibble()
  N <- sum(L$N.m)*nrow(covariates.comb)
  
  # initializing
  data <- data.frame(matrix(NA, N, ncol(L.lpp$data) - 1))
  names(data) <- c("id", "y", "h", names(L.lpp$data)[-(1:4)])
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

dist2V <- function(L){
  d <- rep(NA, sum(L$N.m))
  ind <- 1
  for (m in 1:L$M) {
    d[ind:(ind + L$N.m[m] - 1)] <- pmin(L$z[[m]], L$d[m] - L$z[[m]])  
    ind <- ind + L$N.m[m]
  }
  return(d)
}

fit.lpp = function(L.lpp, r, smooths = NULL, lins = NULL, K.ginv.net = NULL,
                   rho = 10, rho.max = 1e5, eps.rho = 0.01, maxit.rho = 100){

  # get design of the model
  design <- get.design(L.lpp, r, smooths = smooths, lins = lins)
  
  # calculate generalized inverses
  K.ginv <- vector("list", length(design$K))
  if (is.null(K.ginv.net)) {
    K.ginv[[1]] <- ginv(as.matrix(Matrix(design$K[[1]])))
  } else {
    K.ginv[[1]] <- K.ginv.net
  }
  if (length(design$K) > 1) {
    for (a in 2:length(design$K)) {
      K.ginv[[a]] <- ginv(as.matrix(Matrix(design$K[[a]])))
    }
  }
  rho <- rep(rho, length(design$K))

  # determine optimal smoothing parameter rho with Fellner-Schall method
  Delta.rho <- Inf
  it.rho <- 0
  theta <- rep(0, ncol(design$Z))
  while(Delta.rho > eps.rho){
    print(rho)
    it.rho <- it.rho + 1
    fit <- optim(theta, fn = logL, gr = score, design = design, rho = rho, 
                       control = list(fnscale = -1, maxit = 1000),
                       method = "L-BFGS-B")
    theta <- fit$par
    V <- solve(fisher(theta, design, rho))
    
    # update rho 
    rho.new <- rep(NA, length(design$K))
    for (a in 1:length(design$K)) {
      rho.new[a] <- as.vector(rho[a]*(sum(diag(Matrix(K.ginv[[a]]/rho[a], sparse = TRUE)%*%design$K[[a]]))
                                - sum(diag(V[design$ind.smooths[[a]], design$ind.smooths[[a]]]%*%design$K[[a]])))/
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
    Delta.rho <- sqrt(sum((rho.new - rho)^2))/sqrt(sum((rho)^2))
    rho <- rho.new
  }
  return(list(theta.hat = theta, se.hat = sqrt(diag(V)), ind.smooths = design$ind.smooths))
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

intensity.pspline.lpp = function(L.lpp, r, smooths = NULL, lins = NULL, K.ginv.net = NULL){
  # returns an object of class "linim" according to the function density.lpp in the spatstat package
  
  # maximum likelihood estimate for gamma
  fit <- tryCatch(
    fit.lpp(L.lpp, r, smooths = smooths, lins = lins, K.ginv.net = K.ginv.net),
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
  L.im = pixellate(L)
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
  return(L.linim)
}

simulation.chicago <- function(s, delta, h, r, n, varphi, kernel = TRUE, K.ginv.net = NULL) {
  
  # simulate n points with intensity f
  L.lpp <- rlpp(n = n, f = varphi)
  
  # compute ISE for penalized spline based estimate
  # very rarely, singularities occur and cause an error
  ISE.pspline = tryCatch(sum(ise.intensity.splines(L.lpp, r, varphi, K.ginv.net = K.ginv.net))/n^2,
                         error = function(cond){
                           message(paste(cond), "")
                           return(0)
                         })
  if (kernel){
    ISE.kernel= sum(ise.linfun(as.linfun(density.lpp(L.lpp, sigma = as.numeric(bw.lppl(L.lpp)))), varphi, L))/n^2
    return(list(ISE.pspline = ISE.pspline, ISE.kernel = ISE.kernel))
  } else {
    return(list(ISE.pspline = ISE.pspline))
  }
}

intensity.edge = function(z, L, m, gamma.hat, n, varphi){
  # returns the squared difference of the estimated and the true intensity 
  
  B = matrix(0, length(z), L$J.m[m]+2)
  B[, 1] = (1-z/L$delta[m])*(1-z/L$delta[m] > 0)
  B[, 2:(ncol(B)-1)] = splineDesign(knots = L$tau[[m]], x = z, ord = 2, outer.ok = TRUE)
  B[, ncol(B)] = (1-(L$d[m]-z)/L$delta[m])*(1-(L$d[m]-z)/L$delta[m] > 0)
  
  gamma.hat.m = c(gamma.hat[sum(L$J.m)+L$from[m]],
                  gamma.hat[((cumsum(L$J.m)-L$J.m)[m]+1):cumsum(L$J.m)[m]],
                  gamma.hat[sum(L$J.m)+L$to[m]])
  if (min(z) >= 0 & max(z) <= L$d[m]) {
    return((as.vector(exp(B%*%gamma.hat.m)) - varphi(seg = m, tp = z/L$d[m]))^2)
  } else {
    stop("Wrong domain of z on this edge!")
  }
}

ise.intensity.splines = function(L.lpp, r, varphi, K.ginv.net = NULL){
  # returns edge-wise ISE for penalized spline based estimate
  
  n <- nrow(L.lpp$data)
  L <- as.linnet(L.lpp)
  fit <- fit.lpp(L.lpp, r, K.ginv.net = K.ginv.net)
  gamma.hat <- fit$theta.hat
  ISE <- rep(0, L$M)
  for (m in 1:L$M) {
    # very rarely, bad integral behavior causes an error
    ISE[m] <- tryCatch(integrate(f = intensity.edge, lower = 0, upper = L$d[m], L = L, m = m, 
                                gamma.hat = gamma.hat, n = n, varphi = varphi, subdivisions = 1000)$value,
                      error = function(cond){
                        message(paste(cond), "")
                        return(0)
                      })
  }
  return(ISE)
}

squared.diff.linfun = function(z, m, L, f1, f2){
  # returns squared difference of two functions f1 and f2 on L at z of segment m
  
  return((f1(tp = z/L$d[m], seg = m) - f2(tp = z/L$d[m], seg = m))^2)
}

ise.linfun = function(f1, f2, L){
  # returns edge-wise ISE of f1 and f2 on L
  
  ISE = rep(0, L$M)
  for (m in 1:L$M) {
    # very rarely, bad integral behavior causes an error
    ISE[m] = tryCatch(integrate(f = squared.diff.linfun, lower = 0, upper = L$d[m], m = m, 
                                L = L, f1 = f1, f2 = f2, subdivisions = 1000)$value,
                      error = function(cond){
                        message(paste(cond), "")
                        return(0)
                      })
  }
  return(ISE)
}

simulation.study = function(s, delta, h, r, n, varphi, L, kernel = TRUE){
  # function for parallel computing of penalized spline ISE and (if kernel = TRUE) kernel based ISE
  
  # simulate n points with intensity varphi
  L.lpp <- rlpp(n = n, varphi)
  
  # compute ISE for penalized spline based estimate
  # very rarely, singularities occur and cause an error
  ISE.pspline <- tryCatch(sum(ise.intensity.splines(L.lpp, r, varphi))/n^2,
                         error = function(cond){
                           message(paste(cond), "")
                           return(0)
                         })
  if (kernel){
    ISE.kernel <- sum(ise.linfun(as.linfun(density.lpp(L.lpp, sigma = as.numeric(bw.lppl(L.lpp)))), varphi, L))/n^2
    return(list(ISE.pspline = ISE.pspline, ISE.kernel = ISE.kernel))
  } else {
    return(list(ISE.pspline = ISE.pspline))
  }
}


logL <- function(theta, design, rho){
  mu <- exp(as.vector(design$Z%*%theta) + log(design$data$h))
  logL <- sum(design$data$y*log(mu) - mu)
  for (a in 1:length(design$K)) {
    logL <- logL - 0.5*rho[a]*
      as.vector(theta[design$ind.smooths[[a]]]%*%design$K[[a]]%*%theta[design$ind.smooths[[a]]])
  }
  return(logL)
}

score <- function(theta, design, rho){
  mu <- exp(as.vector(design$Z%*%theta) + log(design$data$h))
  score <- as.vector(colSums((design$data$y - mu)*design$Z))
  for (a in 1:length(design$K)) {
    score[design$ind.smooths[[a]]] <- 
      score[design$ind.smooths[[a]]] - rho[a]*design$K[[a]]%*%theta[design$ind.smooths[[a]]]  
  }
  
  return(score)
}

fisher <- function(theta, design, rho){
  mu <- exp(as.vector(design$Z%*%theta) + log(design$data$h))
  Mu <- Matrix(0, nrow(design$Z), nrow(design$Z))
  diag(Mu) <- mu
  fisher <- t(design$Z)%*%Mu%*%design$Z
  for (a in 1:length(design$K)) {
    fisher[design$ind.smooths[[a]], design$ind.smooths[[a]]] <- 
      fisher[design$ind.smooths[[a]], design$ind.smooths[[a]]] + rho[a]*design$K[[a]]
  }
  return(fisher)
}
