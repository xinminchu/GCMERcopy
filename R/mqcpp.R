##############################################################################
# Function for QCPP with ILP
# integer linear programming for Quasi-Clique Partition Problem
# with initialization method Multi-start greedy randomized heuristic (MSH)
# Reference:
# Melo et al. (2022) The minimum quasi-clique partitioning problem:
# Complexity, formulations, and a computational study
# Section 3 Integer programming formulations $3.1 Standard formulation
##############################################################################

mqcpp <- function(edges, gamma, n = NULL, ub = NULL)
{
  gurobi.flag <- require(gurobi)
  highs.flag <- require(highs)
  if (!(gurobi.flag | highs.flag))
  	stop("Please install the GUROBI commercial solver (gurobi.com)",
  		"or the 'highs' package from CRAN before using function 'mqcpp'.")
  		
  ## Preprocessing
  labels <- unique(c(edges))
  if (is.null(n)) 
  	n <- length(labels)
  n <- as.integer(n)
  if (!all(labels %in% (1:n)))
    edges <- matrix(match(edges, labels), ncol = 2)
  if (anyDuplicated(edges))
    edges <- edges[!duplicated(edges),]
  rmv <- which(edges[,1] == edges[,2])
  if (length(rmv) > 0) edges <- edges[-rmv,]
  flip <- which(edges[,1] > edges[,2])
  if (length(flip) > 0) edges[flip,] <- edges[flip,2:1]
  rm(rmv,flip)
  nedges <- nrow(edges)
  stopifnot(all(edges <= n))
  stopifnot(gamma > 0 && gamma <= 1)
  if (is.null(ub))
    ub <- floor(.5 + .5 * sqrt(1 + 8 * nedges / gamma))

  ## Define objective
  ny <- ub # number of quasi-clique (QC) indicator variables
  nx <- n * ub # number of QC membership indicator variables for vertices
  npairs <- choose(n,2)
  nw <- npairs * ub 
  # number of QC membership indicator variables for pairs of vertices
  nvars <- ny + nx + nw
  pairs <- t(combn(n, 2))
  # indexing: y (ub-vector), x (ub x n matrix), w (ub x npairs)
  objective <- numeric(nvars)
  objective[1:ub] <- 1

  ## Define constraints
  ncnstr <- c(n, nx, rep(nw, 3), ub, ub - 1)
  totcnstr <- sum(ncnstr) #csumcnstr[10]
  idxrow <- split(1:totcnstr, rep(1:7, ncnstr))
  A <- matrix(0, totcnstr, nvars) # each row = 1 constraint
  count <- 1
  
  ## Variable indices
  y <- 1:ny
  x <- matrix((ny+1):(ny+nx), ub, n)
  w <- matrix((ny+nx+1):nvars, ub, npairs)

  # Constraints (4) in paper sum_1:ub (x_iv) = 1
  Atmp <- matrix(0, n, nvars)
  for (v in 1:n)
    Atmp[v, x[,v]] <- 1
  A[idxrow[[1]],] <- Atmp

  # Constraints (5) in paper x_iv - y_i <= 0 ## choose '<' for all inequality
  Atmp <- array(0, c(ub, n, nvars))
  for (i in 1:ub) {
    for (v in 1:n) {
      Atmp[i,v, c(x[i,v], y[i])] <- c(1, -1)
    }}
  dim(Atmp) <- c(ub * n, nvars)
  A[idxrow[[2]],] <- Atmp

  # Constraints (6) x_iu + x_iv - w_iuv <= 1
  Atmp <- array(0, c(ub, npairs, nvars))
  for (i in 1:ub) {
    for (uv in 1:npairs) {
    	  u <- pairs[uv,1]
    	  v <- pairs[uv,2]
      Atmp[i, uv, c(x[i,u], x[i,v], w[i, uv])] <- c(1, 1, -1)
    }
  }
  dim(Atmp) <- c(ub * npairs, nvars)
  A[idxrow[[3]],] <- Atmp

  # Constraints (7) w_iuv <= x_iu
  # Constraints (8) w_iuv <= x_iv
  Atmp <- Btmp <- array(0, c(ub, npairs, nvars))
  for (i in 1:ub) {
    for (uv in 1:npairs) {
    	  u <- pairs[uv,1]
   	  v <- pairs[uv,2]
      Atmp[i, uv, c(w[i,uv], x[i,u])] <- c(1, -1)
      Btmp[i, uv, c(w[i,uv], x[i,v])] <- c(1, -1)
    }
  }
  dim(Atmp) <- dim(Btmp) <- c(ub * npairs, nvars)
  A[idxrow[[4]],] <- Atmp
  A[idxrow[[5]],] <- Btmp
  
  # Constraints (9) gamma * sum_{u<v} w_iuv - sum_{u<v, uv \in E} w_iuv <= 0
  Atmp <- matrix(0, ub, nvars)
  edgeschar <- paste(edges[,1], edges[,2], sep = ".")
  pairschar <- paste(pairs[,1], pairs[,2], sep = ".")
  idxedges <- match(edgeschar, pairschar)

  for (i in 1:ub) {
    Atmp[i,w[i,]] <- gamma
    Atmp[i,w[i,idxedges]] <- gamma - 1
  }
  A[idxrow[[6]],] <- Atmp
  
  # Constraints (10) y_i >= y_i+1 dual: y_i+1 - y_i <= 0
  Atmp <- matrix(0, ub-1,nvars)
  for (i in 1:(ub - 1)){
    Atmp[i, y[c(i,i+1)]] <- c(-1, 1)
  }
  A[idxrow[[7]],] <- Atmp
  
  rm(Atmp, edgeschar, pairschar)



  ## Run GUROBI
  if (gurobi.flag) {
	  model <- list(A = A, obj = objective,
	                modelsense = "min",
	                vtype = rep("B", nvars),
	                sense = rep(c("=", "<"), c(n, totcnstr-n)),
	                rhs = rep(c(1,0,1,0,0,0,0), ncnstr)
	                )
	  params <- list(OutputFlag = 0)
	  sol <- gurobi(model, params)
	  if (sol$status == "INFEASIBLE") {
	  	warning(sol$status)
	  	return(list(qc = NULL, nqc = NULL, dens = NULL))
	  }
  } else {
  ## Run 'highs'	
	  sol <- highs_solve(L = objective, A = A, 
		lower = rep(0, nvars), upper = rep(1, nvars), 
		types = rep("I", nvars),	
		lhs = c(rep(1,n), rep(-Inf,totcnstr-n)),
		rhs = rep(c(1,0,1,0,0,0,0), ncnstr))
	  if (sol$status_message != "Optimal") {
	  	warning(sol$status_message)
	  	return(list(qc = NULL, nqc = NULL, dens = NULL))	  	
	  }
  	  names_ <- names(sol)
  	  names_[names_ == "primal_solution"] <- "x"
  	  names(sol) <- names_
  }
  
  ## Format output
  y <- sol$x[y]
  nqc <- sum(y) # number of quasi-cliques
  x <- matrix(sol$x[x], ub, n)
  idx <- arrayInd(which(x == 1), dim(x))
  qc <- idx[,1] # integer vector indicating quasi-clique for each vertex
  w <- matrix(sol$x[w], dim(w))
  nedgeqc <- rowSums(w[1:nqc, idxedges, drop=FALSE])
  nvertqc <- rowSums(x[1:nqc,, drop=FALSE]) 
  dens <-  nedgeqc / choose(nvertqc,2) # quasi-clique density
  list(qc = qc, nqc = nqc, dens = dens) 
  	
}
