#######################################
# Function to find connected components
# from a graph represented by
# symmetric adjacency matrix
#######################################

find_connected_components <- function(x)
{
  stopifnot(is.matrix(x))
  stopifnot(nrow(x) == ncol(x))
  stopifnot(all(x == t(x)))
  n <- ncol(x)
  if (!is.logical(x)) x <- matrix(as.logical(x), n, n)
  if (any(diag(x))) diag(x) <- FALSE
  comps <- integer(n)
  unmarked <- rep(TRUE, n)
  ncomps <- 0
  while (any(comps == 0)) {
    ncomps <- ncomps + 1
    i <- which(unmarked)[1]
    comps[i] <- ncomps
    unmarked[i] <- FALSE	
    pool <- x[,i] & unmarked
    while (any(pool)) {	
      comps[pool] <- ncomps		
      i <- which(pool)[1]
      unmarked[i] <- FALSE
      pool <- (pool | x[,i]) & unmarked
    }
  }
  comps <- match(comps, unique(comps))
  comps
}


#######################################
# Function to create confusion matrix from two 
# integer vectors representing clusterings or 
# from their contingency table
#######################################

get_confusion <- function(x, y = NULL) {
	contingency <- if (is.null(y)) x else table(x,y)
	n <- sum(contingency)
	n11 <- sum(choose(contingency,2))
	n01 <- sum(choose(colSums(contingency),2)) - n11
	n10 <- sum(choose(rowSums(contingency),2)) - n11
	n00 <- choose(n,2) - n01 - n10 - n11
	confusion <- matrix(c(n00,n10,n01,n11), 2, 2)
	confusion
}



#######################################
# Function to retrieve entity/clusters
# from assignment vector
#######################################


get_entity <- function(entvec)
{
  split(seq_along(entvec), entvec)
}


#######################################
# Function to get independent set
# from adjacency matrix of a graph
#######################################

get_independent_set <- function(adjmat)
{
  n <- NCOL(adjmat) 
  if (n == 1) return(1L)
  
  ## Initialize independent set 
  ## with vertex of highest degree 
  degree <- colSums(adjmat == 1)
  colored <- which.max(degree)
  uncolored <- (1:n)[-colored]
  
  ## Partition uncolored vertices into those that have 
  ## no colored neighbor and those that do	
  has_colored_neighbor <- 
    (colSums(adjmat[colored, uncolored, drop = FALSE]) > 0) 
  U1 <- uncolored[!has_colored_neighbor]
  U2 <- uncolored[has_colored_neighbor]
  
  ## Add subsequent vertices 
  while (length(U1) > 0) {
    ## Numbers of neighbors in U2 for vertices of U1
    dU2 <- colSums(adjmat[U2, U1, drop = FALSE] == 1) 
    ## Candidates for coloring
    idx <- which(dU2 == max(dU2))
    if (length(idx) > 1) { # break ties 
      dU1 <- colSums(adjmat[U1, U1[idx]])
      idx <- idx[which.min(dU1)]
    }
    i <- U1[idx]
    uncolored <- setdiff(uncolored, i)
    U1 <- U1[-idx]
    idx <- which(adjmat[i,U1] == 1)
    if (length(idx) > 0) {
      U2 <- c(U2, U1[idx])			
      U1 <- U1[-idx]
    }
  }
  colored <- if (length(uncolored) > 0) {
    (1:n)[-uncolored] } else (1:n)
  colored	
}

##########################################
# Function to generate adjacency matrix
# for weighted or unweighted random graph
##########################################


make_adjmat <- function(n, weight = c("none", "random")){
  adjmat <- matrix(0, n, n) # create an adjacency matrix
  idx <- sample(1:(n^2), 2 * n) # randomly create edges
  if(weight == "random"){
    weight <- runif(2*n)
    adjmat[idx] <- weight
  }else{
    adjmat[idx] <- 1
  }
  diag(adjmat) <- 0
  adjmat <- adjmat + t(adjmat) # make sure it's symmetric
  adjmat[adjmat > 1] <- 1
  return(adjmat)
}



########################################
# Function to generate adjacency matrix
# of random Erdos-Renyi graph
########################################


make_er_graph <- function(n, p = 0.5) {
if (n == 1)
	return(matrix(0, 1, 1))
adjmat <- matrix(0L, n, n)
lower <- which(lower.tri(adjmat))
edge <- as.logical(rbinom(choose(n,2), 1, p))
adjmat[lower[edge]] <- 1L
adjmat + t(adjmat)
}


####################################
# Function to generate random graph
# with given chromatic number
####################################

make_graph_chrom <- function(n, k, p = 0.5)
{
color <- sample.int(k, n, TRUE)
sets <- split(1:n, color) # independent sets
adjmat <- make_er_graph(n, p)
for (i in 1:k)
	adjmat[sets[[i]], sets[[i]]] <- 0L
idx <- t(combn(sapply(sets, "[", 1), 2))
adjmat[idx] <- 1L
adjmat[idx[,2:1]] <- 1L
list(graph = adjmat, k = k, sets = sets)
}



###################################
# Function to get number of colors
# in graph coloring
###################################

get_chromatic <- function(color) length(unique(color))


###################################
# Function to check validity of a
# graph coloring
###################################


is_valid_coloring <- function(adjmat, color)
{
stopifnot(is.list(adjmat) || is.matrix(adjmat))
n <- length(color)
if (is.list(adjmat)) {
	stopifnot(length(adjmat) == length(color))
	for (i in 1:n) {
		if (color[i] %in% color[adjmat[[i]]]) return(FALSE)
	}
	return(TRUE)
}
stopifnot(ncol(adjmat) == length(color))
if (!is.logical(adjmat)) adjmat <- (adjmat == 1)
for (i in 1:n) {
	neighbors <- which(adjmat[,i])
	if (color[i] %in% color[neighbors]) return(FALSE)
}
return(TRUE)
}


##############################################
# Metric learning: function to return weights
# such that linear combination of features
# has minimum sum of squared dissimilarities
# within blocks
##############################################


# Dissimilarity matrices D(i,j), i = 1,...,n (block), j = 1, ..., p (features)
# Weights w(1), ..., w(p)

# Goal: minimize sum(i=1:n) || sum(j=1:p) w(j) D(i,j) ||_F^2
# with respect to w(1), ..., w(p)


# Inputs:
# D:	array (n,n,p) (n = total number of records, p = number of features)
# block: integer vector of length n (block indicator)

# Output:
# weights w
# matrix D' such that objective = || D' w ||^2


learn_metric <- function(D, block)
{
  nblock <- length(unique(block))
  p <- dim(D)[3]
  Dlist <- vector("list", nblock*p)
  dim(Dlist) <- c(nblock, p)
  for (i in 1:nblock) {
    idx <- which(block == i)
    for (j in 1:p) {
      mat <- D[idx,idx,j]
      Dlist[[i,j]] <- mat[lower.tri(mat)]
    }
  }
  Dmat <- matrix(unlist(Dlist), ncol = p)
  svdDmat <- svd(Dmat)
  w <- svdDmat$v[,1] #w <- svdDmat$v[,p] Question: is this p correct? out of bound
  
  return(list(D = Dmat, w = w))
}


##########
# Variant
##########

# min || D' w ||^2  such that w >= 0 sum(w) = 1


learn_metric2 <- function(D, block)
{
  nblock <- length(unique(block))
  p <- dim(D)[3]
  Dlist <- vector("list", nblock*p)
  dim(Dlist) <- c(nblock, p)
  for (i in 1:nblock) {
    idx <- which(block == i)
    for (j in 1:p) {
      mat <- D[idx,idx,j]
      Dlist[[i,j]] <- mat[lower.tri(mat)]
    }
  }
  Dmat <- matrix(unlist(Dlist), ncol = p)
  
  
  
  # use function solve.QP
  # to minimize ||D'w||^2 = w'DD'w = w'Vw
  V <- crossprod(Dmat)
  Amat <- cbind(rep(1,p), diag(1,p))
  # sum weights = 1, all weights >= 0
  dvec <- rep(0, p)
  bvec <- c(1, rep(0,p))
  if(is_positive_definite(V)){
    sol <- solve.QP(V, dvec, Amat, bvec, meq=1)
    w <- sol$solution
  }else{
    w <- rep(0,p)
  }
  
  return(list(D = Dmat, w = w))
  
}



# Check if the matrix is positive definite using chol()
is_positive_definite <- function(mat) {
  tryCatch({
    chol(mat)
    TRUE # If chol() succeeds, the matrix is positive definite
  }, error=function(e) {
    FALSE # If an error occurs, the matrix is not positive definite
  })
}

