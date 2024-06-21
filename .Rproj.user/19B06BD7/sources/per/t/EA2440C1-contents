
## Graph conversion functions between different representations

# Graph Representation: 
# adjacency (square) matrix, 
# adjacency list,
# and edge list (two-column matrix) 

# Graph G(V, E): neighbors list for each vertex v --> N(v)


#######################################
# Function to convert adjacency matrix
# to adjacency list
#######################################


adjmat2adjlist <- function(x)
{
  stopifnot(is.matrix(x))
  stopifnot(nrow(x) == ncol(x))
  stopifnot(all(x == t(x)))  
  xnames <- rownames(x)
  if (!is.logical(x)) 
    x <- matrix(as.logical(x), nrow(x), ncol(x))
  if (any(diag(x))) diag(x) <- FALSE
  out <- apply(x, 2, which, simplify = FALSE)
  names(out) <- xnames
  out
}


#####################################
# Function to convert adjacency list
# to adjacency matrix
#####################################

adjlist2adjmat <- function(x)
{
  stopifnot(is.list(x))
  u <- unique(unlist(x))
  n <- length(x)
  out <- matrix(0L, n, n)
  rownames(out) <- colnames(out) <- names(x)
  for (i in 1:n) 
  	out[i,x[[i]]] <- out[x[[i]],i] <- 1L
  out	
}



################################
# Function to convert edge list
# to adjacency matrix
################################

edges2adjmat <- function(x) 
{
  stopifnot(is.matrix(x) && ncol(x) == 2)
  stopifnot(is.numeric(x) || is.character(x))
  u <- unique(c(x))
  n <- if (is.numeric(u)) max(u) else length(u)
  out <- matrix(0L, n, n)
  if (is.character(x))
    rownames(out) <- colnames(out) <- u
  out[x] <- 1L
  out <- out + t(out)
  diag(out) <- 0L
  out  
}


#######################################
# Function to convert adjacency matrix
# to edge list
#######################################

adjmat2edges <- function(x) 
{
  stopifnot(is.matrix(x) && nrow(x) == ncol(x))
  which(x != 0 & upper.tri(x), TRUE)
}


################################
# Function to convert edge list 
# to adjacency list
################################

edges2adjlist <- function(x)
{
  stopifnot(is.matrix(x) && ncol(x) == 2)
  stopifnot(is.numeric(x) || is.character(x))
  u <- unique(c(x))
  n <- length(u)
  out <- vector("list", n)
  names(out) <- u
  for (i in 1:n) {
  	idx1 <- (x[,1] == u[i])
  	idx2 <- (x[,2] == u[i])
    vals <- c(x[idx1,2], x[idx2,1])
    out[[i]] <- unique(vals)   
  }
  out
}

#######################################
# Function to convert adjacency list
# to edge list
#######################################

adjlist2edges <- function(x)
{
  stopifnot(is.list(x))
  adjmat <- adjlist2adjmat(x)
  out <- adjmat2edges(adjmat)
  out	
}

