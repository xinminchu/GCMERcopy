#######################################
# Function to coloring a given graph
# with specific method
#######################################

graph_coloring <- function(graph, 
	method = c("lf", "sl", "dsatur", "rlf", "msc", "lmxrlf", "tabu"))
{
stopifnot(is.matrix(graph) || is.list(graph))
if (is.matrix(graph)) {
	test <- is.logical(graph) || all(graph %in% c(0,1))
	graph <- if (test) adjmat2adjlist(graph) else edges2adjlist(graph)
} 
method <- match.arg(method, several.ok = TRUE)
n <- length(graph)
out <- matrix(0, n, length(method))
colnames(out) <- method
for (m in method)
	out[,m] <- switch(m, 
		lf = graph_coloring_greedy(graph, "lf"),
		sl = graph_coloring_greedy(graph, "sl"),
		dsatur = graph_coloring_dsatur(graph),
		rlf = graph_coloring_rlf(graph),
		msc = graph_coloring_msc(graph),
		lmxrlf = graph_coloring_lmxrlf(graph),
		tabu = graph_coloring_tabucol(graph, ...))

out
}		


########################
# Function to perform 
# greedy graph coloring 
########################

graph_coloring_greedy <- function(graph, method = c("lf", "sl")) 
{
  stopifnot(is.list(graph))
  method <- match.arg(method)
  n <- length(graph)
  coloring <- integer(n)
  all_colors <- 1:n
  degree <- sapply(graph, length)
  if (method == "lf") {
    sorted_vertices <- order(degree, decreasing = TRUE)
  } else {
    sorted_vertices <- integer(n)
    for (i in n:1) {
      idx <- which.min(degree)
      sorted_vertices[i] <- idx
      if (degree[idx] > 0)
        degree[graph[[idx]]] <- degree[graph[[idx]]] - 1L
      degree[idx] <- NA
    }
  }
  for (vertex in sorted_vertices) {
    used_colors <- unique(coloring[graph[[vertex]]])
    used_colors <- used_colors[used_colors > 0L]
    coloring[vertex] <- 	if (length(used_colors) > 0) {
      all_colors[-used_colors][1] } else 1L
  }  
  coloring
}


########################################
# Recursive Largest First (RLF) method 
# for graph coloring
########################################

## need function "get_independent_set" in helpers.R

graph_coloring_rlf <- function(graph)
{
  stopifnot(is.list(graph) || is.matrix(graph))
  if (is.list(graph)) {
    graph <- adjlist2adjmat(graph)
  } else if (ncol(graph) == 2) {
    graph <- edges2adjmat(graph)
  }
  n <- NCOL(graph)
  if (n == 1) return(1L)
  color <- integer(n)
  current_color <- 0L
  uncolored <- 1:n 
  while (length(uncolored) > 0) {
    current_color <- current_color + 1L
    # cat("new color", current_color, "\n")
    idx <- get_independent_set(graph[uncolored, uncolored])
    color[uncolored[idx]] <- current_color
    uncolored <- uncolored[-idx]
  }	
  color
}
