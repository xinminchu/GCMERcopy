
#######################################
# Function to find entities
# from a given dissimilarity matrix
# with thresholds
#######################################

resolve_entities <- function(D, thresholds, 
	method = c("lf", "sl", "rlf"),
	na.action = c("omit", "fail", "pass"))
{
## Check inputs
stopifnot(is.vector(D) || is.matrix(D))
if (is.matrix(D) && NCOL(D) == 1) dim(D) <- NULL
if (is.matrix(D)) {
	stopifnot(nrow(D) == ncol(D))
	stopifnot(all((D == t(D)) | (is.na(D) & is.na(t(D)))))
}
stopifnot(is.numeric(thresholds))
thresolds <- sort(thresholds)
nthres <- length(thresholds)
na.action <- match.arg(na.action)
if (na.action == "fail") {
	stopifnot(all(!is.na(D) & is.nan(D)))
} else if (na.action == "pass") {
	test <- is.na(D) | is.nan(D)
	if (any(test)) D[test] <- max(thresholds)
}
method <- match.arg(arg = method, 
	choices = c("lf", "sl", "dsatur", "rlf", "msc", "lmxrlf", "tabu"),
	several.ok = TRUE)

## Handle trivial cases
if (length(D) == 0) 
	return(list(ents = NULL, thresholds = thresholds))
if (length(D) == 1) 
	return(list(ents = rep(1, nthres), thresholds = thresholds))

## Determine the number of records
if (is.matrix(D)) {
	nrecs <- ncol(D)
} else { # case: vectorized half-matrix
	nrecs <- (1 + sqrt(1 + 8 * length(D))) / 2
	if (abs(nrecs - round(nrecs)) > 1e-8)
		stop("Argument 'D' does not match the number of elements",
			"of the lower triangular part of a matrix", 
			"(not of the form n * (n-1) / 2).")
}

## Set up output
ents <- matrix(0L, nrecs, nthres) # estimated entities 

## Reconstruct dissimilarity matrix
if (is.vector(D)) {
	mask <- lower.tri(diag(nrecs))
	Dvec <- D
	D <- matrix(0, nrecs, nrecs)
	D[mask] <- Dvec
	D <- D + t(D)
}

## Handle missing values
if (na.action == "omit") {
	test <- is.na(D) | is.nan(D)
	diag(test) <- FALSE
	 if (any(test)) {
		degree <- nrecs - colSums(test)
		ord <- order(degree)
		test <- test[ord,ord]
		count <- 2
		while (count <= nrecs && any(test[count:nrecs,count:nrecs]))
			count <- count + 1
		keep0 <- ord[count:nrecs]
		D <- D[keep0,keep0]
	} else {
		keep0 <- 1:nrecs
	}
}	

## Initialize graph
graph <- NULL

## Loop over thresholds 
for (j in 1:nthres) {
	
	## Coloring vector
	coloring <- rep(NA, nrecs)
	
	## Book keeping
	keep <- keep0
	graph.prev <- graph
	
	## Threshold dissimilarity matrix 
	graph <- (D <= thresholds[j])
	diag(graph) <- FALSE
	
	## Stop if current graph identical to previous one 
	if (identical(graph.prev, graph)) {
		ents[,j] <- ents[,j-1]
		next
	}

	## Determine degree of each vertex in original graph
	degree <- colSums(graph)
	
	## Identify vertices of degree 0 (= cliques of size 1)
	idx0 <- which(degree == 0)
	ncolors <- length(idx0)
	if (ncolors > 0) 
		coloring[keep[idx0]] <- 1:ncolors
			
	## Determine cliques of size 2 (not strictly needed)
	idx1 <- which(degree == 1)
	if (length(idx1) > 0) {
		g <- graph[idx1,idx1]
		g[lower.tri(g)] <- FALSE	
		idx2 <- which(g, TRUE)
		npairs <- NROW(idx2) 
		if (npairs > 0) {
			idx1 <- idx1[c(idx2)]
			coloring[keep[idx1]] <- 
				rep((ncolors+1):(ncolors+npairs), 2)
			ncolors <- ncolors + npairs
		}
	}
	
	## Reduce graph
	idx <- c(idx0,idx1)
	if (length(idx) > 0) {
		graph <- graph[-idx, -idx]
		keep <- keep[-idx]
	}
	
	## Case: empty graph
	if (length(keep) == 0) {
		ents[,j] <- coloring
		next
	}
		
	## Identify connected components
	comps <- find_connected_components(graph)
	ncomps <- max(comps)		
	
	## Graph coloring on connected components
	for (k in 1:ncomps) {
		idx <- which(comps == k)
		compk <- graph[idx,idx]
		
		## Check for trivial case: clique
		if (all(compk[lower.tri(compk)])) { 
			ncolors <- ncolors + 1L
			coloring[keep[idx]] <- ncolors
			next
		} 
		
		## Take the complement of the component
		compk <- !compk
		diag(compk) <- FALSE
		entsk <- graph_coloring(compk, method)

		## Retain a coloring that use the least colors
		ncolorsk <- apply(entsk, 2, get_chromatic)
		entsk <- entsk[,which.min(ncolorsk)]
		ncolorsk <- min(ncolorsk)
		
		## Relabel colors to 1, 2, ... 
		entsk <- match(entsk, unique(entsk))
		
		## Update clique cover
		coloring[keep[idx]] <- ncolors + entsk
		ncolors <- ncolors + ncolorsk			
	}
	ents[,j] <- coloring				
}

list(ents = ents, thresholds = thresholds)

}
