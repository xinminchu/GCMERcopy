\name{resolve_entities}
\alias{resolve_entities}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Entity Resolution by Graph Coloring of Binarized Dissimilarity Matrix}

\description{
This function takes a dissimilarity matrix, turns it into an unweighted similarity graph through thresholding, and then applies graph coloring methods to partition the graph into as independent sets as possible. The vertex sets (which are cliques in the complementary dissimilarity graph) form the resolved entities.
}

\usage{
resolve_entities(D, thresholds, method = c("lf", "sl", "rlf"),
	na.action = c("omit", "fail", "pass"))
}

\arguments{
  \item{D}{
(square and symmetric) dissimilarity matrix or vectorized lower-half-matrix
}
  \item{thresholds}{
numerical thresholds for binarizing \code{D}
}
  \item{method}{
graph coloring method(s). Options are "lf", "sl", "dsatur", "rlf", "msc", "lmxrlf", and "tabu"
}
  \item{na.action}{
method for handling missing values
}
}

\details{

The input \code{D} should either be a symmetric matrix (missing values are allowed) or the lower half of one (diagonal excluded). Each row (and the associated column) corresponds to an instance, record, or object.

The values in \code{thresholds} should be in the range of \code{D} for sensible results, although this is not checked.

The default value for \code{method} includes the fast heuristics \code{lf}, \code{sl}, \code{rlf} for graph coloring that seem in practice to use as many or fewer colors as the slower methods \code{dsatur}, \code{msc}, \code{lmxrlf}, and \code{tabu}.

For \code{na.action}, the value \code{fail} throws an error if any missing values are found in \code{D}. The default value \code{omit} looks for NA's. If it finds any, it sort row (or equivalently column) indices by decreasing number of NA counts and removes rows/columns one by one until no NA remains. This is helpful when \code{D} has a small number of rows and columns that contain only NAs
and all other rows/columns are complete. Finally, the value \code{pass} replace all NAs in \code{D} by \code{max(thresholds)+1} so that the corresponding graph vertices never have an edge between them (and therefore do not belong to the same entity). Other strategies for handling missing values should be carried out manually on \code{D} before calling \code{resolve_entities}.



}

\value{
A list with components
\item{ents}{vector or matrix of resolved entities as integer vectors. Dimensions \code{ncol(D) x length(thresholds)}} (when \code{D} is a matrix)
\item{thresholds}{sorted input argument \code{thresholds}}
}

\references{
%% ~put references to the literature/web site here ~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{graph_coloring}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (D, thresholds, method = c("lf", "sl", "rlf"), na.action = c("omit",
    "fail", "pass"))
{
    stopifnot(is.vector(D) || is.matrix(D))
    if (is.matrix(D) && NCOL(D) == 1)
        dim(D) <- NULL
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
    }
    else if (na.action == "pass") {
        test <- is.na(D) | is.nan(D)
        if (any(test))
            D[test] <- max(thresholds)
    }
    method <- match.arg(arg = method, choices = c("lf", "sl",
        "dsatur", "rlf", "msc", "lmxrlf", "tabu"), several.ok = TRUE)
    if (length(D) == 0)
        return(list(ents = NULL, thresholds = thresholds))
    if (length(D) == 1)
        return(list(ents = rep(1, nthres), thresholds = thresholds))
    if (is.matrix(D)) {
        nrecs <- ncol(D)
    }
    else {
        nrecs <- (1 + sqrt(1 + 8 * length(D)))/2
        if (abs(nrecs - round(nrecs)) > 1e-08)
            stop("Argument 'D' does not match the number of elements",
                "of the lower triangular part of a matrix", "(not of the form n * (n-1) / 2).")
    }
    ents <- matrix(0L, nrecs, nthres)
    if (is.vector(D)) {
        mask <- lower.tri(diag(nrecs))
        Dvec <- D
        D <- matrix(0, nrecs, nrecs)
        D[mask] <- Dvec
        D <- D + t(D)
    }
    if (na.action == "omit") {
        test <- is.na(D) | is.nan(D)
        diag(test) <- FALSE
        if (any(test)) {
            degree <- nrecs - colSums(test)
            ord <- order(degree)
            test <- test[ord, ord]
            count <- 2
            while (count <= nrecs && any(test[count:nrecs, count:nrecs])) count <- count +
                1
            keep0 <- ord[count:nrecs]
            D <- D[keep0, keep0]
        }
        else {
            keep0 <- 1:nrecs
        }
    }
    graph <- NULL
    for (j in 1:nthres) {
        print(thresholds[j])
        coloring <- rep(NA, nrecs)
        keep <- keep0
        graph.prev <- graph
        graph <- (D <= thresholds[j])
        diag(graph) <- FALSE
        if (identical(graph.prev, graph)) {
            ents[, j] <- ents[, j - 1]
            next
        }
        degree <- colSums(graph)
        idx0 <- which(degree == 0)
        ncolors <- length(idx0)
        if (ncolors > 0)
            coloring[keep[idx0]] <- 1:ncolors
        idx1 <- which(degree == 1)
        if (length(idx1) > 0) {
            g <- graph[idx1, idx1]
            g[lower.tri(g)] <- FALSE
            idx2 <- which(g, TRUE)
            npairs <- NROW(idx2)
            if (npairs > 0) {
                idx1 <- idx1[c(idx2)]
                coloring[keep[idx1]] <- rep((ncolors + 1):(ncolors +
                  npairs), 2)
                ncolors <- ncolors + npairs
            }
        }
        idx <- c(idx0, idx1)
        if (length(idx) > 0) {
            graph <- graph[-idx, -idx]
            keep <- keep[-idx]
        }
        if (length(keep) == 0) {
            ents[, j] <- coloring
            next
        }
        comps <- find_connected_components(graph)
        ncomps <- max(comps)
        for (k in 1:ncomps) {
            idx <- which(comps == k)
            compk <- graph[idx, idx]
            if (all(compk[lower.tri(compk)])) {
                ncolors <- ncolors + 1L
                coloring[keep[idx]] <- ncolors
                next
            }
            compk <- !compk
            diag(compk) <- FALSE
            entsk <- graph_coloring(compk, method)
            ncolorsk <- apply(entsk, 2, get_chromatic)
            entsk <- entsk[, which.min(ncolorsk)]
            ncolorsk <- min(ncolorsk)
            entsk <- match(entsk, unique(entsk))
            coloring[keep[idx]] <- ncolors + entsk
            ncolors <- ncolors + ncolorsk
        }
        ents[, j] <- coloring
    }
    list(ents = ents, thresholds = thresholds)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
