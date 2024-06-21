\name{mqcpp}
\alias{mqcpp}

\title{
Solve the Maximal Quasi-Clique Partition Problem (MQCPP)
}
\description{
Formulate the MQCPP as an integer linear program (ILP) and solve it with the GUROBI commerical solver or the R package \code{highs} (depending on availability) 
}

\usage{
mqcpp(edges, gamma, n = NULL, ub = NULL)
}

\arguments{
  \item{edges}{
Two-column matrix of graph edges. The matrix can contain either integers or character strings
}
  \item{gamma}{
Required minimum edge density of quasi-cliques
}
  \item{n}{
Number of vertices in graph. Defaults to the number of unique elements in \code{edges}
}
  \item{ub}{
Upper bound on the number of quasi-cliques. Default value calculated as in Melo et al (2023)
}
}

\details{
A standard ILP formulation of the MQCPP with \eqn{(|V|+1+|E|) U} variables 
is used, where \eqn{|V|} is the number of vertices in the graph, \eqn{|E|} the number of edges, and \eqn{U} the upper bound on the number of quasi-cliques to look for. See Section 3.1 of Melo et al (2022) for details. 

To solve the ILP, the function first tries to use the GUROBI mathematical programming solver by loading its R package interface \code{gurobi}. (GUROBI can be downloaded at \url{https://www.gurobi.com} and used freely under academic license.) If this fails, the function tries to use the R package \code{highs} to solve the ILP. An error message is thrown if neither package is available. 

Note: the \code{gurobi} solver is in general much faster than \code{highs}.
}

\value{
List with components 
\item{qc}{integer vector of quasi-clique membership}
\item{nqc}{number of quasi-cliques found}
\item{dens}{edge density of each quasi-clique}
}

\references{
Melo, Ribeiro, and Riveaux (2022). The minimum quasi-clique partitioning problem: Complexity, formulations, and a computational study. \emph{Information Sciences} 612, 655–674.
}


\examples{
## Example graph of Melo et al (2023)
edges <- matrix(c(1,2,1,5,2,3,3,4,3,5,4,5,5,6,5,8,5,9,6,7,7,8,8,9), 
	ncol = 2, byrow = TRUE)

## Visualize graph with package 'igraph' if installed
test <- require(igraph)
if (test) {
	g <- make_undirected_graph(t(edges))
	plot(g)
}

## Solve MQCPP with edge density threshold 0.6
solution <- mqcpp(edges, .6)

## Display solution 
solution$nqc  # number of quasi-cliques found
solution$dens # edge density of each quasi-clique
solution$qc   # quasi-clique membership vector 
table(solution$qc) # quasi-clique sizes

## Visualize solution
if (test) {
	qc.list <- split(1:vcount(g), solution$qc) # list of quasi-cliques
	plot(g, mark.groups = qc.list, main = "Quasi-Clique Partition of Graph")
}
}
