# GCMER
## Graph Coloring Method for Entity Resolution  

*This R package provides general functions for graph coloring as well as specific functions for entity resolution, aka data deduplication.*  

**To install the package:** 
```
library(devtools)
devtools::install_github("https://github.com/ddegras/GCMER")
```

### Main package functionalities
- **Graph coloring methods:** greedy, recursive largest first (RLF), DSATUR, maximum cardinal search (MCS), lmXRLF, tabu, ...
- **Entity resolution via graph coloring** 
- **Maximal Quasi-clique Partitioning Problem (MQCPP)** 
- **Clustering agreement measures** based on matched pairs, set overlap, mutual information, ...

The graph coloring functions are mostly taken from <https://github.com/saurfang/graphcoloring> which itself is based on <https://github.com/brrcrites/graph-coloring>.

The entity resolution (ER) problem consists in mapping the records of a database to the (unknown) real-world persons or entities they relate to. (It is assumed that the data do not contain uniquely identifying information.) Given a dissimilarity table $D$ between records and a numerical threshold $\gamma$, one can create an unweighted, undirected graph $G$ where two vertices (i.e. records) $i$ and $j$  are connected by an edge if and only their dissimilarity $D_{ij}$ is below the threshold:  $\gamma$.  One can then perform ER by solving the minimum clique partition problem on $G$ (1 clique = 1 entity), that is, find cliques $C_1, ..., C_k$ that partition $G$ with $k$ as small as possible. The clique structure guarantees that the records within an entity are homogeneous (all dissimilarities bounded above by $\gamma$). Computationally, finding the minimum clique partition of $G$ is equivalent to solving the graph coloring problem on the complement graph $\bar{G}$, i.e. finding independent sets $S_1, ..., S_k$ that partition $\bar{G}$ with $k$ as small as possible.  

Details on MQCPP and its integer linear programming (ILP) solutions can be found in Melo et al (2022) (https://doi.org/10.1016/j.ins.2022.08.073)

The clustering agreement measures are described in Wagner and Wagner (2007) (https://publikationen.bibliothek.kit.edu/1000011477/812079)

