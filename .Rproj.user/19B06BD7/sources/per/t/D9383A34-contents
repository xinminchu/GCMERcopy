
##########
# Wrapper 
##########

clustering_agreement <- function(x, y = NULL, 
  method = c("chi2", "rand", "adj_rand", 
  "fowlkes_mallow", "mirkin", "jaccard", 
  "tpr", "fpr", "F_measure", "meila_heckerman", 
  "max_match", "van_dongen", "mutual_info")) {
  	
  method <- match.arg(method, several.ok = TRUE)
  tab <- if (is.null(y)) x else table(x, y) # contingency table
  out <- numeric(length(method))
  names(out) <- method
  for (m in method) 
    out[m] <- switch(method,
      chi2 = chi2(tab), 
      rand = rand(tab), 
      adj_rand = adj_rand(tab), 
      fowlkes_mallow = fowlkes_mallow(tab),
      mirkin = mirkin(tab), 
      jaccard = jaccard(tab), 
      tpr = tpr(tab), 
      fpr = fpr(tab), 
      F_measure = F_measure(tab),
      meila_heckerman = meila_heckerman(tab),
      max_match = max_match(tab), 
      van_dongen = van_dongen(tab), 
      mutual_info = mutual_info(tab)
    )
  out
}


#####################################
### Measures based on Counting Pairs
#####################################

#########################
# Chi squared Coefficient
#########################

# Formula: chi(C1, C2) = sum_i=1:k sum_j=1:l (m.ij - E.ij)^2 / E.ij
# E.ij = |C1i||C2j| / n

# Inputs 'x' and 'y' should be (integer or character) vectors
# indicating cluster membership
# OR 'x' can be a contingency table with 'y' being NULL

chi2 <- function(x, y = NULL) {
  tab <- if (is.null(y)) x else table(x, y) # contingency table
  n <- sum(tab) # number of elements
  rsum <- rowSums(tab)
  csum <- colSums(tab)
  E <- tcrossprod(rsum, csum) / n # expected counts
  sum((m-E)^2 / E)
}



######################
# General Rand index
######################

# Formula: R(C1, C2) = 2(n11+n00) / n(n-1)

rand <- function(x, y = NULL) {
  tab <- if (is.null(y)) x else table(x, y)
  n <- sum(tab)
  rsum <- rowSums(tab)
  csum <- colSums(tab)
  1 - (sum(rsum^2) + sum(csum^2) - 2*sum(tab^2)) / (n*(n-1))
}



######################
# Adjusted rand index
######################

# Formula: R_adj(C1, C2) = sum_i=1^K sum_j=1^L choose(tab.ij, 2) - t3
# t1 = sum choose(|C1|, 2), t2 = sum choose(|C2|, 2)
# t3 = 2*t1*t2 / n(n-1)


adj_rand <- function(x, y = NULL) {
  tab <- if (is.null(y)) x else table(x, y)
  n <- sum(tab)
  rsum <- rowSums(tab)
  csum <- colSums(tab)
  t1 <- sum(choose(rsum, 2))
  t2 <- sum(choose(csum, 2))
  t3 <- t1 * t2 / choose(n, 2)
  num <- sum(choose(tab, 2)) - t3
  den <- (t1 + t2) / 2 - t3
  num / den
}




#######################
# Fowlkes-Mallow Index
#######################

# Formula: FM(C1, C2) = n11 / sqrt((n11+n10)(n11+n01))

fowlkes_mallow <- function(x, y = NULL) {
  tab <- if (is.null(y)) x else table(x, y)
  n <- sum(tab)
  n11 <- sum(choose(tab,2)) 
  rsum <- rowSums(tab)
  csum <- colSums(tab)
  n01 <- sum(choose(csum, 2)) - sum(choose(tab, 2))
  n10 <- sum(choose(rsum, 2)) - sum(choose(tab, 2))
  n11 / sqrt((n11 + n10) * (n11 + n01))
}



#######################
# Mirkin Metric
#######################

# Formula: M(C1, C2) = sum |C1.i|^2 + sum |C2.j|^2 - 2 sum sum m.ij^2
# = 2(n01+n10) = n(n-1) (1-R(C1, C2))

mirkin <- function(x, y = NULL) {
  tab <- if (is.null(y)) x else table(x, y)
  rsum <- rowSums(tab)
  csum <- colSums(tab)
  n01 <- sum(choose(csum, 2)) - sum(choose(tab, 2))
  n10 <- sum(choose(rsum, 2)) - sum(choose(tab, 2))
  2 * (n10 + n01)
}



################
# Jaccard Index
################

# Formula: J(C1, C2) = n11 / (n11 + n10 + n01)

jaccard <- function(x, y = NULL) {
  tab <- if (is.null(y)) x else table(x, y)
  n11 <- sum(choose(tab,2)) # sum_{i,j} C(m_ij, 2)
  rsum <- rowSums(tab)
  csum <- colSums(tab)
  n01 <- sum(choose(csum, 2)) - sum(choose(tab, 2))
  n10 <- sum(choose(rsum, 2)) - sum(choose(tab, 2))
  n11 / (n11 + n10 + n01)
}



########################
# Partition Difference
########################

# Formula: PD(C1, C2) = n00

# partition.diff <- function(x, y = NULL) {
  # tab  <- if (is.null(y)) x else table(x, y)
  # n    <- sum(tab)
  # n11  <- sum(choose(tab, 2)) 
  # rsum <- rowSums(tab)
  # csum <- colSums(tab)
  # n01  <- sum(choose(csum, 2)) - n11
  # n10  <- sum(choose(rsum, 2)) - n11
  # n00  <- choose(n, 2) - (n01 + n10 + n11)
  # n00
# }



#####################
# True positive rate 
#####################

# 'y' (or columns of 'x' if 'y' is NULL) 
# should refer to the true partition

tpr <- function(x, y = NULL) {
  tab <- if (is.null(y)) x else table(x, y)
  csum <- colSums(tab)
  tp <- sum(choose(tab, 2)) # true positive (n11)
  p  <- sum(choose(csum, 2)) # positive (n01 + n11)
  tp / p
}

######################
# False positive rate 
######################

# 'y' (or columns of 'x' if 'y' is NULL) 
# should refer to the true partition

fpr <- function(x, y = NULL) {
  tab <- if (is.null(y)) x else table(x, y)
  rsum <- rowSums(tab)
  csum <- colSums(tab)
  fp <- sum(choose(rsum, 2)) - sum(choose(tab, 2)) # false positive
  n  <- choose(sum(tab), 2) - sum(choose(csum, 2)) # negative
  fp / n
}


###################################
### Measures based on Set Overlaps
###################################

############
# F measure
############

# F(C1, C2) = F(C2) = (1/n) sum_1:K ni max_1:L F(C1i, C2j)
# F(C1i, C2j) = 2 * r.ij * p.ij / (r.ij + p.ij) = 2 m / (|C1i| + |C2j|)

# 'y' (or columns of 'x' if 'y' is NULL) 
# should refer to the true partition

F_measure <- function(x, y = NULL) {
  tab <- if (is.null(y)) x else table(x, y)
  n <- sum(tab)
  rsum <- rowSums(tab)
  csum <- colSums(tab)
  nr <- nrow(tab)
  nc <- ncol(tab)
  den <- outer(rsum, csum, "+")
  Fmat <- 2 * tab / den
  sum(rsum * apply(Fmat, 1, max)) / n
}



######################
# Meila-Heckerman
######################

# paper: Meila (2001): An experimental comparison of model-based clustering methods
# Expectation-Maximization (EM) algorithm (Dempster, Laird, Rubin, 1977)
# Classification EM (CEM) algorithm (Celeux, Govaert, 1992)
# model-based agglomerative clustering (AC) (e.g. Banfield, Raftery, 1993)


# Formula: MH(C1, optC) = (1/n) sum_(i=1:k) max_j m.ij

# 'y' (or columns of 'x' if 'y' is NULL) 
# should refer to the true partition

meila_heckerman <- function(x, y = NULL) {
  tab <- if (is.null(y)) x else table(x, y)
  n <- sum(tab)
  rmax <- apply(tab, 1, max)
  sum(rmax) / n
}


########################
# Maximum-Match measure
########################


max_match <- function(x, y = NULL) {
  tab <- if (is.null(y)) x else table(x, y)
  n <- sum(tab)
  out <- 0
  nr <- nrow(tab)
  nc <- ncol(tab)
  for (i in 1:min(nr,nc)) {
    idx <- arrayInd(which.max(tab), c(nr, nc))
    out <- out + tab[idx[1], idx[2]]
    tab[idx[1],] <- NA
    tab[,idx[2]] <- NA
  }
  out / n
}


######################
# Van Dongen-measure
######################

# Formula: D(C1, C2) = 2n - sum_i=1:k max_j m.ij - sum_j=1:l max_i m.ij

van_dongen <- function(x, y) {
  tab <- if (is.null(y)) x else table(x, y)
  n <- sum(tab)
  rmax <- apply(tab, 1, max)
  cmax <- apply(tab, 2, max)
  2 * n - sum(rmax) - sum(cmax)
}



#########################################
### Measures based on Mutual Information
#########################################

# Entropy associated with clustering
# Formula: H(C) = - sum_i=1:k P(i) log_2 P(i)

# assuming all elements of X have the same probability of being picked
# choosing an element of X at random p = 1/n
# the probability that this element is in cluster Ci \in C is
# P(i) = |Ci| / n


# Mutual information
# Formula: I(C1, C2) = sum_i=1:k sum_j=1:l P(i,j) log2 P(i,j) / (P(i)*P(j))
# P(i,j) = |C1i \cap C2j| / n

mutual_info <- function(x, y = NULL) {
  tab <- if (is.null(y)) x else table(x, y)
  n <- sum(tab)  
  px <- rowSums(tab) / n 
  py <- colSums(tab) / n
  if (any(px == 0) || any(py == 0)) {
  	tab <- tab[px > 0, py > 0, drop = FALSE]
  	px <- px[px > 0]
  	py <- py[py > 0]
  }
  Hx <- - sum(px * log2(px)) 
  Hy <- - sum(py * log2(py))
  pxy <- tab / n 
  I <- sum(pxy * log2(pxy / tcrossprod(px, py))) # mutual information (MI)
  SG <- I / sqrt(Hx * Hy) #normalized MI by Strehl & Ghosh (2002)
  FJ <- 2 * I / (Hx + Hy) #normalized MI by Fred & Jain (2003)
  VI <- Hx + Hy - 2 * I # Variation of Information by Meila (2003)
  c(MI = I, G = SG, FJ = FJ, VI = VI)
}
