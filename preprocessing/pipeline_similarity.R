# PREFEKT
# PREdicting the Functional Effects of Kv muTations
# this script leverages taxonomy-based information on task-relatedness by 
# using the msa-distance matrix to obtaining similarity matrices for the MTL kernel
# one similarity matrix for each value of alpha (see below)
# cf doi.org/10.1007/978-3-642-12683-3_34
# requires the distance matrix provided by pipeline_msa.R

# packages
library("librarian")
librarian::shelf(tidyverse,
                 data.table,
                 openxlsx,
                 Matrix,
                 matrixcalc,
                 ggcorrplot,
                 quiet = TRUE)

# read data
d <- read_csv("distancematrix.csv")
rownames(d) <- colnames(d)

# remember that dist.alignment = sqrt(1 - identity). Let's get the pairwise distance.
d = d^2

# convert to similarity (gamma) by the transformation gamma = alpha − d / dmax
# where alpha ≥ 1 is a hyperparameter to control baseline similarity between tasks
# and dmax is the maximal distance. alpha is set in the master script.
for (a in alpha) {
  gamma <- a - (d / max(d))
  gamma <- as.matrix(gamma)
  
  # matrices need to be positive semi-definite to yield a valid kernel. 
  if (is.positive.semi.definite(as.matrix(gamma)) == FALSE) {
    gamma <- nearPD(as.matrix(gamma), 
                    keepDiag = TRUE)
    
  }
  
  # fix names
  gamma <- data.frame(as.matrix(gamma$mat))
  rownames(gamma) <- rownames(d)
  colnames(gamma) <- colnames(d)
  
  # final check
  if (is.positive.semi.definite(as.matrix(gamma)) == FALSE) {
    stop("Error, similarity matrix is not PDM. Kernel will be invalid.")
    
  }
  
  # export
  write_csv(gamma, paste('similaritymatrix_a', a, '.csv', sep = ''))
  
  # assign alpha value as matrix id
  assign(paste("gamma", a, sep = "_a"), gamma)   
}

# visualization
p <- ggcorrplot(cor(gamma_a1))
plot(p)

tiff("similarity.tiff", units="in", width=8, height=8, res=300)
plot(p)
dev.off()

# clean up
rm(list = grep("^gamma_", ls(), value = TRUE, invert = TRUE)) 