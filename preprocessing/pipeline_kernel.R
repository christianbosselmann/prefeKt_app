# PREFEKT
# PREdicting the Functional Effects of Kv muTations
# get kernel matrix for each set of hyperparameters and scale by similarity 
# matrix to obtain precomputed kernel matrices for later use in the model.
# requires the similarity matrices provided by pipeline_similarity.R

# packages
library("librarian")
librarian::shelf(tidyverse,
                 data.table,
                 Matrix,
                 matrixcalc,
                 kernlab,
                 quiet = TRUE)

set.seed(seed)

# get data
data <- read_csv("data_backup.csv")

# fix the group names
data <- data %>%
  mutate(gene = tolower(gene)) %>%
  mutate(gene = replace(gene, gene == "shaker", "kcnas"))

# set up data without label and group information for later
f <- data %>%
  select(-y) %>%
  select(-gene)

# get standard rbf kernel matrices (Kb) as baseline comparison for later
Kb <- list()
for (s in 1:length(sigma)) {
  Kb[[s]] <- list()
  Kb[[s]] <- as.data.frame(
    kernelMatrix(
      rbfdot(sigma = sigma[s]), 
      as.matrix(f)
    ))
}

# save Kb (list of baseline rbf kernel matrices for each sigma)
saveRDS(Kb, "kernelmatrices_base.rds")

# get similarity matrices from main folder, as output by pipeline_similarity.R
path <- list.files(pattern = "^similaritymatrix_", recursive = FALSE)
path <- str_sort(path, numeric = TRUE) # maintain order of alpha vector

sim_matrices <- lapply(path, read_csv)
names(sim_matrices) <- gsub(path, pattern=".csv$", replacement="")

# fixing rownames
for (m in 1:length(sim_matrices)) {
  rownames(sim_matrices[[m]]) <- colnames(sim_matrices[[m]])
}

# reshape to long format, yielding the pairwise similarity for all combinations of tasks
# loop through all values of alpha, the baseline similarity, then store as list
sim_match <- vector(mode = "list", length = length(path))

for (a in 1:length(path)) {
  sim_match[[a]] <- sim_matrices[[a]] %>%
    rownames_to_column() %>%
    pivot_longer(cols = -c(1)) %>%
    dplyr::rename(k = rowname, l = name)
}

names(sim_match) <- names(sim_matrices)

# constructing Kt, the task similarity matrix
Kt <- data.frame(diag(length(t)))

# get vector of task identity. store for later, used for task-specific model assessment later on..
t <- data$gene
save(t, file = "t_vec.rda")

# nested loop
Km <- list()

for (a in 1:length(sim_match)) {
  Km[[a]] <- list()
  for (s in 1:length(sigma)) {
    Km[[a]][[s]] <- list()
    
    # loop through all Kt(i,j) and enter the corresponding task similarity value s(ti,tj)
    for (i in 1:length(t)) {
      for (j in 1:length(t)) { 
        
        Kt[i,j] # empty index of i,j
        t[i] # task i
        t[j] # task j
        
        Kt[i,j] <- sim_match[[a]] %>%
          filter(k == t[i] & l == t[j]) %>%
          select(value) # similarity(ti,tj) for a given similarity matrix a
      }
    }
    
    # check if Kt is positive semi-definite, otherwise it will not yield a valid kernel
    if (is.positive.semi.definite(as.matrix(Kt)) == FALSE) {
      stop("Error, matrix is not PDM. Kernel will be invalid.")
    }
    
    # constructing Kf, the kernel matrix of the dataset (i.e. similarity between observations)
    f <- data %>%
      select(-y) %>%
      select(-gene) # drop group information and label
    
    Kf <- as.data.frame(
      kernelMatrix(
        rbfdot(sigma = sigma[s]), 
        as.matrix(f)
      ))
    
    Km[[a]][[s]] <- as.matrix(Kf) * as.matrix(Kt) # element-wise multiplication of Kf and Kt to obtain the final kernel metrix Km
  }
}

# save Km (nested list of kernel matrices by alpha*sigma)
saveRDS(Km, "kernelmatrices_mtl.rds")
