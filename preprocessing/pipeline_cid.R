# pipeline: canonical id (cid)
# based on pipeline_msa to generate cid.csv (look-up table of family alignment)

# packages
library("librarian")
librarian::shelf(tidyverse,
                 data.table,
                 quiet = TRUE)

set.seed(7)

# load cid.csv and fix names
cid <- read_csv("features/cid.csv")

names <- c("KCNQ5", "KCNQ4", "KCNQ3", "KCNQ2", "KCNQ1", 
           "KCNH5", "KCNH2", "KCNH1",
           "KCNC3", "KCNC2", "KCNC1", 
           "KCND3", "KCND2",
           "KCNB1", 
           "Shaker", "KCNA5", "KCNA4", "KCNA2", "KCNA1", 
           "cid")

colnames(cid) <- names

# do cid
cid <- cid %>%
  pivot_longer(!cid, names_to = "gene", values_to = "pos") %>%
  group_by(gene) %>%
  arrange(pos) %>%
  distinct(pos, .keep_all = TRUE)

# merge
data$pos <- as.numeric(data$pos)
data <- merge(data, cid)
