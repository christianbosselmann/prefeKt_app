# pipeline: msa

# packages
library("librarian")
librarian::shelf(tidyverse,
                 data.table,
                 msa,
                 seqinr,
                 ape,
                 bios2mds,
                 stringr,
                 quiet = TRUE)

set.seed(seed)

# get paths and fasta files
path <- list.files(path = "fasta", 
                   pattern = "\\.fasta$", 
                   full.names = TRUE)
fasta <- readAAStringSet(path)

# fix names
names <- tolower(str_extract(basename(path), "^([^.]+)"))
names(fasta) <- names

# run alignment
aligned <- msa(fasta, method = "Muscle", cluster = "upgmamax")

# get distance values and draw phylogenetic tree
aligned_conv <- msaConvert(aligned, type = "seqinr::alignment")
distance <- dist.alignment(aligned_conv, "similarity", gap = 1)
tree <- nj(distance)

tiff("tree.tiff", units="in", width=6, height=5, res=300)
plot(tree)
dev.off()

# export distance matrix for later use as similarity matrix
d <- data.frame(as.matrix(distance))
write_csv(d, "distancematrix.csv")

# canonical id
aligned_cw <- msaConvert(aligned, "bios2mds::align")

cid <- data.frame()

for (i in seq_along(names)) {
  cid <- aligned_cw[names[i]] %>%
    unlist(use.names = FALSE) %>%
    rbind(cid)
}

colnames(cid) <- seq_along(cid)
cid <- t(cid)
colnames(cid) <- rev(names)
cid <- as.data.table(cid)

for (i in colnames(cid)) {
  cid[[i]] <- cumsum(cid[[i]] %in% LETTERS)
}

cid$cid <- rownames(cid)

write_csv(cid, "features/cid.csv")

# clean up non-base packages
# otherwise there are some strange Bioconductor-related errors later
# lapply(names(sessionInfo()$otherPkgs), function(pkgs)
#   detach(
#     paste0('package:', pkgs),
#     character.only = TRUE,
#     unload = TRUE,
#     force = TRUE
#   ))
