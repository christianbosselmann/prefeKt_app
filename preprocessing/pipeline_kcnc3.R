# pipeline: kcnc3

# packages
library("librarian")
librarian::shelf(tidyverse,
                 data.table,
                 openxlsx,
                 quiet = TRUE)

# define single amino acid code vector
aa_aa <- c("G","A","L", "M","F","W","K","Q","E","S","P","V","I","C","Y","H","R","N","D","T")

# get all unique missense variants, i.e. remove duplice aa1-aa2 pairs
aa_all_aa <- crossing(aa_aa, aa_aa, .name_repair = ~ c("aa1", "aa2")) %>%
  filter(aa1 != aa2)

# features/sequence/blomap
# standard amino acid encoding
# documentation: https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.71.2602&rep=rep1&type=pdf
aa_blomap_raw <- readr::read_tsv("features/aa_blomap.tsv")

aa_blomap <- aa_all_aa %>%
  merge(aa_blomap_raw, by.x = "aa1", "aa") %>%
  rename(aa_blomap_aa1_dim1 = dim1,
         aa_blomap_aa1_dim2 = dim2,
         aa_blomap_aa1_dim3 = dim3,
         aa_blomap_aa1_dim4 = dim4,
         aa_blomap_aa1_dim5 = dim5) %>%
  merge(aa_blomap_raw, by.x = "aa2", "aa") %>%
  rename(aa_blomap_aa2_dim1 = dim1,
         aa_blomap_aa2_dim2 = dim2,
         aa_blomap_aa2_dim3 = dim3,
         aa_blomap_aa2_dim4 = dim4,
         aa_blomap_aa2_dim5 = dim5)

# features/sequence/braun
# standard amino acid encoding
# documentation: DOI 10.1007/s00894-001-0058-5
aa_braun_raw <- readr::read_tsv("features/aa_braun.tsv") 

aa_braun <- aa_all_aa %>%
  merge(aa_braun_raw, by.x = "aa1", "aa") %>%
  rename(aa_braun_aa1_E1 = E1,
         aa_braun_aa1_E2 = E2,
         aa_braun_aa1_E3 = E3,
         aa_braun_aa1_E4 = E4,
         aa_braun_aa1_E5 = E5) %>%
  merge(aa_braun_raw, by.x = "aa2", "aa") %>%
  rename(aa_braun_aa2_E1 = E1,
         aa_braun_aa2_E2 = E2,
         aa_braun_aa2_E3 = E3,
         aa_braun_aa2_E4 = E4,
         aa_braun_aa2_E5 = E5)

# features/sequence/hydrophobicity
# documentation: DOI 10.1186/s40659-016-0092-5
# method: standard pca to reduce dimensionality across aa pairs
aa_hphob_raw <- readr::read_tsv("features/aa_hphob.tsv")

aa_hphob <- aa_all_aa %>%
  merge(aa_hphob_raw, by.x = "aa1", "aa") %>%
  merge(aa_hphob_raw, by.x = "aa2", "aa")

aa_hphob_pca <- prcomp(aa_hphob[3:198], scale = TRUE)

aa_hphob <- aa_hphob_pca$x[,1:3] %>%
  as_tibble() %>%
  add_column(aa_hphob) %>%
  select(aa1, aa2, PC1, PC2, PC3) %>%
  rename(aa_hphob_pca1 = PC1, aa_hphob_pca2 = PC2, aa_hphob_pca3 = PC3)

# features/sequence/blosum62
# standard substitution matrix
aa_blosum62_raw <- readr::read_tsv("features/aa_blosum62.tsv") %>%
  rename(aa1 = FIRST, aa2 = SECOND, aa_blosum62 = SCORE)

aa_blosum62 <- aa_all_aa %>%
  merge(aa_blosum62_raw, by = c("aa1", "aa2"))

# features/sequence/grantham
# standard substitution matrix
aa_grantham_raw <- readr::read_tsv("features/aa_grantham.tsv") %>%
  rename(aa1 = FIRST, aa2 = SECOND, aa_grantham = SCORE)

aa_grantham <- aa_all_aa %>%
  merge(aa_grantham_raw, by = c("aa1", "aa2"))

# features/structure/predictprotein
# documentation: https://rostlab.org/owiki/index.php/PredictProtein_-_Documentation
# method: input canonical uniprot isoform, download all results
str_kcnc3 <- data.table()

str_kcnc3 <-
  readr::read_delim("features/kcnc3/query.consurf.grades", delim = "\t", skip = 12, col_names = TRUE) %>%
  slice(2:758) %>%
  select(SCORE) %>%
  rename(str_consurf = SCORE) %>%
  add_column(str_kcnc3)

str_kcnc3 <-
  readr::read_delim("features/kcnc3/query.isis", delim = " ", skip = 57, col_names = FALSE) %>%
  select(X3) %>%
  rename(str_isis = X3) %>%
  add_column(str_kcnc3)

str_kcnc3 <-
  readr::read_delim("features/kcnc3/query.mdisorder", delim = "\t", skip = 1, col_names = FALSE) %>%
  slice(1:757) %>%
  select(X3, X5, X7) %>%
  rename(str_nors = X3, str_profbval = X5, str_ucon = X7) %>%
  add_column(str_kcnc3)

str_kcnc3 <-
  readr::read_delim("features/kcnc3/query.phdRdb", delim = "\t", skip = 42, col_names = TRUE) %>%
  slice(2:758) %>%
  select(OtH, OtL, PRHL, PiTo) %>%
  rename(str_helix = OtH, str_loop = OtL, str_prhl = PRHL, str_pito = PiTo) %>%
  add_column(str_kcnc3)

str_kcnc3 <-
  readr::read_delim("features/kcnc3/query.prof1Rdb", delim = "\t", skip = 78, col_names = TRUE) %>%
  select(PACC, PREL, Pbie) %>%
  rename(str_asa = PACC, str_rsa = PREL, str_pbie = Pbie) %>%
  add_column(str_kcnc3)

# features/structure/iupred2a
# documentation: DOI 10.1093/nar/gky384, DOI 10.1002/cpbi.99
# method: input canonical uniprot isoform, select default long disorder, download html
str_kcnc3 <-
  readr::read_delim("features/kcnc3/iupred.html", delim = "\t", skip = 4, col_names = TRUE) %>%
  slice(1:757) %>%
  select("IUPRED SCORE", "ANCHOR SCORE") %>%
  rename(str_iupred = "IUPRED SCORE", str_anchor = "ANCHOR SCORE") %>%
  add_column(str_kcnc3)

# features/structure/uniprot
# method: topology and motifs from uniprot. done manually for now.
# note: make sure to spell-check and use same terms across all tasks.
str_kcnc3$str_top <- NA

str_kcnc3$str_top[1:290] <- "cytoplasmic"
str_kcnc3$str_top[291:309] <- "S1"
str_kcnc3$str_top[310:350] <- "extracellular"
str_kcnc3$str_top[351:370] <- "S2"
str_kcnc3$str_top[371:379] <- "cytoplasmic"
str_kcnc3$str_top[380:398] <- "S3"
str_kcnc3$str_top[399:411] <- "extracellular"
str_kcnc3$str_top[412:434] <- "S4"
str_kcnc3$str_top[435:447] <- "cytoplasmic"
str_kcnc3$str_top[448:469] <- "S5"
str_kcnc3$str_top[470:517] <- "extracellular"
str_kcnc3$str_top[518:539] <- "S6"
str_kcnc3$str_top[540:757] <- "cytoplasmic"

str_kcnc3$str_top[is.na(str_kcnc3$str_top)] <- "none"

# features/structure/motif
# method: additional motifs via expert annotation. done manually for now.
# gating charges ref: doi.org/10.1038/ng1758
str_kcnc3$str_exp <- NA

str_kcnc3$str_exp[503:508] <- "selectivityfilter"
str_kcnc3$str_exp[435:447] <- "S4S5linker"

str_kcnc3$str_exp[414] <- "gating"
str_kcnc3$str_exp[417] <- "gating"
str_kcnc3$str_exp[420] <- "gating"
str_kcnc3$str_exp[423] <- "gating"

str_kcnc3$str_exp[is.na(str_kcnc3$str_exp)] <- "none"

# features/structure/perviewer
# method: download gene family .txt file and rename accordingly
# documentation: DOI 10.1101/gr.252601.119
# problem: no paralog family for KCNC3, hence no parazscore. use dummy value.
str_kcnc3$str_paraz <- 1

# features/structure/netsurfp
# method: run netfsurp2.0 webservice (https://services.healthtech.dtu.dk/service.php?NetSurfP-2.0) and export all files as csv. rename accordingly.
# documentation: DOI 10.1002/prot.25674
str_kcnc3 <-
  readr::read_delim("features/kcnc3/netsurfp.csv", delim = ",", skip = 0, col_names = TRUE) %>%
  select("rsa", "asa", "q3", "q8") %>%
  rename(str_np_rsa = "rsa", str_np_asa = "asa", str_np_q3 = "q3", str_np_q8 = "q8") %>%
  add_column(str_kcnc3)

# features/id/ecdf
# method: absolute relative position of mutation, divided by total transcript length
str_kcnc3$pos <- seq_len(nrow(str_kcnc3))
str_kcnc3$ecdf <- ecdf(str_kcnc3$pos)(str_kcnc3$pos)

# now merge all prepared tables with input
kcnc3 <- read.xlsx(xlsxFile = "input/kcnc3.xlsx", sheet = 1, colNames = TRUE)

kcnc3 <- kcnc3 %>%
  merge(aa_blosum62, by = c("aa1", "aa2")) %>%
  merge(aa_blomap, by = c("aa1", "aa2")) %>%
  merge(aa_braun, by = c("aa1", "aa2")) %>%
  merge(aa_grantham, by = c("aa1", "aa2")) %>%
  merge(aa_hphob, by = c("aa1", "aa2")) %>%
  merge(str_kcnc3, by.x = "pos", by.y = "pos", all = FALSE)

# sanity check
if (sum(is.na(kcnc3)) > 0){
  stop("missing data, check pipeline.")
}

# done