# pipeline: kcnc1

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
str_kcnc1 <- data.table()

str_kcnc1 <-
  readr::read_delim("features/kcnc1/query.consurf.grades", delim = "\t", skip = 12, col_names = TRUE) %>%
  slice(2:512) %>%
  select(SCORE) %>%
  rename(str_consurf = SCORE) %>%
  add_column(str_kcnc1)

str_kcnc1 <-
  readr::read_delim("features/kcnc1/query.isis", delim = " ", skip = 39, col_names = FALSE) %>%
  select(X3) %>%
  rename(str_isis = X3) %>%
  add_column(str_kcnc1)

str_kcnc1 <-
  readr::read_delim("features/kcnc1/query.mdisorder", delim = "\t", skip = 1, col_names = FALSE) %>%
  slice(1:511) %>%
  select(X3, X5, X7) %>%
  rename(str_nors = X3, str_profbval = X5, str_ucon = X7) %>%
  add_column(str_kcnc1)

str_kcnc1 <-
  readr::read_delim("features/kcnc1/query.phdRdb", delim = "\t", skip = 42, col_names = TRUE) %>%
  slice(2:512) %>%
  select(OtH, OtL, PRHL, PiTo) %>%
  rename(str_helix = OtH, str_loop = OtL, str_prhl = PRHL, str_pito = PiTo) %>%
  add_column(str_kcnc1)

str_kcnc1 <-
  readr::read_delim("features/kcnc1/query.prof1Rdb", delim = "\t", skip = 78, col_names = TRUE) %>%
  select(PACC, PREL, Pbie) %>%
  rename(str_asa = PACC, str_rsa = PREL, str_pbie = Pbie) %>%
  add_column(str_kcnc1)

# features/structure/iupred2a
# documentation: DOI 10.1093/nar/gky384, DOI 10.1002/cpbi.99
# method: input canonical uniprot isoform, select default long disorder, download html
str_kcnc1 <-
  readr::read_delim("features/kcnc1/iupred.html", delim = "\t", skip = 4, col_names = TRUE) %>%
  slice(1:511) %>%
  select("IUPRED SCORE", "ANCHOR SCORE") %>%
  rename(str_iupred = "IUPRED SCORE", str_anchor = "ANCHOR SCORE") %>%
  add_column(str_kcnc1)

# features/structure/uniprot
# method: topology and motifs from uniprot. done manually for now.
# note: make sure to spell-check and use same terms across all tasks.
str_kcnc1$str_top <- NA

str_kcnc1$str_top[1:190] <- "cytoplasmic"
str_kcnc1$str_top[191:209] <- "S1"
str_kcnc1$str_top[210:247] <- "extracellular"
str_kcnc1$str_top[248:267] <- "S2"
str_kcnc1$str_top[268:276] <- "cytoplasmic"
str_kcnc1$str_top[277:295] <- "S3"
str_kcnc1$str_top[296:308] <- "extracellular"
str_kcnc1$str_top[309:331] <- "S4"
str_kcnc1$str_top[332:344] <- "cytoplasmic"
str_kcnc1$str_top[345:366] <- "S5"
str_kcnc1$str_top[367:414] <- "extracellular"
str_kcnc1$str_top[415:436] <- "S6"
str_kcnc1$str_top[437:511] <- "cytoplasmic"

str_kcnc1$str_top[is.na(str_kcnc1$str_top)] <- "none"

# features/structure/motif
# method: additional motifs via expert annotation. done manually for now.
# gating charges ref: 10.1085/jgp.114.5.723
# s4-s5 linker ref: 10.1126/science.1116270
str_kcnc1$str_exp <- NA

str_kcnc1$str_exp[400:405] <- "selectivityfilter"
str_kcnc1$str_exp[332:344] <- "S4S5linker"
str_kcnc1$str_exp[311] <- "gating"
str_kcnc1$str_exp[314] <- "gating"
str_kcnc1$str_exp[317] <- "gating"
str_kcnc1$str_exp[320] <- "gating"

str_kcnc1$str_exp[is.na(str_kcnc1$str_exp)] <- "none"

# features/structure/perviewer
# method: download gene family .txt file and rename accordingly
# documentation: DOI 10.1101/gr.252601.119
str_kcnc1 <- readr::read_delim("features/kcnc1/query.per.txt", delim = "\t", skip = 1, col_names = FALSE) %>%
  filter(!grepl('-', X3)) %>%
  select(X6) %>%
  rename(str_paraz = X6) %>%
  slice(1:511) %>%
  add_column(str_kcnc1)

# features/structure/netsurfp
# method: run netfsurp2.0 webservice (https://services.healthtech.dtu.dk/service.php?NetSurfP-2.0) and export all files as csv. rename accordingly.
# documentation: DOI 10.1002/prot.25674
str_kcnc1 <-
  readr::read_delim("features/kcnc1/netsurfp.csv", delim = ",", skip = 0, col_names = TRUE) %>%
  select("rsa", "asa", "q3", "q8") %>%
  rename(str_np_rsa = "rsa", str_np_asa = "asa", str_np_q3 = "q3", str_np_q8 = "q8") %>%
  add_column(str_kcnc1)

# features/id/ecdf
# method: absolute relative position of mutation, divided by total transcript length
str_kcnc1$pos <- seq_len(nrow(str_kcnc1))
str_kcnc1$ecdf <- ecdf(str_kcnc1$pos)(str_kcnc1$pos)

# now merge all prepared tables with input
kcnc1 <- read.xlsx(xlsxFile = "input/kcnc1.xlsx", sheet = 1, colNames = TRUE)

kcnc1 <- kcnc1 %>%
  merge(aa_blosum62, by = c("aa1", "aa2")) %>%
  merge(aa_blomap, by = c("aa1", "aa2")) %>%
  merge(aa_braun, by = c("aa1", "aa2")) %>%
  merge(aa_grantham, by = c("aa1", "aa2")) %>%
  merge(aa_hphob, by = c("aa1", "aa2")) %>%
  merge(str_kcnc1, by.x = "pos", by.y = "pos", all = FALSE)

# sanity check
if (sum(is.na(kcnc1)) > 0){
  stop("missing data, check pipeline.")
}

# done