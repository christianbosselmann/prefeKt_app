# PREFEKT
# PREdicting the Functional Effects of Kv muTations

# packages
library("librarian")
librarian::shelf(tidyverse,
                 tidymodels,
                 data.table,
                 openxlsx,
                 yardstick,
                 caret,
                 bestNormalize,
                 kernlab,
                 e1071,
                 data.table,
                 magic,
                 quiet = TRUE)

# get preprocessed training data set
train <- read_csv("training_data.csv")
y <- as.factor(train$y)

# get similarity matrix
sim_matrices <- read_csv("similaritymatrix.csv")
rownames(sim_matrices) <- colnames(sim_matrices)

sim_match <- sim_matrices %>%
  rownames_to_column() %>%
  pivot_longer(cols = -c(1)) %>%
  rename(k = rowname, l = name)

# get pre-processing recipe
load("recipe.rds")

# get feature lookup tables
vlookup <- read_csv("vlookup.csv") %>% rename(gene = id)
aa_blosum62 <- read_csv("aa_blosum62.csv")
aa_blomap <- read_csv("aa_blomap.csv")
aa_braun <- read_csv("aa_braun.csv")
aa_grantham <- read_csv("aa_grantham.csv")
aa_hphob <- read_csv("aa_hphob.csv")
cid <- read_csv("cid.csv")

# get test data
# input <- read.xlsx("input.xlsx")
input <- df_in

test <- input %>%
  merge(aa_blosum62, by = c("aa1", "aa2")) %>%
  merge(aa_blomap, by = c("aa1", "aa2")) %>%
  merge(aa_braun, by = c("aa1", "aa2")) %>%
  merge(aa_grantham, by = c("aa1", "aa2")) %>%
  merge(aa_hphob, by = c("aa1", "aa2")) %>%
  merge(vlookup, by = c("pos", "gene"), all = FALSE) %>%
  merge(cid)

# apply prep recipe
test <- bake(data_rec, test)

# calculate rbf kernel matrix Kb
df_train <- train %>%
  select(-y) %>%
  select(-gene)

t_train <- train %>%
  mutate(gene = tolower(gene)) %>%
  mutate(gene = replace(gene, gene == "shaker", "kcnas")) %>%
  select(gene) %>%
  as_vector() # training set task vector

df_test <- test %>%
  select(-gene)

t_test <- test %>%
  mutate(gene = tolower(gene)) %>%
  mutate(gene = replace(gene, gene == "shaker", "kcnas")) %>%
  select(gene) %>%
  as_vector() # test set task vector

df_all <- rbind(data.frame(id = "train", df_train),
                data.frame(id = "test", df_test))

t_all <- c(t_train, t_test)

Kb <- as.data.frame(
  kernelMatrix(
    rbfdot(sigma = 0.01), 
    as.matrix(df_all[,-1])
  ))

# calculate task similarity kernel matrix Kt
Kt <- matrix(data = 0, nrow = length(t_all), ncol = length(t_all),
             dimnames = list(as.vector(t_all), as.vector(t_all)))

for (i in 1:length(t_all)) {
  
  for (j in 1:length(t_all)) { 
    
    Kt[i,j] <- sim_match$value[sim_match$k == t_all[[i]] & sim_match$l == t_all[[j]]] # similarity(ti,tj) for a given similarity matrix a
  
    }
}

# calculate MTL kernel Km
Km <- as.matrix(Kb) * as.matrix(Kt)

# set up class weights
weights <- c("GOF" = length(train$y) / sum(train$y == "GOF"),
             "LOF" = length(train$y) / sum(train$y == "LOF"), 
             "Neutral" = length(train$y) / sum(train$y == "Neutral"))

# set up train and test indices
train_indices <- 1:length(t_train)
test_indices <- (length(t_train)+1):(length(t_train)+length(t_test))

# training SVM
model <- e1071::svm(
  x = Km[c(-test_indices), c(-test_indices), drop = FALSE],
  y = y,
  cost = 1,
  probability = TRUE,
  class.weights = weights
)

# predicting test set
test <- as.kernelMatrix(Km[test_indices, -test_indices, drop = FALSE])
prediction <- predict(model, test) %>%
  as_tibble() %>%
  cbind(attr(predict(model, test, probability = TRUE), "probabilities")) %>%
  cbind(attr(predict(model, test, decision.value = T), "decision.values")) %>%
  rename(prediction = value, 
         prob.LOF = LOF, prob.Neutral = Neutral, prob.GOF = GOF, 
         `conf.LOF_Neutral` = `LOF/Neutral`, `conf.LOF_GOF` = `LOF/GOF`, `conf.Neutral_GOF` = `Neutral/GOF`)

# print output
out <- cbind(input, prediction)
write.xlsx(out, "output.xlsx")

# verbose output for app
verb_out <- paste("This variant is predicted to be:", out$prediction, sep = " ")
