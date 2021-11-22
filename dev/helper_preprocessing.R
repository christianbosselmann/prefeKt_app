# PREFEKT
# PREdicting the Functional Effects of Kv muTations
# pre-processing pipeline

# packages
library("librarian")
librarian::shelf(tidyverse,
                 data.table,
                 hablar,
                 janitor,
                 fastDummies,
                 tidymodels,
                 quiet = TRUE)

# remove redundant observations and find conflicts
# conflicts are labeled "unclear" and removed alongside other unclear variants
dupes <- data %>%
  unique() %>%
  get_dupes(aa1, pos, aa2, gene)

data <- anti_join(data, dupes) %>%
  filter(y != "Unclear") %>%
  unique()

# preprocessing
data <- data %>%
  retype()

data <- data %>%
  select(-one_of("pos", 
                 "aa1", 
                 "aa2"))

data_rec <- recipe(y ~ ., data) %>%
  step_dummy(c("str_pbie", 
               "str_prhl", 
               "str_pito", 
               "str_top", 
               "str_exp", 
               "str_np_q3", 
               "str_np_q8"), one_hot = TRUE) %>%
  step_normalize(all_numeric())

data_rec <- prep(data_rec)

data <- bake(data_rec, NULL)

save(data_rec, file = "recipe.rds") # save preprocessing recipe

# export to data backup
write_csv(data, "data_backup.csv")
