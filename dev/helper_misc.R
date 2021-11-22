# PREFEKT
# PREdicting the Functional Effects of Kv muTations
# deployment helper functions

# requires data_backup.csv (pre-processed training data)

# create vlookup table from pipeline output (master.R) by rbinding str_ objects
vlookup <- rbind(data.frame(id = "Shaker", str_kcnas),
                 data.frame(id = "KCNA1", str_kcna1),
                 data.frame(id = "KCNA2", str_kcna2),
                 data.frame(id = "KCNA4", str_kcna4),
                 data.frame(id = "KCNA5", str_kcna5),
                 data.frame(id = "KCNB1", str_kcnb1),
                 data.frame(id = "KCNC1", str_kcnc1),
                 data.frame(id = "KCNC2", str_kcnc2),
                 data.frame(id = "KCNC3", str_kcnc3),
                 data.frame(id = "KCND2", str_kcnd2),
                 data.frame(id = "KCND3", str_kcnd3),
                 data.frame(id = "KCNQ1", str_kcnq1),
                 data.frame(id = "KCNQ2", str_kcnq2),
                 data.frame(id = "KCNQ3", str_kcnq3),
                 data.frame(id = "KCNQ4", str_kcnq4),
                 data.frame(id = "KCNQ5", str_kcnq5),
                 data.frame(id = "KCNH1", str_kcnh1),
                 data.frame(id = "KCNH2", str_kcnh2),
                 data.frame(id = "KCNH5", str_kcnh5))

# Notes:
# recipe.rds comes from a modified preprocessing pipeline script, see helper_preprocessing.R
# aa_ objects come from the feature extraction pipeline and represent aa encoding schemes, substitution matrices etc
# Values for C, alpha and sigma are determined by cv in the development branch (model_mtl.R)

# in case of changes to original data set:
# re-run development branch and model selection entirely
# then update recipe.rds and preprocessed data set (saved as training_data.csv) in deployment branch