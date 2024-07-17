library(dplyr)

args = commandArgs(trailingOnly=TRUE)
dataset <- args[1]

#Read the data
df <- read.csv(file = paste0("../data/",dataset,"/evo_likelihoods_all.csv"))



#List the sample
samples <- list.dirs(path = paste0("../data/",dataset,"/VDJ"), full.names = FALSE, recursive = FALSE)

cor_df <- c()
for(sample in samples){
  #Subset the dataframe per sample
  df_sample <- df[df$original_sample_id == sample,]
  
  Ablang_ESM = c()
  Ablang_ProtBERT = c()
  Ablang_Sapiens = c()
  ProtBERT_ESM = c()
  Sapiens_ESM = c()
  Sapiens_ProtBERT = c()
  
  ## Full VDJ
  Ablang_ESM <- c(Ablang_ESM, cor.test(x = df_sample$evo_likelihood_ablang_full_VDJ,
                                       y = df_sample$evo_likelihood_esm_full_VDJ,
                                       method = "pearson")$estimate)
  #signifanctly correlated
  Ablang_ProtBERT <- c(Ablang_ProtBERT, cor.test(x = df_sample$evo_likelihood_ablang_full_VDJ,
                                                 y = df_sample$evo_likelihood_protbert_full_VDJ,
                                                 method = "pearson")$estimate)
  #signifanctly correlated
  Ablang_Sapiens <- c(Ablang_Sapiens, cor.test(x = df_sample$evo_likelihood_ablang_full_VDJ,
                                               y = df_sample$evo_likelihood_sapiens_full_VDJ,
                                               method = "pearson")$estimate)
  ProtBERT_ESM <- c(ProtBERT_ESM, cor.test(x = df_sample$evo_likelihood_protbert_full_VDJ,
                                           y = df_sample$evo_likelihood_esm_full_VDJ,
                                           method = "pearson")$estimate)
  Sapiens_ESM <- c(Sapiens_ESM, cor.test(x = df_sample$evo_likelihood_sapiens_full_VDJ,
                                         y = df_sample$evo_likelihood_esm_full_VDJ,
                                         method = "pearson")$estimate)
  #signifanctly correlated
  Sapiens_ProtBERT <- c(Sapiens_ProtBERT, cor.test(x = df_sample$evo_likelihood_sapiens_full_VDJ,
                                                   y = df_sample$evo_likelihood_protbert_full_VDJ,
                                                   method = "pearson")$estimate)
  
  ## CDR3 from VDJ
  Ablang_ESM <- c(Ablang_ESM, cor.test(x = df_sample$evo_likelihood_ablang_cdr3_from_VDJ,
                                       y = df_sample$evo_likelihood_esm_cdr3_from_VDJ,
                                       method = "pearson")$estimate)
  Ablang_ProtBERT <- c(Ablang_ProtBERT, cor.test(x = df_sample$evo_likelihood_ablang_cdr3_from_VDJ,
                                                 y = df_sample$evo_likelihood_protbert_cdr3_from_VDJ,
                                                 method = "pearson")$estimate)
  #signifanctly correlated
  Ablang_Sapiens <- c(Ablang_Sapiens, cor.test(x = df_sample$evo_likelihood_ablang_cdr3_from_VDJ,
                                               y = df_sample$evo_likelihood_sapiens_cdr3_from_VDJ,
                                               method = "pearson")$estimate)
  ProtBERT_ESM <- c(ProtBERT_ESM, cor.test(x = df_sample$evo_likelihood_protbert_cdr3_from_VDJ,
                                           y = df_sample$evo_likelihood_esm_cdr3_from_VDJ,
                                           method = "pearson")$estimate)
  Sapiens_ESM <- c(Sapiens_ESM, cor.test(x = df_sample$evo_likelihood_sapiens_cdr3_from_VDJ,
                                         y = df_sample$evo_likelihood_esm_cdr3_from_VDJ,
                                         method = "pearson")$estimate)
  Sapiens_ProtBERT <- c(Sapiens_ProtBERT, cor.test(x = df_sample$evo_likelihood_sapiens_cdr3_from_VDJ,
                                                   y = df_sample$evo_likelihood_protbert_cdr3_from_VDJ,
                                                   method = "pearson")$estimate)
  
  ## CDR3 only
  Ablang_ESM <- c(Ablang_ESM, cor.test(x = df_sample$evo_likelihood_ablang_cdr3_only,
                                       y = df_sample$evo_likelihood_esm_cdr3_only,
                                       method = "pearson")$estimate)
  Ablang_ProtBERT <- c(Ablang_ProtBERT, cor.test(x = df_sample$evo_likelihood_ablang_cdr3_only,
                                                 y = df_sample$evo_likelihood_protbert_cdr3_only,
                                                 method = "pearson")$estimate)
  Ablang_Sapiens <- c(Ablang_Sapiens, cor.test(x = df_sample$evo_likelihood_ablang_cdr3_only,
                                               y = df_sample$evo_likelihood_sapiens_cdr3_only,
                                               method = "pearson")$estimate)
  ProtBERT_ESM <- c(ProtBERT_ESM, cor.test(x = df_sample$evo_likelihood_protbert_cdr3_only,
                                           y = df_sample$evo_likelihood_esm_cdr3_only,
                                           method = "pearson")$estimate)
  Sapiens_ESM <- c(Sapiens_ESM, cor.test(x = df_sample$evo_likelihood_sapiens_cdr3_only,
                                         y = df_sample$evo_likelihood_esm_cdr3_only,
                                         method = "pearson")$estimate)
  Sapiens_ProtBERT <- c(Sapiens_ProtBERT, cor.test(x = df_sample$evo_likelihood_sapiens_cdr3_only,
                                                   y = df_sample$evo_likelihood_protbert_cdr3_only,
                                                   method = "pearson")$estimate)
  
  #Combine into a dataframe
  cor_df_sample <- data.frame(Ablang_ESM, Ablang_ProtBERT, Ablang_Sapiens, ProtBERT_ESM, Sapiens_ESM, Sapiens_ProtBERT,
                              sample, source = c("Full VDJ", "CDR3 from VDJ", "CDR3 only"))
  #Append to the main dataframe
  cor_df <- rbind(cor_df, cor_df_sample)

}

write.csv(cor_df, file = paste0("../data/",dataset,"/PLMCorrelation.csv"), row.names = F)

