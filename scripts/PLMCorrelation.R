library(dplyr)

args = commandArgs(trailingOnly=TRUE)
dataset <- args[1]

#Read the data
#Read the data
load(paste0("../data/",dataset,"/VDJ_PLL_",dataset,".RData")) #vdj_likelihood
vdj <- vdj_likelihood

#List the sample
samples <- list.dirs(path = paste0("../data/",dataset,"/VDJ"), full.names = FALSE, recursive = FALSE)

cor_df <- c()
for(sample in samples){
  #Subset the dataframe per sample
  df_sample <- vdj[vdj$sample_id == sample,]
  
  Ablang1_ESM1b = c()
  Ablang1_ProtBERT = c()
  Ablang1_Sapiens = c()
  ProtBERT_ESM1b = c()
  Sapiens_ESM1b = c()
  Sapiens_ProtBERT = c()
  Ablang2_ESM1b = c()
  Ablang2_ProtBERT = c()
  Ablang2_Sapiens = c()
  Ablang2_Ablang1 = c()
  Ablang2_ESMc = c()
  ESMc_ESM1b = c()
  ESMc_ProtBERT = c()
  ESMc_Sapiens = c()
  ESMc_Ablang1 = c()
  
  ## Full VDJ
  Ablang1_ESM1b <- c(Ablang1_ESM1b, cor.test(x = c(df_sample$HC_evo_likelihood_ablang1_full_VDJ,
                                                   df_sample$LC_evo_likelihood_ablang1_full_VDJ),
                                             y = c(df_sample$HC_evo_likelihood_esm1b_full_VDJ,
                                                   df_sample$LC_evo_likelihood_esm1b_full_VDJ),
                                       method = "pearson")$estimate)
  Ablang1_ProtBERT <- c(Ablang1_ProtBERT, cor.test(x = c(df_sample$HC_evo_likelihood_ablang1_full_VDJ,
                                                         df_sample$LC_evo_likelihood_ablang1_full_VDJ),
                                                 y = c(df_sample$HC_evo_likelihood_protbert_full_VDJ,
                                                       df_sample$LC_evo_likelihood_protbert_full_VDJ),
                                                 method = "pearson")$estimate)
  Ablang1_Sapiens <- c(Ablang1_Sapiens, cor.test(x = c(df_sample$HC_evo_likelihood_ablang1_full_VDJ,
                                                       df_sample$LC_evo_likelihood_ablang1_full_VDJ),
                                               y = c(df_sample$HC_evo_likelihood_sapiens_full_VDJ,
                                                     df_sample$LC_evo_likelihood_sapiens_full_VDJ),
                                               method = "pearson")$estimate)
  ProtBERT_ESM1b <- c(ProtBERT_ESM1b, cor.test(x = c(df_sample$HC_evo_likelihood_protbert_full_VDJ,
                                                   df_sample$LC_evo_likelihood_protbert_full_VDJ),
                                           y = c(df_sample$HC_evo_likelihood_esm1b_full_VDJ,
                                                 df_sample$LC_evo_likelihood_esm1b_full_VDJ),
                                           method = "pearson")$estimate)
  Sapiens_ESM1b <- c(Sapiens_ESM1b, cor.test(x = c(df_sample$HC_evo_likelihood_sapiens_full_VDJ,
                                                   df_sample$LC_evo_likelihood_sapiens_full_VDJ),
                                         y = c(df_sample$HC_evo_likelihood_esm1b_full_VDJ,
                                               df_sample$LC_evo_likelihood_esm1b_full_VDJ),
                                         method = "pearson")$estimate)
  Sapiens_ProtBERT <- c(Sapiens_ProtBERT, cor.test(x = c(df_sample$HC_evo_likelihood_sapiens_full_VDJ,
                                                         df_sample$LC_evo_likelihood_sapiens_full_VDJ),
                                                   y = c(df_sample$HC_evo_likelihood_protbert_full_VDJ,
                                                         df_sample$LC_evo_likelihood_protbert_full_VDJ),
                                                   method = "pearson")$estimate)
  Ablang2_ESM1b <- c(Ablang2_ESM1b, cor.test(x = c(df_sample$HC_evo_likelihood_ablang2_full_VDJ,
                                                   df_sample$LC_evo_likelihood_ablang2_full_VDJ),
                                           y = c(df_sample$HC_evo_likelihood_esm1b_full_VDJ,
                                                 df_sample$LC_evo_likelihood_esm1b_full_VDJ),
                                           method = "pearson")$estimate)
  Ablang2_ProtBERT <- c(Ablang2_ProtBERT, cor.test(x = c(df_sample$HC_evo_likelihood_ablang2_full_VDJ,
                                                         df_sample$LC_evo_likelihood_ablang2_full_VDJ),
                                                 y = c(df_sample$HC_evo_likelihood_protbert_full_VDJ,
                                                       df_sample$LC_evo_likelihood_protbert_full_VDJ),
                                                 method = "pearson")$estimate)
  Ablang2_Sapiens <- c(Ablang2_Sapiens, cor.test(x = c(df_sample$HC_evo_likelihood_ablang2_full_VDJ,
                                                       df_sample$LC_evo_likelihood_ablang2_full_VDJ),
                                                   y = c(df_sample$HC_evo_likelihood_sapiens_full_VDJ,
                                                         df_sample$LC_evo_likelihood_sapiens_full_VDJ),
                                                   method = "pearson")$estimate)
  Ablang2_Ablang1 <- c(Ablang2_Ablang1, cor.test(x = c(df_sample$HC_evo_likelihood_ablang2_full_VDJ,
                                                       df_sample$LC_evo_likelihood_ablang2_full_VDJ),
                                                   y = c(df_sample$HC_evo_likelihood_ablang1_full_VDJ,
                                                         df_sample$LC_evo_likelihood_ablang1_full_VDJ),
                                                   method = "pearson")$estimate)
  Ablang2_ESMc <- c(Ablang2_ESMc, cor.test(x = c(df_sample$HC_evo_likelihood_ablang2_full_VDJ,
                                                 df_sample$LC_evo_likelihood_ablang2_full_VDJ),
                                         y = c(df_sample$HC_evo_likelihood_esmc_full_VDJ,
                                               df_sample$LC_evo_likelihood_esmc_full_VDJ),
                                         method = "pearson")$estimate)
  ESMc_ESM1b <- c(ESMc_ESM1b, cor.test(x = c(df_sample$HC_evo_likelihood_esmc_full_VDJ,
                                             df_sample$LC_evo_likelihood_esmc_full_VDJ),
                                         y = c(df_sample$HC_evo_likelihood_esm1b_full_VDJ,
                                               df_sample$LC_evo_likelihood_esm1b_full_VDJ),
                                         method = "pearson")$estimate)
  ESMc_ProtBERT <- c(ESMc_ProtBERT, cor.test(x = c(df_sample$HC_evo_likelihood_esmc_full_VDJ,
                                                   df_sample$LC_evo_likelihood_esmc_full_VDJ),
                                           y = c(df_sample$HC_evo_likelihood_protbert_full_VDJ,
                                                 df_sample$LC_evo_likelihood_protbert_full_VDJ),
                                           method = "pearson")$estimate)
  ESMc_Sapiens <- c(ESMc_Sapiens, cor.test(x = c(df_sample$HC_evo_likelihood_esmc_full_VDJ,
                                                   df_sample$LC_evo_likelihood_esmc_full_VDJ),
                                               y = c(df_sample$HC_evo_likelihood_sapiens_full_VDJ,
                                                     df_sample$LC_evo_likelihood_sapiens_full_VDJ),
                                               method = "pearson")$estimate)
  ESMc_Ablang1 <- c(ESMc_Ablang1, cor.test(x = c(df_sample$HC_evo_likelihood_esmc_full_VDJ,
                                                   df_sample$LC_evo_likelihood_esmc_full_VDJ),
                                               y = c(df_sample$HC_evo_likelihood_ablang1_full_VDJ,
                                                     df_sample$LC_evo_likelihood_ablang1_full_VDJ),
                                               method = "pearson")$estimate)
  
  ## CDR3 from VDJ
  Ablang1_ESM1b <- c(Ablang1_ESM1b, cor.test(x = c(df_sample$HC_evo_likelihood_ablang1_cdr3_from_VDJ,
                                                   df_sample$LC_evo_likelihood_ablang1_cdr3_from_VDJ),
                                             y = c(df_sample$HC_evo_likelihood_esm1b_cdr3_from_VDJ,
                                                   df_sample$LC_evo_likelihood_esm1b_cdr3_from_VDJ),
                                             method = "pearson")$estimate)
  Ablang1_ProtBERT <- c(Ablang1_ProtBERT, cor.test(x = c(df_sample$HC_evo_likelihood_ablang1_cdr3_from_VDJ,
                                                         df_sample$LC_evo_likelihood_ablang1_cdr3_from_VDJ),
                                                 y = c(df_sample$HC_evo_likelihood_protbert_cdr3_from_VDJ,
                                                       df_sample$LC_evo_likelihood_protbert_cdr3_from_VDJ),
                                                 method = "pearson")$estimate)
  Ablang1_Sapiens <- c(Ablang1_Sapiens, cor.test(x = c(df_sample$HC_evo_likelihood_ablang1_cdr3_from_VDJ,
                                                       df_sample$LC_evo_likelihood_ablang1_cdr3_from_VDJ),
                                               y = c(df_sample$HC_evo_likelihood_sapiens_cdr3_from_VDJ,
                                                     df_sample$LC_evo_likelihood_sapiens_cdr3_from_VDJ),
                                               method = "pearson")$estimate)
  ProtBERT_ESM1b <- c(ProtBERT_ESM1b, cor.test(x = c(df_sample$HC_evo_likelihood_protbert_cdr3_from_VDJ,
                                                   df_sample$LC_evo_likelihood_protbert_cdr3_from_VDJ),
                                           y = c(df_sample$HC_evo_likelihood_esm1b_cdr3_from_VDJ,
                                                 df_sample$LC_evo_likelihood_esm1b_cdr3_from_VDJ),
                                           method = "pearson")$estimate)
  Sapiens_ESM1b <- c(Sapiens_ESM1b, cor.test(x = c(df_sample$HC_evo_likelihood_sapiens_cdr3_from_VDJ,
                                                   df_sample$LC_evo_likelihood_sapiens_cdr3_from_VDJ),
                                         y = c(df_sample$HC_evo_likelihood_esm1b_cdr3_from_VDJ,
                                               df_sample$LC_evo_likelihood_esm1b_cdr3_from_VDJ),
                                         method = "pearson")$estimate)
  Sapiens_ProtBERT <- c(Sapiens_ProtBERT, cor.test(x = c(df_sample$HC_evo_likelihood_sapiens_cdr3_from_VDJ,
                                                         df_sample$LC_evo_likelihood_sapiens_cdr3_from_VDJ),
                                                   y = c(df_sample$HC_evo_likelihood_protbert_cdr3_from_VDJ,
                                                         df_sample$LC_evo_likelihood_protbert_cdr3_from_VDJ),
                                                   method = "pearson")$estimate)
  Ablang2_ESM1b <- c(Ablang2_ESM1b, cor.test(x = c(df_sample$HC_evo_likelihood_ablang2_cdr3_from_VDJ,
                                                   df_sample$LC_evo_likelihood_ablang2_cdr3_from_VDJ),
                                           y = c(df_sample$HC_evo_likelihood_esm1b_cdr3_from_VDJ,
                                                 df_sample$LC_evo_likelihood_esm1b_cdr3_from_VDJ),
                                           method = "pearson")$estimate)
  Ablang2_ProtBERT <- c(Ablang2_ProtBERT, cor.test(x = c(df_sample$HC_evo_likelihood_ablang2_cdr3_from_VDJ,
                                                         df_sample$LC_evo_likelihood_ablang2_cdr3_from_VDJ),
                                                 y = c(df_sample$HC_evo_likelihood_protbert_cdr3_from_VDJ,
                                                       df_sample$LC_evo_likelihood_protbert_cdr3_from_VDJ),
                                                 method = "pearson")$estimate)
  Ablang2_Sapiens <- c(Ablang2_Sapiens, cor.test(x = c(df_sample$HC_evo_likelihood_ablang2_cdr3_from_VDJ,
                                                       df_sample$LC_evo_likelihood_ablang2_cdr3_from_VDJ),
                                                   y = c(df_sample$HC_evo_likelihood_sapiens_cdr3_from_VDJ,
                                                         df_sample$LC_evo_likelihood_sapiens_cdr3_from_VDJ),
                                                   method = "pearson")$estimate)
  Ablang2_Ablang1 <- c(Ablang2_Ablang1, cor.test(x = c(df_sample$HC_evo_likelihood_ablang2_cdr3_from_VDJ,
                                                       df_sample$LC_evo_likelihood_ablang2_cdr3_from_VDJ),
                                                   y = c(df_sample$HC_evo_likelihood_ablang1_cdr3_from_VDJ,
                                                         df_sample$LC_evo_likelihood_ablang1_cdr3_from_VDJ),
                                                   method = "pearson")$estimate)
  Ablang2_ESMc <- c(Ablang2_ESMc, cor.test(x = c(df_sample$HC_evo_likelihood_ablang2_cdr3_from_VDJ,
                                                 df_sample$LC_evo_likelihood_ablang2_cdr3_from_VDJ),
                                         y = c(df_sample$HC_evo_likelihood_esmc_cdr3_from_VDJ,
                                               df_sample$LC_evo_likelihood_esmc_cdr3_from_VDJ),
                                         method = "pearson")$estimate)
  ESMc_ESM1b <- c(ESMc_ESM1b, cor.test(x = c(df_sample$HC_evo_likelihood_esmc_cdr3_from_VDJ,
                                             df_sample$LC_evo_likelihood_esmc_cdr3_from_VDJ),
                                         y = c(df_sample$HC_evo_likelihood_esm1b_cdr3_from_VDJ,
                                               df_sample$LC_evo_likelihood_esm1b_cdr3_from_VDJ),
                                         method = "pearson")$estimate)
  ESMc_ProtBERT <- c(ESMc_ProtBERT, cor.test(x = c(df_sample$HC_evo_likelihood_esmc_cdr3_from_VDJ,
                                                   df_sample$LC_evo_likelihood_esmc_cdr3_from_VDJ),
                                           y = c(df_sample$HC_evo_likelihood_protbert_cdr3_from_VDJ,
                                                 df_sample$LC_evo_likelihood_protbert_cdr3_from_VDJ),
                                           method = "pearson")$estimate)
  ESMc_Sapiens <- c(ESMc_Sapiens, cor.test(x = c(df_sample$HC_evo_likelihood_esmc_cdr3_from_VDJ,
                                                   df_sample$LC_evo_likelihood_esmc_cdr3_from_VDJ),
                                               y = c(df_sample$HC_evo_likelihood_sapiens_cdr3_from_VDJ,
                                                     df_sample$LC_evo_likelihood_sapiens_cdr3_from_VDJ),
                                               method = "pearson")$estimate)
  ESMc_Ablang1 <- c(ESMc_Ablang1, cor.test(x = c(df_sample$HC_evo_likelihood_esmc_cdr3_from_VDJ,
                                                   df_sample$LC_evo_likelihood_esmc_cdr3_from_VDJ),
                                               y = c(df_sample$HC_evo_likelihood_ablang1_cdr3_from_VDJ,
                                                     df_sample$LC_evo_likelihood_ablang1_cdr3_from_VDJ),
                                               method = "pearson")$estimate)
  
  ## CDR3 only
  Ablang1_ESM1b <- c(Ablang1_ESM1b, cor.test(x = c(df_sample$HC_evo_likelihood_ablang1_cdr3_only,
                                                   df_sample$LC_evo_likelihood_ablang1_cdr3_only),
                                             y = c(df_sample$HC_evo_likelihood_esm1b_cdr3_only,
                                                   df_sample$LC_evo_likelihood_esm1b_cdr3_only),
                                             method = "pearson")$estimate)
  Ablang1_ProtBERT <- c(Ablang1_ProtBERT, cor.test(x = c(df_sample$HC_evo_likelihood_ablang1_cdr3_only,
                                                         df_sample$LC_evo_likelihood_ablang1_cdr3_only),
                                                 y = c(df_sample$HC_evo_likelihood_protbert_cdr3_only,
                                                       df_sample$LC_evo_likelihood_protbert_cdr3_only),
                                                 method = "pearson")$estimate)
  Ablang1_Sapiens <- c(Ablang1_Sapiens, cor.test(x = c(df_sample$HC_evo_likelihood_ablang1_cdr3_only,
                                                       df_sample$LC_evo_likelihood_ablang1_cdr3_only),
                                               y = c(df_sample$HC_evo_likelihood_sapiens_cdr3_only,
                                                     df_sample$LC_evo_likelihood_sapiens_cdr3_only),
                                               method = "pearson")$estimate)
  ProtBERT_ESM1b <- c(ProtBERT_ESM1b, cor.test(x = c(df_sample$HC_evo_likelihood_protbert_cdr3_only,
                                                   df_sample$LC_evo_likelihood_protbert_cdr3_only),
                                           y = c(df_sample$HC_evo_likelihood_esm1b_cdr3_only,
                                                 df_sample$LC_evo_likelihood_esm1b_cdr3_only),
                                           method = "pearson")$estimate)
  Sapiens_ESM1b <- c(Sapiens_ESM1b, cor.test(x = c(df_sample$HC_evo_likelihood_sapiens_cdr3_only,
                                                   df_sample$LC_evo_likelihood_sapiens_cdr3_only),
                                         y = c(df_sample$HC_evo_likelihood_esm1b_cdr3_only,
                                               df_sample$LC_evo_likelihood_esm1b_cdr3_only),
                                         method = "pearson")$estimate)
  Sapiens_ProtBERT <- c(Sapiens_ProtBERT, cor.test(x = c(df_sample$HC_evo_likelihood_sapiens_cdr3_only,
                                                         df_sample$LC_evo_likelihood_sapiens_cdr3_only),
                                                   y = c(df_sample$HC_evo_likelihood_protbert_cdr3_only,
                                                         df_sample$LC_evo_likelihood_protbert_cdr3_only),
                                                   method = "pearson")$estimate)
  Ablang2_ESM1b <- c(Ablang2_ESM1b, cor.test(x = c(df_sample$HC_evo_likelihood_ablang2_cdr3_only,
                                                   df_sample$LC_evo_likelihood_ablang2_cdr3_only),
                                           y = c(df_sample$HC_evo_likelihood_esm1b_cdr3_only,
                                                 df_sample$LC_evo_likelihood_esm1b_cdr3_only),
                                           method = "pearson")$estimate)
  Ablang2_ProtBERT <- c(Ablang2_ProtBERT, cor.test(x = c(df_sample$HC_evo_likelihood_ablang2_cdr3_only,
                                                         df_sample$LC_evo_likelihood_ablang2_cdr3_only),
                                                 y = c(df_sample$HC_evo_likelihood_protbert_cdr3_only,
                                                       df_sample$LC_evo_likelihood_protbert_cdr3_only),
                                                 method = "pearson")$estimate)
  Ablang2_Sapiens <- c(Ablang2_Sapiens, cor.test(x = c(df_sample$HC_evo_likelihood_ablang2_cdr3_only,
                                                       df_sample$LC_evo_likelihood_ablang2_cdr3_only),
                                                   y = c(df_sample$HC_evo_likelihood_sapiens_cdr3_only,
                                                         df_sample$LC_evo_likelihood_sapiens_cdr3_only),
                                                   method = "pearson")$estimate)
  Ablang2_Ablang1 <- c(Ablang2_Ablang1, cor.test(x = c(df_sample$HC_evo_likelihood_ablang2_cdr3_only,
                                                       df_sample$LC_evo_likelihood_ablang2_cdr3_only),
                                                   y = c(df_sample$HC_evo_likelihood_ablang1_cdr3_only,
                                                         df_sample$LC_evo_likelihood_ablang1_cdr3_only),
                                                   method = "pearson")$estimate)
  Ablang2_ESMc <- c(Ablang2_ESMc, cor.test(x = c(df_sample$HC_evo_likelihood_ablang2_cdr3_only,
                                                 df_sample$LC_evo_likelihood_ablang2_cdr3_only),
                                         y = c(df_sample$HC_evo_likelihood_esmc_cdr3_only,
                                               df_sample$LC_evo_likelihood_esmc_cdr3_only),
                                         method = "pearson")$estimate)
  ESMc_ESM1b <- c(ESMc_ESM1b, cor.test(x = c(df_sample$HC_evo_likelihood_esmc_cdr3_only,
                                             df_sample$LC_evo_likelihood_esmc_cdr3_only),
                                         y = c(df_sample$HC_evo_likelihood_esm1b_cdr3_only,
                                               df_sample$LC_evo_likelihood_esm1b_cdr3_only),
                                         method = "pearson")$estimate)
  ESMc_ProtBERT <- c(ESMc_ProtBERT, cor.test(x = c(df_sample$HC_evo_likelihood_esmc_cdr3_only,
                                                   df_sample$LC_evo_likelihood_esmc_cdr3_only),
                                           y = c(df_sample$HC_evo_likelihood_protbert_cdr3_only,
                                                 df_sample$LC_evo_likelihood_protbert_cdr3_only),
                                           method = "pearson")$estimate)
  ESMc_Sapiens <- c(ESMc_Sapiens, cor.test(x = c(df_sample$HC_evo_likelihood_esmc_cdr3_only,
                                                   df_sample$LC_evo_likelihood_esmc_cdr3_only),
                                               y = c(df_sample$HC_evo_likelihood_sapiens_cdr3_only,
                                                     df_sample$LC_evo_likelihood_sapiens_cdr3_only),
                                               method = "pearson")$estimate)
  ESMc_Ablang1 <- c(ESMc_Ablang1, cor.test(x = c(df_sample$HC_evo_likelihood_esmc_cdr3_only,
                                                   df_sample$LC_evo_likelihood_esmc_cdr3_only),
                                               y = c(df_sample$HC_evo_likelihood_ablang1_cdr3_only,
                                                     df_sample$LC_evo_likelihood_ablang1_cdr3_only),
                                               method = "pearson")$estimate)
  
  #Combine into a dataframe
  cor_df_sample <- data.frame(Ablang1_ESM1b, Ablang1_ProtBERT, Ablang1_Sapiens, ProtBERT_ESM1b, Sapiens_ESM1b, Sapiens_ProtBERT,
                             Ablang2_ESM1b, Ablang2_ProtBERT, Ablang2_Sapiens, Ablang2_Ablang1, Ablang2_ESMc, ESMc_ESM1b, 
                             ESMc_ProtBERT, ESMc_Sapiens, ESMc_Ablang1, sample, source = c("Full VDJ", "CDR3 from VDJ", "CDR3 only"))
  #Append to the main dataframe
  cor_df <- rbind(cor_df, cor_df_sample)

}

write.csv(cor_df, file = paste0("../data/",dataset,"/PLMCorrelation.csv"), row.names = F)

