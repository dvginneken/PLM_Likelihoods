library(dplyr)

args = commandArgs(trailingOnly=TRUE)
dataset <- args[1]

#Read the data
load(paste0("../data/",dataset,"/VDJ_PLL_",dataset,".RData")) #vdj_likelihood
vdj <- vdj_likelihood

#List the sample
samples <- list.dirs(path = paste0("../data/",dataset,"/VDJ"), full.names = FALSE, recursive = FALSE)

HC_cor_df <- c()
LC_cor_df <- c()
for(sample in samples){
  #Subset the dataframe per sample
  df_sample <- vdj[vdj$sample_id == sample,]
  
  ###HC
  HC_Ablang1_ESM1b = c()
  HC_Ablang1_ProtBERT = c()
  HC_Ablang1_Sapiens = c()
  HC_ProtBERT_ESM1b = c()
  HC_Sapiens_ESM1b = c()
  HC_Sapiens_ProtBERT = c()
  HC_Ablang2_ESM1b = c()
  HC_Ablang2_ProtBERT = c()
  HC_Ablang2_Sapiens = c()
  HC_Ablang2_Ablang1 = c()
  HC_Ablang2_ESMc = c()
  HC_ESMc_ESM1b = c()
  HC_ESMc_ProtBERT = c()
  HC_ESMc_Sapiens = c()
  HC_ESMc_Ablang1 = c()
  
  ## Full VDJ
  HC_Ablang1_ESM1b <- c(HC_Ablang1_ESM1b, cor.test(x = df_sample$HC_evo_likelihood_ablang1_full_VDJ,
                                             y = df_sample$HC_evo_likelihood_esm1b_full_VDJ,
                                             method = "pearson")$estimate)
  HC_Ablang1_ProtBERT <- c(HC_Ablang1_ProtBERT, cor.test(x = df_sample$HC_evo_likelihood_ablang1_full_VDJ,
                                                   y = df_sample$HC_evo_likelihood_protbert_full_VDJ,
                                                   method = "pearson")$estimate)
  HC_Ablang1_Sapiens <- c(HC_Ablang1_Sapiens, cor.test(x = df_sample$HC_evo_likelihood_ablang1_full_VDJ,
                                                 y = df_sample$HC_evo_likelihood_sapiens_full_VDJ,
                                                 method = "pearson")$estimate)
  HC_ProtBERT_ESM1b <- c(HC_ProtBERT_ESM1b, cor.test(x = df_sample$HC_evo_likelihood_protbert_full_VDJ,
                                               y = df_sample$HC_evo_likelihood_esm1b_full_VDJ,
                                               method = "pearson")$estimate)
  HC_Sapiens_ESM1b <- c(HC_Sapiens_ESM1b, cor.test(x = df_sample$HC_evo_likelihood_sapiens_full_VDJ,
                                             y = df_sample$HC_evo_likelihood_esm1b_full_VDJ,
                                             method = "pearson")$estimate)
  HC_Sapiens_ProtBERT <- c(HC_Sapiens_ProtBERT, cor.test(x = df_sample$HC_evo_likelihood_sapiens_full_VDJ,
                                                   y = df_sample$HC_evo_likelihood_protbert_full_VDJ,
                                                   method = "pearson")$estimate)
  HC_Ablang2_ESM1b <- c(HC_Ablang2_ESM1b, cor.test(x = df_sample$HC_evo_likelihood_ablang2_full_VDJ,
                                             y = df_sample$HC_evo_likelihood_esm1b_full_VDJ,
                                             method = "pearson")$estimate)
  HC_Ablang2_ProtBERT <- c(HC_Ablang2_ProtBERT, cor.test(x = df_sample$HC_evo_likelihood_ablang2_full_VDJ,
                                                   y = df_sample$HC_evo_likelihood_protbert_full_VDJ,
                                                   method = "pearson")$estimate)
  HC_Ablang2_Sapiens <- c(HC_Ablang2_Sapiens, cor.test(x = df_sample$HC_evo_likelihood_ablang2_full_VDJ,
                                                 y = df_sample$HC_evo_likelihood_sapiens_full_VDJ,
                                                 method = "pearson")$estimate)
  HC_Ablang2_Ablang1 <- c(HC_Ablang2_Ablang1, cor.test(x = df_sample$HC_evo_likelihood_ablang2_full_VDJ,
                                                 y = df_sample$HC_evo_likelihood_ablang1_full_VDJ,
                                                 method = "pearson")$estimate)
  HC_Ablang2_ESMc <- c(HC_Ablang2_ESMc, cor.test(x = df_sample$HC_evo_likelihood_ablang2_full_VDJ,
                                           y = df_sample$HC_evo_likelihood_esmc_full_VDJ,
                                           method = "pearson")$estimate)
  HC_ESMc_ESM1b <- c(HC_ESMc_ESM1b, cor.test(x = df_sample$HC_evo_likelihood_esmc_full_VDJ,
                                       y = df_sample$HC_evo_likelihood_esm1b_full_VDJ,
                                       method = "pearson")$estimate)
  HC_ESMc_ProtBERT <- c(HC_ESMc_ProtBERT, cor.test(x = df_sample$HC_evo_likelihood_esmc_full_VDJ,
                                             y = df_sample$HC_evo_likelihood_protbert_full_VDJ,
                                             method = "pearson")$estimate)
  HC_ESMc_Sapiens <- c(HC_ESMc_Sapiens, cor.test(x = df_sample$HC_evo_likelihood_esmc_full_VDJ,
                                           y = df_sample$HC_evo_likelihood_sapiens_full_VDJ,
                                           method = "pearson")$estimate)
  HC_ESMc_Ablang1 <- c(HC_ESMc_Ablang1, cor.test(x = df_sample$HC_evo_likelihood_esmc_full_VDJ,
                                           y = df_sample$HC_evo_likelihood_ablang1_full_VDJ,
                                           method = "pearson")$estimate)
  
  ## CDR3 from VDJ
  HC_Ablang1_ESM1b <- c(HC_Ablang1_ESM1b, cor.test(x = df_sample$HC_evo_likelihood_ablang1_cdr3_from_VDJ,
                                                   y = df_sample$HC_evo_likelihood_esm1b_cdr3_from_VDJ,
                                                   method = "pearson")$estimate)
  HC_Ablang1_ProtBERT <- c(HC_Ablang1_ProtBERT, cor.test(x = df_sample$HC_evo_likelihood_ablang1_cdr3_from_VDJ,
                                                         y = df_sample$HC_evo_likelihood_protbert_cdr3_from_VDJ,
                                                         method = "pearson")$estimate)
  HC_Ablang1_Sapiens <- c(HC_Ablang1_Sapiens, cor.test(x = df_sample$HC_evo_likelihood_ablang1_cdr3_from_VDJ,
                                                       y = df_sample$HC_evo_likelihood_sapiens_cdr3_from_VDJ,
                                                       method = "pearson")$estimate)
  HC_ProtBERT_ESM1b <- c(HC_ProtBERT_ESM1b, cor.test(x = df_sample$HC_evo_likelihood_protbert_cdr3_from_VDJ,
                                                     y = df_sample$HC_evo_likelihood_esm1b_cdr3_from_VDJ,
                                                     method = "pearson")$estimate)
  HC_Sapiens_ESM1b <- c(HC_Sapiens_ESM1b, cor.test(x = df_sample$HC_evo_likelihood_sapiens_cdr3_from_VDJ,
                                                   y = df_sample$HC_evo_likelihood_esm1b_cdr3_from_VDJ,
                                                   method = "pearson")$estimate)
  HC_Sapiens_ProtBERT <- c(HC_Sapiens_ProtBERT, cor.test(x = df_sample$HC_evo_likelihood_sapiens_cdr3_from_VDJ,
                                                         y = df_sample$HC_evo_likelihood_protbert_cdr3_from_VDJ,
                                                         method = "pearson")$estimate)
  HC_Ablang2_ESM1b <- c(HC_Ablang2_ESM1b, cor.test(x = df_sample$HC_evo_likelihood_ablang2_cdr3_from_VDJ,
                                                   y = df_sample$HC_evo_likelihood_esm1b_cdr3_from_VDJ,
                                                   method = "pearson")$estimate)
  HC_Ablang2_ProtBERT <- c(HC_Ablang2_ProtBERT, cor.test(x = df_sample$HC_evo_likelihood_ablang2_cdr3_from_VDJ,
                                                         y = df_sample$HC_evo_likelihood_protbert_cdr3_from_VDJ,
                                                         method = "pearson")$estimate)
  HC_Ablang2_Sapiens <- c(HC_Ablang2_Sapiens, cor.test(x = df_sample$HC_evo_likelihood_ablang2_cdr3_from_VDJ,
                                                       y = df_sample$HC_evo_likelihood_sapiens_cdr3_from_VDJ,
                                                       method = "pearson")$estimate)
  HC_Ablang2_Ablang1 <- c(HC_Ablang2_Ablang1, cor.test(x = df_sample$HC_evo_likelihood_ablang2_cdr3_from_VDJ,
                                                       y = df_sample$HC_evo_likelihood_ablang1_cdr3_from_VDJ,
                                                       method = "pearson")$estimate)
  HC_Ablang2_ESMc <- c(HC_Ablang2_ESMc, cor.test(x = df_sample$HC_evo_likelihood_ablang2_cdr3_from_VDJ,
                                                 y = df_sample$HC_evo_likelihood_esmc_cdr3_from_VDJ,
                                                 method = "pearson")$estimate)
  HC_ESMc_ESM1b <- c(HC_ESMc_ESM1b, cor.test(x = df_sample$HC_evo_likelihood_esmc_cdr3_from_VDJ,
                                             y = df_sample$HC_evo_likelihood_esm1b_cdr3_from_VDJ,
                                             method = "pearson")$estimate)
  HC_ESMc_ProtBERT <- c(HC_ESMc_ProtBERT, cor.test(x = df_sample$HC_evo_likelihood_esmc_cdr3_from_VDJ,
                                                   y = df_sample$HC_evo_likelihood_protbert_cdr3_from_VDJ,
                                                   method = "pearson")$estimate)
  HC_ESMc_Sapiens <- c(HC_ESMc_Sapiens, cor.test(x = df_sample$HC_evo_likelihood_esmc_cdr3_from_VDJ,
                                                 y = df_sample$HC_evo_likelihood_sapiens_cdr3_from_VDJ,
                                                 method = "pearson")$estimate)
  HC_ESMc_Ablang1 <- c(HC_ESMc_Ablang1, cor.test(x = df_sample$HC_evo_likelihood_esmc_cdr3_from_VDJ,
                                                 y = df_sample$HC_evo_likelihood_ablang1_cdr3_from_VDJ,
                                                 method = "pearson")$estimate)
  
  ## CDR3 only
  HC_Ablang1_ESM1b <- c(HC_Ablang1_ESM1b, cor.test(x = df_sample$HC_evo_likelihood_ablang1_cdr3_only,
                                                   y = df_sample$HC_evo_likelihood_esm1b_cdr3_only,
                                                   method = "pearson")$estimate)
  HC_Ablang1_ProtBERT <- c(HC_Ablang1_ProtBERT, cor.test(x = df_sample$HC_evo_likelihood_ablang1_cdr3_only,
                                                         y = df_sample$HC_evo_likelihood_protbert_cdr3_only,
                                                         method = "pearson")$estimate)
  HC_Ablang1_Sapiens <- c(HC_Ablang1_Sapiens, cor.test(x = df_sample$HC_evo_likelihood_ablang1_cdr3_only,
                                                       y = df_sample$HC_evo_likelihood_sapiens_cdr3_only,
                                                       method = "pearson")$estimate)
  HC_ProtBERT_ESM1b <- c(HC_ProtBERT_ESM1b, cor.test(x = df_sample$HC_evo_likelihood_protbert_cdr3_only,
                                                     y = df_sample$HC_evo_likelihood_esm1b_cdr3_only,
                                                     method = "pearson")$estimate)
  HC_Sapiens_ESM1b <- c(HC_Sapiens_ESM1b, cor.test(x = df_sample$HC_evo_likelihood_sapiens_cdr3_only,
                                                   y = df_sample$HC_evo_likelihood_esm1b_cdr3_only,
                                                   method = "pearson")$estimate)
  HC_Sapiens_ProtBERT <- c(HC_Sapiens_ProtBERT, cor.test(x = df_sample$HC_evo_likelihood_sapiens_cdr3_only,
                                                         y = df_sample$HC_evo_likelihood_protbert_cdr3_only,
                                                         method = "pearson")$estimate)
  HC_Ablang2_ESM1b <- c(HC_Ablang2_ESM1b, cor.test(x = df_sample$HC_evo_likelihood_ablang2_cdr3_only,
                                                   y = df_sample$HC_evo_likelihood_esm1b_cdr3_only,
                                                   method = "pearson")$estimate)
  HC_Ablang2_ProtBERT <- c(HC_Ablang2_ProtBERT, cor.test(x = df_sample$HC_evo_likelihood_ablang2_cdr3_only,
                                                         y = df_sample$HC_evo_likelihood_protbert_cdr3_only,
                                                         method = "pearson")$estimate)
  HC_Ablang2_Sapiens <- c(HC_Ablang2_Sapiens, cor.test(x = df_sample$HC_evo_likelihood_ablang2_cdr3_only,
                                                       y = df_sample$HC_evo_likelihood_sapiens_cdr3_only,
                                                       method = "pearson")$estimate)
  HC_Ablang2_Ablang1 <- c(HC_Ablang2_Ablang1, cor.test(x = df_sample$HC_evo_likelihood_ablang2_cdr3_only,
                                                       y = df_sample$HC_evo_likelihood_ablang1_cdr3_only,
                                                       method = "pearson")$estimate)
  HC_Ablang2_ESMc <- c(HC_Ablang2_ESMc, cor.test(x = df_sample$HC_evo_likelihood_ablang2_cdr3_only,
                                                 y = df_sample$HC_evo_likelihood_esmc_cdr3_only,
                                                 method = "pearson")$estimate)
  HC_ESMc_ESM1b <- c(HC_ESMc_ESM1b, cor.test(x = df_sample$HC_evo_likelihood_esmc_cdr3_only,
                                             y = df_sample$HC_evo_likelihood_esm1b_cdr3_only,
                                             method = "pearson")$estimate)
  HC_ESMc_ProtBERT <- c(HC_ESMc_ProtBERT, cor.test(x = df_sample$HC_evo_likelihood_esmc_cdr3_only,
                                                   y = df_sample$HC_evo_likelihood_protbert_cdr3_only,
                                                   method = "pearson")$estimate)
  HC_ESMc_Sapiens <- c(HC_ESMc_Sapiens, cor.test(x = df_sample$HC_evo_likelihood_esmc_cdr3_only,
                                                 y = df_sample$HC_evo_likelihood_sapiens_cdr3_only,
                                                 method = "pearson")$estimate)
  HC_ESMc_Ablang1 <- c(HC_ESMc_Ablang1, cor.test(x = df_sample$HC_evo_likelihood_esmc_cdr3_only,
                                                 y = df_sample$HC_evo_likelihood_ablang1_cdr3_only,
                                                 method = "pearson")$estimate)
  
  #Combine into a dataframe
  HC_cor_df_sample <- data.frame(HC_Ablang1_ESM1b, HC_Ablang1_ProtBERT, HC_Ablang1_Sapiens,
                                 HC_ProtBERT_ESM1b, HC_Sapiens_ESM1b, HC_Sapiens_ProtBERT,
                                 HC_Ablang2_ESM1b, HC_Ablang2_ProtBERT, HC_Ablang2_Sapiens,
                                 HC_Ablang2_Ablang1, HC_Ablang2_ESMc, HC_ESMc_ESM1b,
                                 HC_ESMc_ProtBERT, HC_ESMc_Sapiens, HC_ESMc_Ablang1,
                                 sample = sample, chain = "HC",
                                 source = c("Full VDJ", "CDR3 from VDJ", "CDR3 only"))
  
  ###LC
  LC_Ablang1_ESM1b = c()
  LC_Ablang1_ProtBERT = c()
  LC_Ablang1_Sapiens = c()
  LC_ProtBERT_ESM1b = c()
  LC_Sapiens_ESM1b = c()
  LC_Sapiens_ProtBERT = c()
  LC_Ablang2_ESM1b = c()
  LC_Ablang2_ProtBERT = c()
  LC_Ablang2_Sapiens = c()
  LC_Ablang2_Ablang1 = c()
  LC_Ablang2_ESMc = c()
  LC_ESMc_ESM1b = c()
  LC_ESMc_ProtBERT = c()
  LC_ESMc_Sapiens = c()
  LC_ESMc_Ablang1 = c()
  
  ## Full VDJ
  LC_Ablang1_ESM1b <- c(LC_Ablang1_ESM1b, cor.test(x = df_sample$LC_evo_likelihood_ablang1_full_VDJ,
                                                   y = df_sample$LC_evo_likelihood_esm1b_full_VDJ,
                                                   method = "pearson")$estimate)
  LC_Ablang1_ProtBERT <- c(LC_Ablang1_ProtBERT, cor.test(x = df_sample$LC_evo_likelihood_ablang1_full_VDJ,
                                                         y = df_sample$LC_evo_likelihood_protbert_full_VDJ,
                                                         method = "pearson")$estimate)
  LC_Ablang1_Sapiens <- c(LC_Ablang1_Sapiens, cor.test(x = df_sample$LC_evo_likelihood_ablang1_full_VDJ,
                                                       y = df_sample$LC_evo_likelihood_sapiens_full_VDJ,
                                                       method = "pearson")$estimate)
  LC_ProtBERT_ESM1b <- c(LC_ProtBERT_ESM1b, cor.test(x = df_sample$LC_evo_likelihood_protbert_full_VDJ,
                                                     y = df_sample$LC_evo_likelihood_esm1b_full_VDJ,
                                                     method = "pearson")$estimate)
  LC_Sapiens_ESM1b <- c(LC_Sapiens_ESM1b, cor.test(x = df_sample$LC_evo_likelihood_sapiens_full_VDJ,
                                                   y = df_sample$LC_evo_likelihood_esm1b_full_VDJ,
                                                   method = "pearson")$estimate)
  LC_Sapiens_ProtBERT <- c(LC_Sapiens_ProtBERT, cor.test(x = df_sample$LC_evo_likelihood_sapiens_full_VDJ,
                                                         y = df_sample$LC_evo_likelihood_protbert_full_VDJ,
                                                         method = "pearson")$estimate)
  LC_Ablang2_ESM1b <- c(LC_Ablang2_ESM1b, cor.test(x = df_sample$LC_evo_likelihood_ablang2_full_VDJ,
                                                   y = df_sample$LC_evo_likelihood_esm1b_full_VDJ,
                                                   method = "pearson")$estimate)
  LC_Ablang2_ProtBERT <- c(LC_Ablang2_ProtBERT, cor.test(x = df_sample$LC_evo_likelihood_ablang2_full_VDJ,
                                                         y = df_sample$LC_evo_likelihood_protbert_full_VDJ,
                                                         method = "pearson")$estimate)
  LC_Ablang2_Sapiens <- c(LC_Ablang2_Sapiens, cor.test(x = df_sample$LC_evo_likelihood_ablang2_full_VDJ,
                                                       y = df_sample$LC_evo_likelihood_sapiens_full_VDJ,
                                                       method = "pearson")$estimate)
  LC_Ablang2_Ablang1 <- c(LC_Ablang2_Ablang1, cor.test(x = df_sample$LC_evo_likelihood_ablang2_full_VDJ,
                                                       y = df_sample$LC_evo_likelihood_ablang1_full_VDJ,
                                                       method = "pearson")$estimate)
  LC_Ablang2_ESMc <- c(LC_Ablang2_ESMc, cor.test(x = df_sample$LC_evo_likelihood_ablang2_full_VDJ,
                                                 y = df_sample$LC_evo_likelihood_esmc_full_VDJ,
                                                 method = "pearson")$estimate)
  LC_ESMc_ESM1b <- c(LC_ESMc_ESM1b, cor.test(x = df_sample$LC_evo_likelihood_esmc_full_VDJ,
                                             y = df_sample$LC_evo_likelihood_esm1b_full_VDJ,
                                             method = "pearson")$estimate)
  LC_ESMc_ProtBERT <- c(LC_ESMc_ProtBERT, cor.test(x = df_sample$LC_evo_likelihood_esmc_full_VDJ,
                                                   y = df_sample$LC_evo_likelihood_protbert_full_VDJ,
                                                   method = "pearson")$estimate)
  LC_ESMc_Sapiens <- c(LC_ESMc_Sapiens, cor.test(x = df_sample$LC_evo_likelihood_esmc_full_VDJ,
                                                 y = df_sample$LC_evo_likelihood_sapiens_full_VDJ,
                                                 method = "pearson")$estimate)
  LC_ESMc_Ablang1 <- c(LC_ESMc_Ablang1, cor.test(x = df_sample$LC_evo_likelihood_esmc_full_VDJ,
                                                 y = df_sample$LC_evo_likelihood_ablang1_full_VDJ,
                                                 method = "pearson")$estimate)
  
  ## CDR3 from VDJ
  LC_Ablang1_ESM1b <- c(LC_Ablang1_ESM1b, cor.test(x = df_sample$LC_evo_likelihood_ablang1_cdr3_from_VDJ,
                                                   y = df_sample$LC_evo_likelihood_esm1b_cdr3_from_VDJ,
                                                   method = "pearson")$estimate)
  LC_Ablang1_ProtBERT <- c(LC_Ablang1_ProtBERT, cor.test(x = df_sample$LC_evo_likelihood_ablang1_cdr3_from_VDJ,
                                                         y = df_sample$LC_evo_likelihood_protbert_cdr3_from_VDJ,
                                                         method = "pearson")$estimate)
  LC_Ablang1_Sapiens <- c(LC_Ablang1_Sapiens, cor.test(x = df_sample$LC_evo_likelihood_ablang1_cdr3_from_VDJ,
                                                       y = df_sample$LC_evo_likelihood_sapiens_cdr3_from_VDJ,
                                                       method = "pearson")$estimate)
  LC_ProtBERT_ESM1b <- c(LC_ProtBERT_ESM1b, cor.test(x = df_sample$LC_evo_likelihood_protbert_cdr3_from_VDJ,
                                                     y = df_sample$LC_evo_likelihood_esm1b_cdr3_from_VDJ,
                                                     method = "pearson")$estimate)
  LC_Sapiens_ESM1b <- c(LC_Sapiens_ESM1b, cor.test(x = df_sample$LC_evo_likelihood_sapiens_cdr3_from_VDJ,
                                                   y = df_sample$LC_evo_likelihood_esm1b_cdr3_from_VDJ,
                                                   method = "pearson")$estimate)
  LC_Sapiens_ProtBERT <- c(LC_Sapiens_ProtBERT, cor.test(x = df_sample$LC_evo_likelihood_sapiens_cdr3_from_VDJ,
                                                         y = df_sample$LC_evo_likelihood_protbert_cdr3_from_VDJ,
                                                         method = "pearson")$estimate)
  LC_Ablang2_ESM1b <- c(LC_Ablang2_ESM1b, cor.test(x = df_sample$LC_evo_likelihood_ablang2_cdr3_from_VDJ,
                                                   y = df_sample$LC_evo_likelihood_esm1b_cdr3_from_VDJ,
                                                   method = "pearson")$estimate)
  LC_Ablang2_ProtBERT <- c(LC_Ablang2_ProtBERT, cor.test(x = df_sample$LC_evo_likelihood_ablang2_cdr3_from_VDJ,
                                                         y = df_sample$LC_evo_likelihood_protbert_cdr3_from_VDJ,
                                                         method = "pearson")$estimate)
  LC_Ablang2_Sapiens <- c(LC_Ablang2_Sapiens, cor.test(x = df_sample$LC_evo_likelihood_ablang2_cdr3_from_VDJ,
                                                       y = df_sample$LC_evo_likelihood_sapiens_cdr3_from_VDJ,
                                                       method = "pearson")$estimate)
  LC_Ablang2_Ablang1 <- c(LC_Ablang2_Ablang1, cor.test(x = df_sample$LC_evo_likelihood_ablang2_cdr3_from_VDJ,
                                                       y = df_sample$LC_evo_likelihood_ablang1_cdr3_from_VDJ,
                                                       method = "pearson")$estimate)
  LC_Ablang2_ESMc <- c(LC_Ablang2_ESMc, cor.test(x = df_sample$LC_evo_likelihood_ablang2_cdr3_from_VDJ,
                                                 y = df_sample$LC_evo_likelihood_esmc_cdr3_from_VDJ,
                                                 method = "pearson")$estimate)
  LC_ESMc_ESM1b <- c(LC_ESMc_ESM1b, cor.test(x = df_sample$LC_evo_likelihood_esmc_cdr3_from_VDJ,
                                             y = df_sample$LC_evo_likelihood_esm1b_cdr3_from_VDJ,
                                             method = "pearson")$estimate)
  LC_ESMc_ProtBERT <- c(LC_ESMc_ProtBERT, cor.test(x = df_sample$LC_evo_likelihood_esmc_cdr3_from_VDJ,
                                                   y = df_sample$LC_evo_likelihood_protbert_cdr3_from_VDJ,
                                                   method = "pearson")$estimate)
  LC_ESMc_Sapiens <- c(LC_ESMc_Sapiens, cor.test(x = df_sample$LC_evo_likelihood_esmc_cdr3_from_VDJ,
                                                 y = df_sample$LC_evo_likelihood_sapiens_cdr3_from_VDJ,
                                                 method = "pearson")$estimate)
  LC_ESMc_Ablang1 <- c(LC_ESMc_Ablang1, cor.test(x = df_sample$LC_evo_likelihood_esmc_cdr3_from_VDJ,
                                                 y = df_sample$LC_evo_likelihood_ablang1_cdr3_from_VDJ,
                                                 method = "pearson")$estimate)
  
  ## CDR3 only
  LC_Ablang1_ESM1b <- c(LC_Ablang1_ESM1b, cor.test(x = df_sample$LC_evo_likelihood_ablang1_cdr3_only,
                                                   y = df_sample$LC_evo_likelihood_esm1b_cdr3_only,
                                                   method = "pearson")$estimate)
  LC_Ablang1_ProtBERT <- c(LC_Ablang1_ProtBERT, cor.test(x = df_sample$LC_evo_likelihood_ablang1_cdr3_only,
                                                         y = df_sample$LC_evo_likelihood_protbert_cdr3_only,
                                                         method = "pearson")$estimate)
  LC_Ablang1_Sapiens <- c(LC_Ablang1_Sapiens, cor.test(x = df_sample$LC_evo_likelihood_ablang1_cdr3_only,
                                                       y = df_sample$LC_evo_likelihood_sapiens_cdr3_only,
                                                       method = "pearson")$estimate)
  LC_ProtBERT_ESM1b <- c(LC_ProtBERT_ESM1b, cor.test(x = df_sample$LC_evo_likelihood_protbert_cdr3_only,
                                                     y = df_sample$LC_evo_likelihood_esm1b_cdr3_only,
                                                     method = "pearson")$estimate)
  LC_Sapiens_ESM1b <- c(LC_Sapiens_ESM1b, cor.test(x = df_sample$LC_evo_likelihood_sapiens_cdr3_only,
                                                   y = df_sample$LC_evo_likelihood_esm1b_cdr3_only,
                                                   method = "pearson")$estimate)
  LC_Sapiens_ProtBERT <- c(LC_Sapiens_ProtBERT, cor.test(x = df_sample$LC_evo_likelihood_sapiens_cdr3_only,
                                                         y = df_sample$LC_evo_likelihood_protbert_cdr3_only,
                                                         method = "pearson")$estimate)
  LC_Ablang2_ESM1b <- c(LC_Ablang2_ESM1b, cor.test(x = df_sample$LC_evo_likelihood_ablang2_cdr3_only,
                                                   y = df_sample$LC_evo_likelihood_esm1b_cdr3_only,
                                                   method = "pearson")$estimate)
  LC_Ablang2_ProtBERT <- c(LC_Ablang2_ProtBERT, cor.test(x = df_sample$LC_evo_likelihood_ablang2_cdr3_only,
                                                         y = df_sample$LC_evo_likelihood_protbert_cdr3_only,
                                                         method = "pearson")$estimate)
  LC_Ablang2_Sapiens <- c(LC_Ablang2_Sapiens, cor.test(x = df_sample$LC_evo_likelihood_ablang2_cdr3_only,
                                                       y = df_sample$LC_evo_likelihood_sapiens_cdr3_only,
                                                       method = "pearson")$estimate)
  LC_Ablang2_Ablang1 <- c(LC_Ablang2_Ablang1, cor.test(x = df_sample$LC_evo_likelihood_ablang2_cdr3_only,
                                                       y = df_sample$LC_evo_likelihood_ablang1_cdr3_only,
                                                       method = "pearson")$estimate)
  LC_Ablang2_ESMc <- c(LC_Ablang2_ESMc, cor.test(x = df_sample$LC_evo_likelihood_ablang2_cdr3_only,
                                                 y = df_sample$LC_evo_likelihood_esmc_cdr3_only,
                                                 method = "pearson")$estimate)
  LC_ESMc_ESM1b <- c(LC_ESMc_ESM1b, cor.test(x = df_sample$LC_evo_likelihood_esmc_cdr3_only,
                                             y = df_sample$LC_evo_likelihood_esm1b_cdr3_only,
                                             method = "pearson")$estimate)
  LC_ESMc_ProtBERT <- c(LC_ESMc_ProtBERT, cor.test(x = df_sample$LC_evo_likelihood_esmc_cdr3_only,
                                                   y = df_sample$LC_evo_likelihood_protbert_cdr3_only,
                                                   method = "pearson")$estimate)
  LC_ESMc_Sapiens <- c(LC_ESMc_Sapiens, cor.test(x = df_sample$LC_evo_likelihood_esmc_cdr3_only,
                                                 y = df_sample$LC_evo_likelihood_sapiens_cdr3_only,
                                                 method = "pearson")$estimate)
  LC_ESMc_Ablang1 <- c(LC_ESMc_Ablang1, cor.test(x = df_sample$LC_evo_likelihood_esmc_cdr3_only,
                                                 y = df_sample$LC_evo_likelihood_ablang1_cdr3_only,
                                                 method = "pearson")$estimate)
  
  #Combine into a dataframe
  LC_cor_df_sample <- data.frame(LC_Ablang1_ESM1b, LC_Ablang1_ProtBERT, LC_Ablang1_Sapiens,
                                 LC_ProtBERT_ESM1b, LC_Sapiens_ESM1b, LC_Sapiens_ProtBERT,
                                 LC_Ablang2_ESM1b, LC_Ablang2_ProtBERT, LC_Ablang2_Sapiens,
                                 LC_Ablang2_Ablang1, LC_Ablang2_ESMc, LC_ESMc_ESM1b,
                                 LC_ESMc_ProtBERT, LC_ESMc_Sapiens, LC_ESMc_Ablang1,
                                 sample = sample, chain = "LC",
                                 source = c("Full VDJ", "CDR3 from VDJ", "CDR3 only"))  
  
  
  #Append to the main dataframe
  HC_cor_df <- rbind(HC_cor_df, HC_cor_df_sample)
  LC_cor_df <- rbind(LC_cor_df, LC_cor_df_sample)
}

write.csv(HC_cor_df, file = paste0("../data/",dataset,"/PLMCorrelation_HC.csv"), row.names = F)
write.csv(LC_cor_df, file = paste0("../data/",dataset,"/PLMCorrelation_LC.csv"), row.names = F)

