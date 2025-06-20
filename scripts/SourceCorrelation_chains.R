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
  HC_full_VDJ__CDR3_only = c()
  HC_full_VDJ__CDR3_from_VDJ = c()
  HC_CDR3_only__CDR3_from_VDJ = c()
  
  
  ## Ablang1
  HC_CDR3_only__CDR3_from_VDJ <- c(HC_CDR3_only__CDR3_from_VDJ, cor.test(x = df_sample$HC_evo_likelihood_ablang1_cdr3_from_VDJ,
                                                                   y = df_sample$HC_evo_likelihood_ablang1_cdr3_only,
                                                                   method = "pearson")$estimate)
  
  HC_full_VDJ__CDR3_from_VDJ <- c(HC_full_VDJ__CDR3_from_VDJ, cor.test(x = df_sample$HC_evo_likelihood_ablang1_cdr3_from_VDJ,
                                                                 y = df_sample$HC_evo_likelihood_ablang1_full_VDJ,
                                                                 method = "pearson")$estimate)
  
  HC_full_VDJ__CDR3_only <- c(HC_full_VDJ__CDR3_only, cor.test(x = df_sample$HC_evo_likelihood_ablang1_full_VDJ,
                                                         y = df_sample$HC_evo_likelihood_ablang1_cdr3_only,
                                                         method = "pearson")$estimate)
  
  ## Ablang2
  HC_CDR3_only__CDR3_from_VDJ <- c(HC_CDR3_only__CDR3_from_VDJ, cor.test(x = df_sample$HC_evo_likelihood_ablang2_cdr3_from_VDJ,
                                                                   y = df_sample$HC_evo_likelihood_ablang2_cdr3_only,
                                                                   method = "pearson")$estimate)
  
  HC_full_VDJ__CDR3_from_VDJ <- c(HC_full_VDJ__CDR3_from_VDJ, cor.test(x = df_sample$HC_evo_likelihood_ablang2_cdr3_from_VDJ,
                                                                 y = df_sample$HC_evo_likelihood_ablang2_full_VDJ,
                                                                 method = "pearson")$estimate)
  
  HC_full_VDJ__CDR3_only <- c(HC_full_VDJ__CDR3_only, cor.test(x = df_sample$HC_evo_likelihood_ablang2_full_VDJ,
                                                         y = df_sample$HC_evo_likelihood_ablang2_cdr3_only,
                                                         method = "pearson")$estimate)
  
  ## Sapiens
  HC_CDR3_only__CDR3_from_VDJ <- c(HC_CDR3_only__CDR3_from_VDJ, cor.test(x = df_sample$HC_evo_likelihood_sapiens_cdr3_from_VDJ,
                                                                   y = df_sample$HC_evo_likelihood_sapiens_cdr3_only,
                                                                   method = "pearson")$estimate)
  
  HC_full_VDJ__CDR3_from_VDJ <- c(HC_full_VDJ__CDR3_from_VDJ, cor.test(x = df_sample$HC_evo_likelihood_sapiens_cdr3_from_VDJ,
                                                                 y = df_sample$HC_evo_likelihood_sapiens_full_VDJ,
                                                                 method = "pearson")$estimate)
  
  HC_full_VDJ__CDR3_only <- c(HC_full_VDJ__CDR3_only, cor.test(x = df_sample$HC_evo_likelihood_sapiens_full_VDJ,
                                                         y = df_sample$HC_evo_likelihood_sapiens_cdr3_only,
                                                         method = "pearson")$estimate)
  
  ## ESM1b
  HC_CDR3_only__CDR3_from_VDJ <- c(HC_CDR3_only__CDR3_from_VDJ, cor.test(x = df_sample$HC_evo_likelihood_esm1b_cdr3_from_VDJ,
                                                                   y = df_sample$HC_evo_likelihood_esm1b_cdr3_only,
                                                                   method = "pearson")$estimate)
  
  HC_full_VDJ__CDR3_from_VDJ <- c(HC_full_VDJ__CDR3_from_VDJ, cor.test(x = df_sample$HC_evo_likelihood_esm1b_cdr3_from_VDJ,
                                                                 y = df_sample$HC_evo_likelihood_esm1b_full_VDJ,
                                                                 method = "pearson")$estimate)
  
  HC_full_VDJ__CDR3_only <- c(HC_full_VDJ__CDR3_only, cor.test(x = df_sample$HC_evo_likelihood_esm1b_full_VDJ,
                                                         y = df_sample$HC_evo_likelihood_esm1b_cdr3_only,
                                                         method = "pearson")$estimate)
  
  ## ESMc
  HC_CDR3_only__CDR3_from_VDJ <- c(HC_CDR3_only__CDR3_from_VDJ, cor.test(x = df_sample$HC_evo_likelihood_esmc_cdr3_from_VDJ,
                                                                   y = df_sample$HC_evo_likelihood_esmc_cdr3_only,
                                                                   method = "pearson")$estimate)
  
  HC_full_VDJ__CDR3_from_VDJ <- c(HC_full_VDJ__CDR3_from_VDJ, cor.test(x = df_sample$HC_evo_likelihood_esmc_cdr3_from_VDJ,
                                                                 y = df_sample$HC_evo_likelihood_esmc_full_VDJ,
                                                                 method = "pearson")$estimate)
  
  HC_full_VDJ__CDR3_only <- c(HC_full_VDJ__CDR3_only, cor.test(x = df_sample$HC_evo_likelihood_esmc_full_VDJ,
                                                         y = df_sample$HC_evo_likelihood_esmc_cdr3_only,
                                                         method = "pearson")$estimate)
  
  ## ProtBERT
  HC_CDR3_only__CDR3_from_VDJ <- c(HC_CDR3_only__CDR3_from_VDJ, cor.test(x = df_sample$HC_evo_likelihood_protbert_cdr3_from_VDJ,
                                                                   y = df_sample$HC_evo_likelihood_protbert_cdr3_only,
                                                                   method = "pearson")$estimate)
  
  HC_full_VDJ__CDR3_from_VDJ <- c(HC_full_VDJ__CDR3_from_VDJ, cor.test(x = df_sample$HC_evo_likelihood_protbert_cdr3_from_VDJ,
                                                                 y = df_sample$HC_evo_likelihood_protbert_full_VDJ,
                                                                 method = "pearson")$estimate)
  
  HC_full_VDJ__CDR3_only <- c(HC_full_VDJ__CDR3_only, cor.test(x = df_sample$HC_evo_likelihood_protbert_full_VDJ,
                                                         y = df_sample$HC_evo_likelihood_protbert_cdr3_only,
                                                         method = "pearson")$estimate)
  
  #Combine into a dataframe
  HC_cor_df_sample <- data.frame(HC_full_VDJ__CDR3_only, HC_full_VDJ__CDR3_from_VDJ, HC_CDR3_only__CDR3_from_VDJ, 
                              sample, chain = "HC", model = c("Ablang1", "Ablang2", "Sapiens", "ESM-1b", "ESM-C", "ProtBERT"))
  
  
  ###LC
  LC_full_VDJ__CDR3_only = c()
  LC_full_VDJ__CDR3_from_VDJ = c()
  LC_CDR3_only__CDR3_from_VDJ = c()
  
  ## Ablang1
  LC_CDR3_only__CDR3_from_VDJ <- c(LC_CDR3_only__CDR3_from_VDJ, cor.test(x = df_sample$LC_evo_likelihood_ablang1_cdr3_from_VDJ,
                                                                         y = df_sample$LC_evo_likelihood_ablang1_cdr3_only,
                                                                         method = "pearson")$estimate)
  
  LC_full_VDJ__CDR3_from_VDJ <- c(LC_full_VDJ__CDR3_from_VDJ, cor.test(x = df_sample$LC_evo_likelihood_ablang1_cdr3_from_VDJ,
                                                                       y = df_sample$LC_evo_likelihood_ablang1_full_VDJ,
                                                                       method = "pearson")$estimate)
  
  LC_full_VDJ__CDR3_only <- c(LC_full_VDJ__CDR3_only, cor.test(x = df_sample$LC_evo_likelihood_ablang1_full_VDJ,
                                                               y = df_sample$LC_evo_likelihood_ablang1_cdr3_only,
                                                               method = "pearson")$estimate)
  
  ## Ablang2
  LC_CDR3_only__CDR3_from_VDJ <- c(LC_CDR3_only__CDR3_from_VDJ, cor.test(x = df_sample$LC_evo_likelihood_ablang2_cdr3_from_VDJ,
                                                                         y = df_sample$LC_evo_likelihood_ablang2_cdr3_only,
                                                                         method = "pearson")$estimate)
  
  LC_full_VDJ__CDR3_from_VDJ <- c(LC_full_VDJ__CDR3_from_VDJ, cor.test(x = df_sample$LC_evo_likelihood_ablang2_cdr3_from_VDJ,
                                                                       y = df_sample$LC_evo_likelihood_ablang2_full_VDJ,
                                                                       method = "pearson")$estimate)
  
  LC_full_VDJ__CDR3_only <- c(LC_full_VDJ__CDR3_only, cor.test(x = df_sample$LC_evo_likelihood_ablang2_full_VDJ,
                                                               y = df_sample$LC_evo_likelihood_ablang2_cdr3_only,
                                                               method = "pearson")$estimate)
  
  ## Sapiens
  LC_CDR3_only__CDR3_from_VDJ <- c(LC_CDR3_only__CDR3_from_VDJ, cor.test(x = df_sample$LC_evo_likelihood_sapiens_cdr3_from_VDJ,
                                                                         y = df_sample$LC_evo_likelihood_sapiens_cdr3_only,
                                                                         method = "pearson")$estimate)
  
  LC_full_VDJ__CDR3_from_VDJ <- c(LC_full_VDJ__CDR3_from_VDJ, cor.test(x = df_sample$LC_evo_likelihood_sapiens_cdr3_from_VDJ,
                                                                       y = df_sample$LC_evo_likelihood_sapiens_full_VDJ,
                                                                       method = "pearson")$estimate)
  
  LC_full_VDJ__CDR3_only <- c(LC_full_VDJ__CDR3_only, cor.test(x = df_sample$LC_evo_likelihood_sapiens_full_VDJ,
                                                               y = df_sample$LC_evo_likelihood_sapiens_cdr3_only,
                                                               method = "pearson")$estimate)
  
  ## ESM1b
  LC_CDR3_only__CDR3_from_VDJ <- c(LC_CDR3_only__CDR3_from_VDJ, cor.test(x = df_sample$LC_evo_likelihood_esm1b_cdr3_from_VDJ,
                                                                         y = df_sample$LC_evo_likelihood_esm1b_cdr3_only,
                                                                         method = "pearson")$estimate)
  
  LC_full_VDJ__CDR3_from_VDJ <- c(LC_full_VDJ__CDR3_from_VDJ, cor.test(x = df_sample$LC_evo_likelihood_esm1b_cdr3_from_VDJ,
                                                                       y = df_sample$LC_evo_likelihood_esm1b_full_VDJ,
                                                                       method = "pearson")$estimate)
  
  LC_full_VDJ__CDR3_only <- c(LC_full_VDJ__CDR3_only, cor.test(x = df_sample$LC_evo_likelihood_esm1b_full_VDJ,
                                                               y = df_sample$LC_evo_likelihood_esm1b_cdr3_only,
                                                               method = "pearson")$estimate)
  
  ## ESMc
  LC_CDR3_only__CDR3_from_VDJ <- c(LC_CDR3_only__CDR3_from_VDJ, cor.test(x = df_sample$LC_evo_likelihood_esmc_cdr3_from_VDJ,
                                                                         y = df_sample$LC_evo_likelihood_esmc_cdr3_only,
                                                                         method = "pearson")$estimate)
  
  LC_full_VDJ__CDR3_from_VDJ <- c(LC_full_VDJ__CDR3_from_VDJ, cor.test(x = df_sample$LC_evo_likelihood_esmc_cdr3_from_VDJ,
                                                                       y = df_sample$LC_evo_likelihood_esmc_full_VDJ,
                                                                       method = "pearson")$estimate)
  
  LC_full_VDJ__CDR3_only <- c(LC_full_VDJ__CDR3_only, cor.test(x = df_sample$LC_evo_likelihood_esmc_full_VDJ,
                                                               y = df_sample$LC_evo_likelihood_esmc_cdr3_only,
                                                               method = "pearson")$estimate)
  
  ## ProtBERT
  LC_CDR3_only__CDR3_from_VDJ <- c(LC_CDR3_only__CDR3_from_VDJ, cor.test(x = df_sample$LC_evo_likelihood_protbert_cdr3_from_VDJ,
                                                                         y = df_sample$LC_evo_likelihood_protbert_cdr3_only,
                                                                         method = "pearson")$estimate)
  
  LC_full_VDJ__CDR3_from_VDJ <- c(LC_full_VDJ__CDR3_from_VDJ, cor.test(x = df_sample$LC_evo_likelihood_protbert_cdr3_from_VDJ,
                                                                       y = df_sample$LC_evo_likelihood_protbert_full_VDJ,
                                                                       method = "pearson")$estimate)
  
  LC_full_VDJ__CDR3_only <- c(LC_full_VDJ__CDR3_only, cor.test(x = df_sample$LC_evo_likelihood_protbert_full_VDJ,
                                                               y = df_sample$LC_evo_likelihood_protbert_cdr3_only,
                                                               method = "pearson")$estimate)
  
  #Combine into a dataframe
  LC_cor_df_sample <- data.frame(LC_full_VDJ__CDR3_only, LC_full_VDJ__CDR3_from_VDJ, LC_CDR3_only__CDR3_from_VDJ, 
                                 sample, chain = "LC", model = c("Ablang1", "Ablang2", "Sapiens", "ESM-1b", "ESM-C", "ProtBERT"))
  
  
  #Append to the main dataframe
  HC_cor_df <- rbind(HC_cor_df, HC_cor_df_sample)
  LC_cor_df <- rbind(LC_cor_df, LC_cor_df_sample)

}

write.csv(HC_cor_df, file = paste0("../data/",dataset,"/SourceCorrelation_HC.csv"), row.names = F)
write.csv(LC_cor_df, file = paste0("../data/",dataset,"/SourceCorrelation_LC.csv"), row.names = F)
