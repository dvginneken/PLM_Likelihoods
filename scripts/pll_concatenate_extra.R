# This script concatenates the evo_likelihoods from all models and all input sources (cdr3_from_VDJ, cdr3_only, full_VDJ) of all samples in a dataset into a single file.
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
dataset <- args[1]

samples <- list.dirs(path = paste0("../data/",dataset,"/VDJ"), full.names = FALSE, recursive = FALSE)
load(paste0("../data/",dataset,"/VDJ_",dataset,".RData"))
vdj_likelihood <- data.frame()

#for each sample
for (sample in samples){
  vdj_sample <- vdj[vdj$sample_id == sample,]
  #ablang1
  evo_likelihood_ablang1_full_VDJ <- read.csv(paste0("../data/",dataset,"/VDJ/",sample,"/evo_likelihoods/full_VDJ/evo_likelihood_ablang.csv"), header = TRUE)[,c("barcode", "HC_evo_likelihood", "LC_evo_likelihood")]
  colnames(evo_likelihood_ablang1_full_VDJ) <- c("barcode", "HC_evo_likelihood_ablang1_full_VDJ", "LC_evo_likelihood_ablang1_full_VDJ")
  
  #ablang2
  evo_likelihood_ablang2_full_VDJ <- read.csv(paste0("../data/",dataset,"/VDJ/",sample,"/evo_likelihoods/full_VDJ/evo_likelihood_ablang2.csv"), header = TRUE)[,c("barcode", "HC_evo_likelihood", "LC_evo_likelihood", "pair_evo_likelihood")]
  colnames(evo_likelihood_ablang2_full_VDJ) <- c("barcode", "HC_evo_likelihood_ablang2_full_VDJ", "LC_evo_likelihood_ablang2_full_VDJ", "paired_evo_likelihood_ablang2_full_VDJ")
  
  #sapiens
  evo_likelihood_sapiens_full_VDJ <- read.csv(paste0("../data/",dataset,"/VDJ/",sample,"/evo_likelihoods/full_VDJ/evo_likelihood_sapiens.csv"), header = TRUE)[,c("barcode", "HC_evo_likelihood", "LC_evo_likelihood")]
  colnames(evo_likelihood_sapiens_full_VDJ) <- c("barcode", "HC_evo_likelihood_sapiens_full_VDJ", "LC_evo_likelihood_sapiens_full_VDJ")
  
  #esm 1b
  evo_likelihood_esm1b_full_VDJ <- read.csv(paste0("../data/",dataset,"/VDJ/",sample,"/evo_likelihoods/full_VDJ/evo_likelihood_esm1b.csv"), header = TRUE)[,c("barcode", "HC_evo_likelihood", "LC_evo_likelihood")]
  colnames(evo_likelihood_esm1b_full_VDJ) <- c("barcode", "HC_evo_likelihood_esm1b_full_VDJ", "LC_evo_likelihood_esm1b_full_VDJ")
  
  #esm c
  evo_likelihood_esmc_full_VDJ <- read.csv(paste0("../data/",dataset,"/VDJ/",sample,"/evo_likelihoods/full_VDJ/evo_likelihood_esmc.csv"), header = TRUE)[,c("barcode", "HC_evo_likelihood", "LC_evo_likelihood")]
  colnames(evo_likelihood_esmc_full_VDJ) <- c("barcode", "HC_evo_likelihood_esmc_full_VDJ", "LC_evo_likelihood_esmc_full_VDJ")
  
  #protbert
  evo_likelihood_protbert_full_VDJ <- read.csv(paste0("../data/",dataset,"/VDJ/",sample,"/evo_likelihoods/full_VDJ/evo_likelihood_protbert.csv"), header = TRUE)[,c("barcode", "HC_evo_likelihood", "LC_evo_likelihood")]
  colnames(evo_likelihood_protbert_full_VDJ) <- c("barcode", "HC_evo_likelihood_protbert_full_VDJ", "LC_evo_likelihood_protbert_full_VDJ")
  
  # Join all likelihoods to the VDJ dataframe
  vdj_sample <- left_join(vdj_sample, evo_likelihood_ablang1_full_VDJ, by = "barcode")
  vdj_sample <- left_join(vdj_sample, evo_likelihood_ablang2_full_VDJ, by = "barcode")
  vdj_sample <- left_join(vdj_sample, evo_likelihood_sapiens_full_VDJ, by = "barcode")
  vdj_sample <- left_join(vdj_sample, evo_likelihood_esm1b_full_VDJ, by = "barcode")
  vdj_sample <- left_join(vdj_sample, evo_likelihood_esmc_full_VDJ, by = "barcode")
  vdj_sample <- left_join(vdj_sample, evo_likelihood_protbert_full_VDJ, by = "barcode")
  
  vdj_likelihood <- rbind(vdj_likelihood, vdj_sample)
  
}

save(vdj_likelihood, file = paste0("../data/",dataset,"/VDJ_PLL_",dataset,".RData"))
write.csv(vdj_likelihood, paste0("../data/",dataset,"/vdj_evolike_combine.csv"), row.names=F)
