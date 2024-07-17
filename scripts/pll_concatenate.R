# This script concatenates the evo_likelihoods from all models (Ablang, Sapiens, ESM, Protbert) and all input sources (cdr3_from_VDJ, cdr3_only, full_VDJ) of all samples in a dataset into a single file.
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
dataset <- args[1]

samples <- list.dirs(path = paste0("../data/",dataset,"/VDJ"), full.names = FALSE, recursive = FALSE)
evo_likelihood <- c()
#for each sample
for (sample in samples){
  #ablang
  evo_likelihood_ablang_cdr3_from_VDJ <- read.csv(paste0("../data/",dataset,"/VDJ/",sample,"/evo_likelihoods/cdr3_from_VDJ/evo_likelihood_ablang.csv"), header = TRUE)
  evo_likelihood_ablang_cdr3_only <- read.csv(paste0("../data/",dataset,"/VDJ/",sample,"/evo_likelihoods/cdr3_only/evo_likelihood_ablang.csv"), header = TRUE)
  evo_likelihood_ablang_full_VDJ <- read.csv(paste0("../data/",dataset,"/VDJ/",sample,"/evo_likelihoods/full_VDJ/evo_likelihood_ablang.csv"), header = TRUE)
  
  #sapiens
  evo_likelihood_sapiens_cdr3_from_VDJ <- read.csv(paste0("../data/",dataset,"/VDJ/",sample,"/evo_likelihoods/cdr3_from_VDJ/evo_likelihood_sapiens.csv"), header = TRUE)
  evo_likelihood_sapiens_cdr3_only <- read.csv(paste0("../data/",dataset,"/VDJ/",sample,"/evo_likelihoods/cdr3_only/evo_likelihood_sapiens.csv"), header = TRUE)
  evo_likelihood_sapiens_full_VDJ <- read.csv(paste0("../data/",dataset,"/VDJ/",sample,"/evo_likelihoods/full_VDJ/evo_likelihood_sapiens.csv"), header = TRUE)
  
  #esm
  evo_likelihood_esm_cdr3_from_VDJ <- read.csv(paste0("../data/",dataset,"/VDJ/",sample,"/evo_likelihoods/cdr3_from_VDJ/evo_likelihood_esm.csv"), header = TRUE)
  evo_likelihood_esm_cdr3_only <- read.csv(paste0("../data/",dataset,"/VDJ/",sample,"/evo_likelihoods/cdr3_only/evo_likelihood_esm.csv"), header = TRUE)
  evo_likelihood_esm_full_VDJ <- read.csv(paste0("../data/",dataset,"/VDJ/",sample,"/evo_likelihoods/full_VDJ/evo_likelihood_esm.csv"), header = TRUE)
  
  #protbert
  evo_likelihood_protbert_cdr3_from_VDJ <- read.csv(paste0("../data/",dataset,"/VDJ/",sample,"/evo_likelihoods/cdr3_from_VDJ/evo_likelihood_protbert.csv"), header = TRUE)
  evo_likelihood_protbert_cdr3_only <- read.csv(paste0("../data/",dataset,"/VDJ/",sample,"/evo_likelihoods/cdr3_only/evo_likelihood_protbert.csv"), header = TRUE)
  evo_likelihood_protbert_full_VDJ <- read.csv(paste0("../data/",dataset,"/VDJ/",sample,"/evo_likelihoods/full_VDJ/evo_likelihood_protbert.csv"), header = TRUE)
  
  evo_likelihood_sample <- inner_join(evo_likelihood_ablang_cdr3_from_VDJ,
                                 evo_likelihood_ablang_cdr3_only[,c("contig_id","evo_likelihood")], by = "contig_id") %>%
                                 inner_join(evo_likelihood_ablang_full_VDJ[,c("contig_id","evo_likelihood")], by = "contig_id") %>%
                                 inner_join(evo_likelihood_sapiens_cdr3_from_VDJ[,c("contig_id","evo_likelihood")], by = "contig_id") %>%
                                 inner_join(evo_likelihood_sapiens_cdr3_only[,c("contig_id","evo_likelihood")], by = "contig_id") %>%
                                 inner_join(evo_likelihood_sapiens_full_VDJ[,c("contig_id","evo_likelihood")], by = "contig_id") %>%
                                 inner_join(evo_likelihood_esm_cdr3_from_VDJ[,c("contig_id","evo_likelihood")], by = "contig_id") %>%
                                 inner_join(evo_likelihood_esm_cdr3_only[,c("contig_id","evo_likelihood")], by = "contig_id") %>%
                                 inner_join(evo_likelihood_esm_full_VDJ[,c("contig_id","evo_likelihood")], by = "contig_id") %>%
                                 inner_join(evo_likelihood_protbert_cdr3_from_VDJ[,c("contig_id","evo_likelihood")], by = "contig_id") %>%
                                 inner_join(evo_likelihood_protbert_cdr3_only[,c("contig_id","evo_likelihood")], by = "contig_id") %>%
                                 inner_join(evo_likelihood_protbert_full_VDJ[,c("contig_id","evo_likelihood")], by = "contig_id")
  evo_likelihood_sample <- cbind(evo_likelihood_sample, sample)
  
  colnames(evo_likelihood_sample) <- c("barcode",
                                       "contig_id",
                                       "chain",
                                       "v_gene",
                                       "d_gene",
                                       "j_gene",
                                       "c_gene",
                                       "raw_clonotype_id",
                                       "raw_consensus_id",
                                       "evo_likelihood_ablang_cdr3_from_VDJ",
                                       "evo_likelihood_ablang_cdr3_only",
                                       "evo_likelihood_ablang_full_VDJ",
                                       "evo_likelihood_sapiens_cdr3_from_VDJ",
                                       "evo_likelihood_sapiens_cdr3_only",
                                       "evo_likelihood_sapiens_full_VDJ",
                                       "evo_likelihood_esm_cdr3_from_VDJ",
                                       "evo_likelihood_esm_cdr3_only",
                                       "evo_likelihood_esm_full_VDJ",
                                       "evo_likelihood_protbert_cdr3_from_VDJ",
                                       "evo_likelihood_protbert_cdr3_only",
                                       "evo_likelihood_protbert_full_VDJ",
                                       "original_sample_id")
  
  evo_likelihood <- rbind(evo_likelihood, evo_likelihood_sample)
}

write.csv(evo_likelihood, paste0("../data/",dataset,"/evo_likelihoods_all.csv"), row.names = FALSE)