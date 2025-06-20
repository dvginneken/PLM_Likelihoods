library(dplyr)

args = commandArgs(trailingOnly=TRUE)
dataset <- args[1]

#Read the data
load(paste0("../data/",dataset,"/VDJ_PLL_",dataset,".RData")) #vdj_likelihood
vdj <- vdj_likelihood

#List the sample
samples <- list.dirs(path = paste0("../data/",dataset,"/VDJ"), full.names = FALSE, recursive = FALSE)

cor_vdj <- c()
for(sample in samples){
  #Subset the dataframe per sample
  vdj_sample <- vdj[vdj$sample_id == sample,]

  full_VDJ__CDR3_only = c()
  full_VDJ__CDR3_from_VDJ = c()
  CDR3_only__CDR3_from_VDJ = c()
  
  
  ## Ablang1
  CDR3_only__CDR3_from_VDJ <- c(CDR3_only__CDR3_from_VDJ, cor.test(x = c(vdj_sample$HC_evo_likelihood_ablang1_cdr3_from_VDJ,
                                                                         vdj_sample$LC_evo_likelihood_ablang1_cdr3_from_VDJ),
                                                                   y = c(vdj_sample$HC_evo_likelihood_ablang1_cdr3_only,
                                                                         vdj_sample$LC_evo_likelihood_ablang1_cdr3_only),
                                                                   method = "pearson")$estimate)
  
  full_VDJ__CDR3_from_VDJ <- c(full_VDJ__CDR3_from_VDJ, cor.test(x = c(vdj_sample$HC_evo_likelihood_ablang1_cdr3_from_VDJ,
                                                                       vdj_sample$LC_evo_likelihood_ablang1_cdr3_from_VDJ),
                                                                 y = c(vdj_sample$HC_evo_likelihood_ablang1_full_VDJ,
                                                                       vdj_sample$LC_evo_likelihood_ablang1_full_VDJ),
                                                                 method = "pearson")$estimate)
  
  full_VDJ__CDR3_only <- c(full_VDJ__CDR3_only, cor.test(x = c(vdj_sample$HC_evo_likelihood_ablang1_full_VDJ,
                                                               vdj_sample$LC_evo_likelihood_ablang1_full_VDJ),
                                                         y = c(vdj_sample$HC_evo_likelihood_ablang1_cdr3_only,
                                                               vdj_sample$LC_evo_likelihood_ablang1_cdr3_only),
                                                         method = "pearson")$estimate)
  
  ## Ablang2
  CDR3_only__CDR3_from_VDJ <- c(CDR3_only__CDR3_from_VDJ, cor.test(x = c(vdj_sample$HC_evo_likelihood_ablang2_cdr3_from_VDJ,
                                                                         vdj_sample$LC_evo_likelihood_ablang2_cdr3_from_VDJ),
                                                                   y = c(vdj_sample$HC_evo_likelihood_ablang2_cdr3_only,
                                                                         vdj_sample$LC_evo_likelihood_ablang2_cdr3_only),
                                                                   method = "pearson")$estimate)
  
  full_VDJ__CDR3_from_VDJ <- c(full_VDJ__CDR3_from_VDJ, cor.test(x = c(vdj_sample$HC_evo_likelihood_ablang2_cdr3_from_VDJ,
                                                                       vdj_sample$LC_evo_likelihood_ablang2_cdr3_from_VDJ),
                                                                 y = c(vdj_sample$HC_evo_likelihood_ablang2_full_VDJ,
                                                                       vdj_sample$LC_evo_likelihood_ablang2_full_VDJ),
                                                                 method = "pearson")$estimate)
  
  full_VDJ__CDR3_only <- c(full_VDJ__CDR3_only, cor.test(x = c(vdj_sample$HC_evo_likelihood_ablang2_full_VDJ,
                                                               vdj_sample$LC_evo_likelihood_ablang2_full_VDJ),
                                                         y = c(vdj_sample$HC_evo_likelihood_ablang2_cdr3_only,
                                                               vdj_sample$LC_evo_likelihood_ablang2_cdr3_only),
                                                         method = "pearson")$estimate)
  
  ## Sapiens
  CDR3_only__CDR3_from_VDJ <- c(CDR3_only__CDR3_from_VDJ, cor.test(x = c(vdj_sample$HC_evo_likelihood_sapiens_cdr3_from_VDJ,
                                                                         vdj_sample$LC_evo_likelihood_sapiens_cdr3_from_VDJ),
                                                                   y = c(vdj_sample$HC_evo_likelihood_sapiens_cdr3_only,
                                                                         vdj_sample$LC_evo_likelihood_sapiens_cdr3_only),
                                                                   method = "pearson")$estimate)
  
  full_VDJ__CDR3_from_VDJ <- c(full_VDJ__CDR3_from_VDJ, cor.test(x = c(vdj_sample$HC_evo_likelihood_sapiens_cdr3_from_VDJ,
                                                                       vdj_sample$LC_evo_likelihood_sapiens_cdr3_from_VDJ),
                                                                 y = c(vdj_sample$HC_evo_likelihood_sapiens_full_VDJ,
                                                                       vdj_sample$LC_evo_likelihood_sapiens_full_VDJ),
                                                                 method = "pearson")$estimate)
  
  full_VDJ__CDR3_only <- c(full_VDJ__CDR3_only, cor.test(x = c(vdj_sample$HC_evo_likelihood_sapiens_full_VDJ,
                                                               vdj_sample$LC_evo_likelihood_sapiens_full_VDJ),
                                                         y = c(vdj_sample$HC_evo_likelihood_sapiens_cdr3_only,
                                                               vdj_sample$LC_evo_likelihood_sapiens_cdr3_only),
                                                         method = "pearson")$estimate)
  
  ## ESM 1b
  CDR3_only__CDR3_from_VDJ <- c(CDR3_only__CDR3_from_VDJ, cor.test(x = c(vdj_sample$HC_evo_likelihood_esm1b_cdr3_from_VDJ,
                                                                         vdj_sample$LC_evo_likelihood_esm1b_cdr3_from_VDJ),
                                                                   y = c(vdj_sample$HC_evo_likelihood_esm1b_cdr3_only,
                                                                         vdj_sample$LC_evo_likelihood_esm1b_cdr3_only),
                                                                   method = "pearson")$estimate)
  
  full_VDJ__CDR3_from_VDJ <- c(full_VDJ__CDR3_from_VDJ, cor.test(x = c(vdj_sample$HC_evo_likelihood_esm1b_cdr3_from_VDJ,
                                                                       vdj_sample$LC_evo_likelihood_esm1b_cdr3_from_VDJ),
                                                                 y = c(vdj_sample$HC_evo_likelihood_esm1b_full_VDJ,
                                                                       vdj_sample$LC_evo_likelihood_esm1b_full_VDJ),
                                                                 method = "pearson")$estimate)
  
  full_VDJ__CDR3_only <- c(full_VDJ__CDR3_only, cor.test(x = c(vdj_sample$HC_evo_likelihood_esm1b_full_VDJ,
                                                               vdj_sample$LC_evo_likelihood_esm1b_full_VDJ),
                                                         y = c(vdj_sample$HC_evo_likelihood_esm1b_cdr3_only,
                                                               vdj_sample$LC_evo_likelihood_esm1b_cdr3_only),
                                                         method = "pearson")$estimate)
  
  ## ESM C
  CDR3_only__CDR3_from_VDJ <- c(CDR3_only__CDR3_from_VDJ, cor.test(x = c(vdj_sample$HC_evo_likelihood_esmc_cdr3_from_VDJ,
                                                                         vdj_sample$LC_evo_likelihood_esmc_cdr3_from_VDJ),
                                                                   y = c(vdj_sample$HC_evo_likelihood_esmc_cdr3_only,
                                                                         vdj_sample$LC_evo_likelihood_esmc_cdr3_only),
                                                                   method = "pearson")$estimate)
  
  full_VDJ__CDR3_from_VDJ <- c(full_VDJ__CDR3_from_VDJ, cor.test(x = c(vdj_sample$HC_evo_likelihood_esmc_cdr3_from_VDJ,
                                                                       vdj_sample$LC_evo_likelihood_esmc_cdr3_from_VDJ),
                                                                 y = c(vdj_sample$HC_evo_likelihood_esmc_full_VDJ,
                                                                       vdj_sample$LC_evo_likelihood_esmc_full_VDJ),
                                                                 method = "pearson")$estimate)
  
  full_VDJ__CDR3_only <- c(full_VDJ__CDR3_only, cor.test(x = c(vdj_sample$HC_evo_likelihood_esmc_full_VDJ,
                                                               vdj_sample$LC_evo_likelihood_esmc_full_VDJ),
                                                         y = c(vdj_sample$HC_evo_likelihood_esmc_cdr3_only,
                                                               vdj_sample$LC_evo_likelihood_esmc_cdr3_only),
                                                         method = "pearson")$estimate)
  
  ## ProtBERT
  CDR3_only__CDR3_from_VDJ <- c(CDR3_only__CDR3_from_VDJ, cor.test(x = c(vdj_sample$HC_evo_likelihood_protbert_cdr3_from_VDJ,
                                                                         vdj_sample$LC_evo_likelihood_protbert_cdr3_from_VDJ),
                                                                   y = c(vdj_sample$HC_evo_likelihood_protbert_cdr3_only,
                                                                         vdj_sample$LC_evo_likelihood_protbert_cdr3_only),
                                                                   method = "pearson")$estimate)
  
  full_VDJ__CDR3_from_VDJ <- c(full_VDJ__CDR3_from_VDJ, cor.test(x = c(vdj_sample$HC_evo_likelihood_protbert_cdr3_from_VDJ,
                                                                       vdj_sample$LC_evo_likelihood_protbert_cdr3_from_VDJ),
                                                                 y = c(vdj_sample$HC_evo_likelihood_protbert_full_VDJ,
                                                                       vdj_sample$LC_evo_likelihood_protbert_full_VDJ),
                                                                 method = "pearson")$estimate)
  
  full_VDJ__CDR3_only <- c(full_VDJ__CDR3_only, cor.test(x = c(vdj_sample$HC_evo_likelihood_protbert_full_VDJ,
                                                               vdj_sample$LC_evo_likelihood_protbert_full_VDJ),
                                                         y = c(vdj_sample$HC_evo_likelihood_protbert_cdr3_only,
                                                               vdj_sample$LC_evo_likelihood_protbert_cdr3_only),
                                                         method = "pearson")$estimate)
  
  #Combine into a dataframe
  cor_vdj_sample <- data.frame(full_VDJ__CDR3_only, full_VDJ__CDR3_from_VDJ, CDR3_only__CDR3_from_VDJ, 
                              sample, model = c("Ablang1", "Ablang2", "Sapiens", "ESM-1b", "ESM-C", "ProtBERT"))

  #Append to the main dataframe
  cor_vdj <- rbind(cor_vdj, cor_vdj_sample)

  
}
write.csv(cor_vdj, file = paste0("../data/",dataset,"/SourceCorrelation.csv"), row.names = F)
