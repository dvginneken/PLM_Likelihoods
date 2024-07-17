library(dplyr)

args = commandArgs(trailingOnly=TRUE)
dataset <- args[1]

#Read the data
df <- read.csv(file = paste0("../data/",dataset,"/evo_likelihoods_all.csv"))

#Combine kappa and lampda light chains
df$chain <- case_match(df$chain,
                       "IGH" ~ "HC",
                       "IGK" ~ "LC",
                       "IGL" ~ "LC")

#List the sample
samples <- list.dirs(path = paste0("../data/",dataset,"/VDJ"), full.names = FALSE, recursive = FALSE)

cor_df <- c()
for(sample in samples){
  #Subset the dataframe per sample
  df_sample <- df[df$original_sample_id == sample,]
  
  for(chain in c("HC", "LC")){
    
    #Subset the dataframe per chain
    df_chain <- df_sample[df_sample$chain == chain,]
    
    full_VDJ__CDR3_only = c()
    full_VDJ__CDR3_from_VDJ = c()
    CDR3_only__CDR3_from_VDJ = c()
    
    
    ## Ablang
    CDR3_only__CDR3_from_VDJ <- c(CDR3_only__CDR3_from_VDJ, cor.test(x = df_chain$evo_likelihood_ablang_cdr3_from_VDJ,
                                                                     y = df_chain$evo_likelihood_ablang_cdr3_only,
                                                                     method = "pearson")$estimate)
    
    #signifanctly correlated
    full_VDJ__CDR3_from_VDJ <- c(full_VDJ__CDR3_from_VDJ, cor.test(x = df_chain$evo_likelihood_ablang_cdr3_from_VDJ,
                                                                   y = df_chain$evo_likelihood_ablang_full_VDJ,
                                                                   method = "pearson")$estimate)
    
    full_VDJ__CDR3_only <- c(full_VDJ__CDR3_only, cor.test(x = df_chain$evo_likelihood_ablang_full_VDJ,
                                                           y = df_chain$evo_likelihood_ablang_cdr3_only,
                                                           method = "pearson")$estimate)
    
    ## Sapiens
    CDR3_only__CDR3_from_VDJ <- c(CDR3_only__CDR3_from_VDJ, cor.test(x = df_chain$evo_likelihood_sapiens_cdr3_from_VDJ,
                                                                     y = df_chain$evo_likelihood_sapiens_cdr3_only,
                                                                     method = "pearson")$estimate)
    
    #signifanctly correlated
    full_VDJ__CDR3_from_VDJ <- c(full_VDJ__CDR3_from_VDJ, cor.test(x = df_chain$evo_likelihood_sapiens_cdr3_from_VDJ,
                                                                   y = df_chain$evo_likelihood_sapiens_full_VDJ,
                                                                   method = "pearson")$estimate)
    
    full_VDJ__CDR3_only <- c(full_VDJ__CDR3_only, cor.test(x = df_chain$evo_likelihood_sapiens_full_VDJ,
                                                           y = df_chain$evo_likelihood_sapiens_cdr3_only,
                                                           method = "pearson")$estimate)
    
    ## ESM
    CDR3_only__CDR3_from_VDJ <- c(CDR3_only__CDR3_from_VDJ, cor.test(x = df_chain$evo_likelihood_esm_cdr3_from_VDJ,
                                                                     y = df_chain$evo_likelihood_esm_cdr3_only,
                                                                     method = "pearson")$estimate)
    
    full_VDJ__CDR3_from_VDJ <- c(full_VDJ__CDR3_from_VDJ, cor.test(x = df_chain$evo_likelihood_esm_cdr3_from_VDJ,
                                                                   y = df_chain$evo_likelihood_esm_full_VDJ,
                                                                   method = "pearson")$estimate)
    
    full_VDJ__CDR3_only <- c(full_VDJ__CDR3_only, cor.test(x = df_chain$evo_likelihood_esm_full_VDJ,
                                                           y = df_chain$evo_likelihood_esm_cdr3_only,
                                                           method = "pearson")$estimate)
    
    ## ProtBERT
    CDR3_only__CDR3_from_VDJ <- c(CDR3_only__CDR3_from_VDJ, cor.test(x = df_chain$evo_likelihood_protbert_cdr3_from_VDJ,
                                                                     y = df_chain$evo_likelihood_protbert_cdr3_only,
                                                                     method = "pearson")$estimate)
    
    full_VDJ__CDR3_from_VDJ <- c(full_VDJ__CDR3_from_VDJ, cor.test(x = df_chain$evo_likelihood_protbert_cdr3_from_VDJ,
                                                                   y = df_chain$evo_likelihood_protbert_full_VDJ,
                                                                   method = "pearson")$estimate)
    
    full_VDJ__CDR3_only <- c(full_VDJ__CDR3_only, cor.test(x = df_chain$evo_likelihood_protbert_full_VDJ,
                                                           y = df_chain$evo_likelihood_protbert_cdr3_only,
                                                           method = "pearson")$estimate)
    
    #Combine into a dataframe
    cor_df_sample <- data.frame(full_VDJ__CDR3_only, full_VDJ__CDR3_from_VDJ, CDR3_only__CDR3_from_VDJ, 
                                sample, chain, model = c("Ablang", "Sapiens", "ESM-1b", "ProtBERT"))
    
    #Append to the main dataframe
    cor_df <- rbind(cor_df, cor_df_sample)
  }

}

write.csv(cor_df, file = paste0("../data/",dataset,"/SourceCorrelation_chains.csv"), row.names = F)
