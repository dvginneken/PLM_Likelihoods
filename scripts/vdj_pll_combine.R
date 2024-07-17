library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
dataset <- args[1]

#Read the data
df <- read.csv(file = paste0("../data/",dataset,"/evo_likelihoods_all.csv"))
load(paste0("../data/",dataset,"/VDJ_",dataset,".RData")) #vdj

#Check for overlapping barcodes and remove these cells
duplicate_bc <- vdj[duplicated(vdj$barcode),"barcode"]
paste0("Removed ",nrow(vdj[vdj$barcode %in% duplicate_bc,])," cells from VDJ due to duplicated barcode.")
vdj <- vdj[!(vdj$barcode %in% duplicate_bc),]

df$barcode <- substr(df$barcode, 1, nchar(df$barcode)-2)

for(chain in c("IGH", "IGK", "IGL")){
  df_chain <- df[df$chain == chain,]
  
  #Check for overlapping barcodes and remove these cells
  duplicate_bc <- df_chain[duplicated(df_chain$barcode),"barcode"]
  paste0("Removed ",nrow(df_chain[df_chain$barcode %in% duplicate_bc,])," cells from evo-likelihood dataframe due to duplicated barcode.")
  df_chain <- df_chain[!(df_chain$barcode %in% duplicate_bc),]
  
  #Join the vdj and evo-likelihoods
  vdj <- left_join(vdj, df_chain[,grepl("^evo|barcode", names(df_chain))], by = "barcode") %>% 
    rename_at(vars(starts_with("evo")), ~ str_c(paste0(chain,"_"), .))
  
}

save(vdj, file = paste0("../data/",dataset,"/VDJ_PLL_",dataset,".RData"))
write.csv(vdj, paste0("../data/",dataset,"/vdj_evolike_combine.csv"), row.names=F)
