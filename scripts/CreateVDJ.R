#This script is used to build the VDJ dataframe for a given dataset and count somatic hypermutation

library(Platypus)
library(Biostrings)
library(stringdist)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
dataset <- args[1]

#List the sample
samples <- list.dirs(path = paste0("../data/",dataset,"/VDJ"), full.names = TRUE, recursive = FALSE)
print(samples)
#Build the VDJ dataframe
vdj <- VDJ_build(VDJ.sample.list = samples,
                 remove.divergent.cells = T,
                 complete.cells.only = T,
                 trim.germlines = T,
                 parallel = T,
                 num.cores = 5)
print("finished VDJ: ")

#Count somatic hypermutation (hamming distance, ignore gaps)
source("SHM_functions.R")
#Remove NA germlines
vdj <- vdj[is.na(vdj$VDJ_germline_nt_trimmed)==FALSE,]
vdj <- SHM_calculator(vdj)

save(vdj, file = paste0("../data/",dataset,"/VDJ_",dataset,".RData"))
write.csv(vdj, file = paste0("../data/",dataset,"/VDJ_",dataset,".csv"), row.names = FALSE)