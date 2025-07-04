#This script is used to generate lineage trees using AntibodyForests
library(stringr)
library(AntibodyForests)

args = commandArgs(trailingOnly=TRUE)
dataset <- args[1]
method <- args[2]

load(paste0("../data/",dataset,"/VDJ_",dataset,".RData")) #vdj

#Source AntibodyForests from https://github.com/alexyermanos/AntibodyForests
source("AntibodyForests.R")
af <- Af_build(VDJ = vdj,
               sequence.columns = c("VDJ_sequence_aa_trimmed"),
               germline.columns = c("VDJ_germline_aa_trimmed"),
               construction.method = method,
               parallel = F)

#Save AB-forests object
method = str_split_1(method, ".")[3]
save(af, file = paste0("../data/",dataset,"/AF_",dataset,"_",method,"_HC.RData"))