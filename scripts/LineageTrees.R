#This script is used to generate lineage trees using AntibodyForests
library(stringr)

args = commandArgs(trailingOnly=TRUE)
dataset <- args[1]
method <- args[2]

load(paste0("../data/",dataset,"/VDJ_PLL_",dataset,".RData")) #vdj

#Source AntibodyForests from https://github.com/alexyermanos/AntibodyForests
source("AntibodyForests.R")
af <- AntibodyForests(VDJ = vdj,
                            sequence.columns = "VDJ_sequence_aa_trimmed",
                            germline.columns = "VDJ_germline_aa_trimmed",
                            node.features = c("isotype", "VDJ_vgene", "SHM_count",
                                              "IGH_evo_likelihood_ablang_full_VDJ",
                                              "IGH_evo_likelihood_sapiens_full_VDJ",
                                              "IGH_evo_likelihood_protbert_full_VDJ",
                                              "IGH_evo_likelihood_esm_full_VDJ"),
                            construction.method = method,
                            parallel = F)

#Save AB-forests object
method = str_split_1(method, ".")[3]
save(af, file = paste0("../data/",dataset,"/AF_",dataset,"_",method,"_HC.RData"))