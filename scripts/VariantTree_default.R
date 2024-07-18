#This script is used to generate a lineage tree using AntibodyForests

library(Platypus)

vdj <- read.delim("../data/OVA_V7/evo_likelihoods/evo_likelihood_esm.csv", sep = ",", header = T)

#Source AntibodyForests from https://github.com/alexyermanos/AntibodyForests
source("AntibodyForests.R")
af <- AntibodyForests(VDJ = vdj,
                            sequence.columns = c("VDJ_sequence_aa","VJ_sequence_aa"),
                            germline.columns = c("VDJ_germline_aa","VJ_germline_aa"),
                            node.features = c("evo_likelihood", "octet.affinity"),
                            construction.method = "phylo.network.default",
                            parallel = F)
save(af, file = "../data/OVA_V7/AF_VariantTree.RData")