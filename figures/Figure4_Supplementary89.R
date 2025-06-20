##Figure 5 - Phylogenetic analysis of clonal lineages
library(igraph)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(AntibodyForests)

source("~/OneDrive - UMC Utrecht/Documenten/Steropodon_theme.R")

#Figure 5A+B
#Read data
load("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/OVA_V7/AF_OVA_V7_default_HC.RData") #af
af_ova <- af
rm(af)
load("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/horns2020a__VDJ_RAW/AF_horns2020a__VDJ_RAW_default_HC.RData")
af_horns <- af
rm(af)
load("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Bruhn/AF_Bruhn_default_HC.RData")
af_bruhn <- af
rm(af)
load("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim/AF_Kim_default_HC.RData")
af_kim <- af
rm(af)

#Read VDJ dataframes
load("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/OVA_V7/VDJ_PLL_OVA_V7.RData")
vdj_ova <- vdj_likelihood
load("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/horns2020a__VDJ_RAW/VDJ_PLL_horns2020a__VDJ_RAW.RData")
vdj_horns <- vdj_likelihood
load("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Bruhn/VDJ_PLL_Bruhn.RData")
vdj_bruhn <- vdj_likelihood
load("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim/VDJ_PLL_Kim.RData")
vdj_kim <- vdj_likelihood
rm(vdj_likelihood)

#Add pseudolikelihoods
af_ova <- Af_add_node_feature(af_ova, vdj_ova, c("HC_evo_likelihood_esmc_full_VDJ",
                                                 "paired_evo_likelihood_ablang2_full_VDJ"))
af_kim <- Af_add_node_feature(af_kim, vdj_kim, c("HC_evo_likelihood_esmc_full_VDJ",
                                                 "paired_evo_likelihood_ablang2_full_VDJ"))
af_horns <- Af_add_node_feature(af_horns, vdj_horns, c("HC_evo_likelihood_esmc_full_VDJ",
                                                         "paired_evo_likelihood_ablang2_full_VDJ"))
af_bruhn <- Af_add_node_feature(af_bruhn, vdj_bruhn, c("HC_evo_likelihood_esmc_full_VDJ",
                                                         "paired_evo_likelihood_ablang2_full_VDJ"))
#A - Lineage Trees
#Plot Figure 5A
Af_plot_tree(af_ova,
                     sample = "S1",
                     clonotype = "clonotype4",
                     color.by = "HC_evo_likelihood_esmc_full_VDJ",
                     label.by = "size",
                     show.size.legend = F,
                     main.title = "Mouse 1",
                     sub.title = "Clonotype 4",
                     arrow.size = 0.5,
             output.file = "/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figure4/Mouse1Clonotype4_esmc.pdf")



#Plot Figure 5B
Af_plot_tree(af_kim,
                     sample = "SRR17729703",
                     clonotype = "clonotype1",
                     color.by = "HC_evo_likelihood_esmc_full_VDJ",
                     label.by = "size",
                     show.size.legend = F,
                     main.title = "Individual 3",
                     sub.title = "Clonotype 1",
             arrow.size = 0.5,
             output.file = "/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figure4/Human3Clonotype1_esmc.pdf")

#B - Distance to germline
#Function to get a dataframe with the distance to the germline for each node and the pseudolikelihoods
distance_germline <- function(input){
  #Go over each tree in the AntibodyForests object and create a metric dataframe
  metric_df <- lapply(seq_along(input), function(sample){
    lapply(seq_along(input[[sample]]), function(clonotype){
      sample_name <- names(input)[sample]
      clonotype_name <- names(input[[sample_name]])[clonotype]
      
      tree = input[[sample_name]][[clonotype_name]][["igraph"]]
      nodes = igraph::V(tree)[names(igraph::V(tree)) != "germline"]
      node_features = input[[sample_name]][[clonotype_name]][["nodes"]][names(input[[sample_name]][[clonotype_name]][["nodes"]]) != "germline"]
      
      #Get the total length of shortest paths between each node and the germline
      distance <- igraph::distances(tree, v = "germline", to = nodes, algorithm = "dijkstra",
                                    weights = edge_attr(tree)$edge.length)
      distance = t(as.data.frame(distance))
      
      
      #Evo likelihood
      esm1b = unlist(lapply(node_features,function(x){unique(x[["IGH_evo_likelihood_esm_full_VDJ"]])[!is.na(unique(x[["IGH_evo_likelihood_esm_full_VDJ"]]))]}))[rownames(distance)]
      protbert = unlist(lapply(node_features,function(x){unique(x[["IGH_evo_likelihood_protbert_full_VDJ"]])[!is.na(unique(x[["IGH_evo_likelihood_protbert_full_VDJ"]]))]}))[rownames(distance)]
      ablang1 = unlist(lapply(node_features,function(x){unique(x[["IGH_evo_likelihood_ablang_full_VDJ"]])[!is.na(unique(x[["IGH_evo_likelihood_ablang_full_VDJ"]]))]}))[rownames(distance)]
      sapiens = unlist(lapply(node_features,function(x){unique(x[["IGH_evo_likelihood_sapiens_full_VDJ"]])[!is.na(unique(x[["IGH_evo_likelihood_sapiens_full_VDJ"]]))]}))[rownames(distance)]
      esmc = unlist(lapply(node_features,function(x){unique(x[["HC_evo_likelihood_esmc_full_VDJ"]])[!is.na(unique(x[["HC_evo_likelihood_esmc_full_VDJ"]]))]}))[rownames(distance)]
      ablang2 = unlist(lapply(node_features,function(x){unique(x[["paired_evo_likelihood_ablang2_full_VDJ"]])[!is.na(unique(x[["paired_evo_likelihood_ablang2_full_VDJ"]]))]}))[rownames(distance)]
      
      distance = data.frame("distance" = distance, "sample" = sample_name, "clonotype" = clonotype_name, 
                            "node" = rownames(distance), "esm1b" = esm1b, "protbert" = protbert, 
                            "ablang1" = ablang1, "sapiens" = sapiens, "esmc" = esmc, "ablang2" = ablang2)
      return(distance)
    })
  })
  #Transform list into dataframe
  dfs <- lapply(metric_df, function(x){do.call(rbind, x)})
  df <- do.call(rbind, dfs)
  rownames(df) <- paste0(df$sample, "_", df$clonotype, "_", df$node)
  
  return(df)
}
#Get distance to germline for all samples
distance_ova <- distance_germline(af_ova)
distance_horns <- distance_germline(af_horns)
distance_bruhn <- distance_germline(af_bruhn)
distance_kim <- distance_germline(af_kim)

#change sample names
distance_ova$sample <- case_match(distance_ova$sample,
                                  "S1" ~ "Mouse1",
                                  "S2" ~ "Mouse2",
                                  "S3" ~ "Mouse3",
                                  "S4" ~ "Mouse4",
                                  "S5" ~ "Mouse5")
distance_horns$sample <- case_match(distance_horns$sample,
                                    "Influenza.vac.11.12.human.S1" ~ "Individual1",
                                    "Influenza.vac.11.12.human.S2" ~ "Human2",
                                    "Influenza.vac.11.12.human.S3" ~ "Human3",
                                    "Influenza.vac.11.12.human.S4" ~ "Human4")
#Only keep 1 replicate
distance_horns <- distance_horns[distance_horns$sample == "Individual1",]
distance_bruhn$sample <- "Individual2"
distance_kim$sample <- case_match(distance_kim$sample,
                                  "SRR17729703" ~ "Individual3",
                                  "SRR17729692" ~ "Individual4",
                                  "SRR17729726" ~ "Individual5")
#Only keep specific samples
distance_kim <- distance_kim[distance_kim$sample %in% paste0("Individual",3:5),]
distance_all <- rbind(distance_ova, distance_horns, distance_bruhn, distance_kim)
distance_human <- rbind(distance_horns, distance_bruhn, distance_kim)

#Plot the distance to the germline against the pseudolikelihoods
#esmc
cor_human <- cor.test(distance_human$germline, distance_human$esmc)$estimate
cor_mouse <- cor.test(distance_ova$germline, distance_ova$esmc)$estimate
p <- ggplot(distance_all, aes(x = germline, y = esmc, colour = sample)) + 
  geom_point(size = 1) + 
  geom_smooth(method = "lm", se=F) +
  scale_color_manual(values = c("Individual1" = "#99d8c9",
                                "Individual2" = "#66c2a4",
                                "Individual3" = "#41ae76",
                                "Individual4" = "#238b45",
                                "Individual5" = "#005824",
                                "Mouse1" = "#fcc5c0",
                                "Mouse2" = "#fa9fb5",
                                "Mouse3" = "#f768a1",
                                "Mouse4" = "#c51b8a",
                                "Mouse5" = "#7a0177"),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 18)) +
  xlab("Distance to germline") + ylab("ESM-C SP") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), "\nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#Plot Figure 4C
pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figure4/Distance_esmc.pdf")
print(p)
dev.off()

#ablang2
cor_human <- cor.test(distance_human$germline, distance_human$ablang2)$estimate
cor_mouse <- cor.test(distance_ova$germline, distance_ova$ablang2)$estimate
q <- ggplot(distance_all, aes(x = germline, y = ablang2, colour = sample)) + 
  geom_point(size = 1) + 
  geom_smooth(method = "lm", se=F) +
  scale_color_manual(values = c("Individual1" = "#99d8c9",
                                "Individual2" = "#66c2a4",
                                "Individual3" = "#41ae76",
                                "Individual4" = "#238b45",
                                "Individual5" = "#005824",
                                "Mouse1" = "#fcc5c0",
                                "Mouse2" = "#fa9fb5",
                                "Mouse3" = "#f768a1",
                                "Mouse4" = "#c51b8a",
                                "Mouse5" = "#7a0177"),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 18)) +
  xlab("Distance to germline") + ylab("Ablang2 SP") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), "\nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))

pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figure4/Distance_ablang2.pdf")
print(p)
dev.off()

pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figure4/Distance_main.pdf", width = 10, height = 5)
ggarrange(p, q, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()

#Supplementary Figure 6
#protbert
cor_human <- cor.test(distance_human$germline, distance_human$protbert)$estimate
cor_mouse <- cor.test(distance_ova$germline, distance_ova$protbert)$estimate
a <- ggplot(distance_all, aes(x = germline, y = protbert, colour = sample)) + 
  geom_point(size = 1) + 
  geom_smooth(method = "lm", se=F) +
  scale_color_manual(values = c("Individual1" = "#99d8c9",
                                "Individual2" = "#66c2a4",
                                "Individual3" = "#41ae76",
                                "Individual4" = "#238b45",
                                "Individual5" = "#005824",
                                "Mouse1" = "#fcc5c0",
                                "Mouse2" = "#fa9fb5",
                                "Mouse3" = "#f768a1",
                                "Mouse4" = "#c51b8a",
                                "Mouse5" = "#7a0177"),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 18)) +
  xlab("Distance to germline") + ylab("ProtBERT SP") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), "\nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))
#ablang
cor_human <- cor.test(distance_human$germline, distance_human$ablang1)$estimate
cor_mouse <- cor.test(distance_ova$germline, distance_ova$ablang1)$estimate
b <- ggplot(distance_all, aes(x = germline, y = ablang1, colour = sample)) + 
  geom_point(size = 1) + 
  geom_smooth(method = "lm", se=F) +
  scale_color_manual(values = c("Individual1" = "#99d8c9",
                                "Individual2" = "#66c2a4",
                                "Individual3" = "#41ae76",
                                "Individual4" = "#238b45",
                                "Individual5" = "#005824",
                                "Mouse1" = "#fcc5c0",
                                "Mouse2" = "#fa9fb5",
                                "Mouse3" = "#f768a1",
                                "Mouse4" = "#c51b8a",
                                "Mouse5" = "#7a0177"),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 18)) +
  xlab("Distance to germline") + ylab("Ablang SP") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), "\nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))
#sapiens
cor_human <- cor.test(distance_human$germline, distance_human$sapiens)$estimate
cor_mouse <- cor.test(distance_ova$germline, distance_ova$sapiens)$estimate
c <- ggplot(distance_all, aes(x = germline, y = sapiens, colour = sample)) + 
  geom_point(size = 1) + 
  geom_smooth(method = "lm", se=F) +
  scale_color_manual(values = c("Individual1" = "#99d8c9",
                                "Individual2" = "#66c2a4",
                                "Individual3" = "#41ae76",
                                "Individual4" = "#238b45",
                                "Individual5" = "#005824",
                                "Mouse1" = "#fcc5c0",
                                "Mouse2" = "#fa9fb5",
                                "Mouse3" = "#f768a1",
                                "Mouse4" = "#c51b8a",
                                "Mouse5" = "#7a0177"),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 18)) +
  xlab("Distance to germline") + ylab("Sapiens SP") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), "\nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))
#esm-1b
cor_human <- cor.test(distance_human$germline, distance_human$esm1b)$estimate
cor_mouse <- cor.test(distance_ova$germline, distance_ova$esm1b)$estimate
d <- ggplot(distance_all, aes(x = germline, y = esm1b, colour = sample)) + 
  geom_point(size = 1) + 
  geom_smooth(method = "lm", se=F) +
  scale_color_manual(values = c("Individual1" = "#99d8c9",
                                "Individual2" = "#66c2a4",
                                "Individual3" = "#41ae76",
                                "Individual4" = "#238b45",
                                "Individual5" = "#005824",
                                "Mouse1" = "#fcc5c0",
                                "Mouse2" = "#fa9fb5",
                                "Mouse3" = "#f768a1",
                                "Mouse4" = "#c51b8a",
                                "Mouse5" = "#7a0177"),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 18)) +
  xlab("Distance to germline") + ylab("ESM-1b SP") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), "\nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#Plot supplementary Figure 6
pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figure4/DistanceGermline_sup.pdf", width = 20, height = 6)
ggarrange(a, b, c, d, ncol = 4, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()

#Figure 4DE
#Read data
mut_ova <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/OVA_V7/mutational_rank_table_OVA_V7_default.csv")
mut_ova2 <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/OVA_V7/mutational_rank_table_OVA_V7_default_2.csv")
mut_ova <- rbind(mut_ova, mut_ova2)
mut_horns <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/horns2020a__VDJ_RAW/mutational_rank_table_horns2020a__VDJ_RAW_default.csv")
mut_horns2 <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/horns2020a__VDJ_RAW/mutational_rank_table_horns2020a__VDJ_RAW_default_2.csv")
mut_horns <- rbind(mut_horns, mut_horns2)
mut_bruhn <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Bruhn/mutational_rank_table_Bruhn_default.csv")
mut_bruhn2 <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Bruhn/mutational_rank_table_Bruhn_default_2.csv")
mut_bruhn <- rbind(mut_bruhn, mut_bruhn2)
mut_kim <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim/mutational_rank_table_Kim_default.csv")
mut_kim2 <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim/mutational_rank_table_Kim_default_2.csv")
mut_kim <- rbind(mut_kim, mut_kim2)
mut_kim_extra <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim_extra/mutational_rank_table_Kim_extra_default.csv")
mut_kim_extra2 <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim_extra/mutational_rank_table_Kim_extra_default_2.csv")
mut_kim_extra <- rbind(mut_kim_extra, mut_kim_extra2)
mut_mathew <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Mathew/mutational_rank_table_Mathew_default.csv")
mut_mathew2 <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Mathew/mutational_rank_table_Mathew_default_2.csv")
mut_mathew <- rbind(mut_mathew, mut_mathew2)

#change sample names
mut_ova$sample <- case_match(mut_ova$sample,
                             "S1" ~ "Mouse1",
                             "S2" ~ "Mouse2",
                             "S3" ~ "Mouse3",
                             "S4" ~ "Mouse4",
                             "S5" ~ "Mouse5")
mut_horns$sample <- case_match(mut_horns$sample,
                               "Influenza.vac.11.12.human.S1" ~ "Individual1",
                               "Influenza.vac.11.12.human.S2" ~ "Human2",
                               "Influenza.vac.11.12.human.S3" ~ "Human3",
                               "Influenza.vac.11.12.human.S4" ~ "Human4")
#Only keep 1 replicate
mut_horns <- mut_horns[mut_horns$sample == "Individual1",]
mut_bruhn$sample <- "Individual2"
mut_kim$sample <- case_match(mut_kim$sample,
                             "SRR17729703" ~ "Individual3",
                             "SRR17729692" ~ "Individual4",
                             "SRR17729726" ~ "Individual5")
mut_kim <- mut_kim[mut_kim$sample %in% paste0("Individual",2:5),]
mut_human <- rbind(mut_horns, mut_bruhn, mut_kim)
mut_all <- rbind(mut_ova, mut_human)

mut_all_esmc <- mut_all[mut_all$model == "esmc",]
mut_all_ablang2 <- mut_all[mut_all$model == "ablang2",]
mut_all_esm1b <- mut_all[mut_all$model == "esm1b",]
mut_all_protbert <- mut_all[mut_all$model == "protbert",]
mut_all_sapiens <- mut_all[mut_all$model == "sapiens",]
mut_all_ablang1 <- mut_all[mut_all$model == "ablang1",]

# C - Distribution of likelihood rank of the substitution along the edge
a <- ggplot(mut_all_esmc, aes(mean_sub_rank)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 1, linewidth = 1) +
  scale_color_manual(values = c("Individual1" = "#99d8c9",
                                "Individual2" = "#66c2a4",
                                "Individual3" = "#41ae76",
                                "Individual4" = "#238b45",
                                "Individual5" = "#005824",
                                "Mouse1" = "#fcc5c0",
                                "Mouse2" = "#fa9fb5",
                                "Mouse3" = "#f768a1",
                                "Mouse4" = "#c51b8a",
                                "Mouse5" = "#7a0177"),
                     name = "Sample") +
  theme_steropodon() +
  ylim(0,300)+
  scale_x_continuous(breaks = c(1,5,10,15,20), limits = c(1,NA)) +
  theme(text = element_text(size = 16)) +
  xlab("Substitution Rank") + ylab("Number of edges") +
  ggtitle("ESM-C")
b <- ggplot(mut_all_esm1b, aes(mean_sub_rank)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 1, linewidth = 1) +
  scale_color_manual(values = c("Individual1" = "#99d8c9",
                                "Individual2" = "#66c2a4",
                                "Individual3" = "#41ae76",
                                "Individual4" = "#238b45",
                                "Individual5" = "#005824",
                                "Mouse1" = "#fcc5c0",
                                "Mouse2" = "#fa9fb5",
                                "Mouse3" = "#f768a1",
                                "Mouse4" = "#c51b8a",
                                "Mouse5" = "#7a0177"),
                     name = "Sample") +
  theme_steropodon() +
  ylim(0,300)+
  scale_x_continuous(breaks = c(1,5,10,15,20), limits = c(1,NA)) +
  theme(text = element_text(size = 16)) +
  xlab("Substitution Rank") + ylab("Number of edges") +
  ggtitle("ESM-1b")
c <- ggplot(mut_all_ablang2, aes(mean_sub_rank)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 1, linewidth = 1) +
  scale_color_manual(values = c("Individual1" = "#99d8c9",
                                "Individual2" = "#66c2a4",
                                "Individual3" = "#41ae76",
                                "Individual4" = "#238b45",
                                "Individual5" = "#005824",
                                "Mouse1" = "#fcc5c0",
                                "Mouse2" = "#fa9fb5",
                                "Mouse3" = "#f768a1",
                                "Mouse4" = "#c51b8a",
                                "Mouse5" = "#7a0177"),
                     name = "Sample") +
  theme_steropodon() +
  ylim(0,300)+
  scale_x_continuous(breaks = c(1,5,10,15,20), limits = c(1,NA)) +
  theme(text = element_text(size = 16)) +
  xlab("Substitution Rank") + ylab("Number of edges") +
  ggtitle("Ablang2")
d <- ggplot(mut_all_ablang1, aes(mean_sub_rank)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 1, linewidth = 1) +
  scale_color_manual(values = c("Individual1" = "#99d8c9",
                                "Individual2" = "#66c2a4",
                                "Individual3" = "#41ae76",
                                "Individual4" = "#238b45",
                                "Individual5" = "#005824",
                                "Mouse1" = "#fcc5c0",
                                "Mouse2" = "#fa9fb5",
                                "Mouse3" = "#f768a1",
                                "Mouse4" = "#c51b8a",
                                "Mouse5" = "#7a0177"),
                     name = "Sample") +
  theme_steropodon() +
  ylim(0,300)+
  scale_x_continuous(breaks = c(1,5,10,15,20), limits = c(1,NA)) +
  theme(text = element_text(size = 16)) +
  xlab("Substitution Rank") + ylab("Number of edges") +
  ggtitle("Ablang1")

#Plot Figure 5C
pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figure4/MutationRanks_main.pdf", width = 20, height = 4.5)
print(ggarrange(a, b, c, d, ncol = 4, common.legend = T, legend = "right"))
dev.off()

# ggplot(mut_all_esmc, aes(mean_sub_prob)) +
#   #geom_histogram(color = "white", fill = "black", binwidth = 1) +
#   geom_freqpoly(aes(colour = sample), linewidth = 2) +
#   scale_color_manual(values = c("Individual1" = "#99d8c9",
#                                 "Individual2" = "#66c2a4",
#                                 "Individual3" = "#41ae76",
#                                 "Individual4" = "#238b45",
#                                 "Individual5" = "#005824",
#                                 "Mouse1" = "#fcc5c0",
#                                 "Mouse2" = "#fa9fb5",
#                                 "Mouse3" = "#f768a1",
#                                 "Mouse4" = "#c51b8a",
#                                 "Mouse5" = "#7a0177"),
#                      name = "Sample") +
#   theme_minimal() +
#   theme(text = element_text(size = 16)) +
#   xlab("Substitution Likelihood") + ylab("Number of edges") +
#   ggtitle("ESM-C")
# 
# ggplot(mut_all_ablang2, aes(mean_sub_prob)) +
#   #geom_histogram(color = "white", fill = "black", binwidth = 1) +
#   geom_freqpoly(aes(colour = sample), linewidth = 2) +
#   scale_color_manual(values = c("Individual1" = "#99d8c9",
#                                 "Individual2" = "#66c2a4",
#                                 "Individual3" = "#41ae76",
#                                 "Individual4" = "#238b45",
#                                 "Individual5" = "#005824",
#                                 "Mouse1" = "#fcc5c0",
#                                 "Mouse2" = "#fa9fb5",
#                                 "Mouse3" = "#f768a1",
#                                 "Mouse4" = "#c51b8a",
#                                 "Mouse5" = "#7a0177"),
#                      name = "Sample") +
#   theme_minimal() +
#   theme(text = element_text(size = 16)) +
#   xlab("Substitution Likelihood") + ylab("Number of edges") +
#   ggtitle("Ablang2")

#E - quantification of mean likelihood rank per sample for each model
mut_all <- mut_all[mut_all$model != "",]
mut_all_esm1b %>% group_by(sample) %>% mutate(average = mean(as.numeric(mean_sub_rank), na.rm = T)) -> mut_all_esm1b
mut_all_ablang1 %>% group_by(sample) %>% mutate(average = mean(as.numeric(mean_sub_rank), na.rm = T)) -> mut_all_ablang1
mut_all_sapiens %>% group_by(sample) %>% mutate(average = mean(as.numeric(mean_sub_rank), na.rm = T)) -> mut_all_sapiens
mut_all_protbert %>% group_by(sample) %>% mutate(average = mean(as.numeric(mean_sub_rank), na.rm = T)) -> mut_all_protbert
mut_all_ablang2 %>% group_by(sample) %>% mutate(average = mean(as.numeric(mean_sub_rank), na.rm = T)) -> mut_all_ablang2
mut_all_esmc %>% group_by(sample) %>% mutate(average = mean(as.numeric(mean_sub_rank), na.rm = T)) -> mut_all_esmc
mut_all <- rbind(mut_all_esm1b, mut_all_ablang1, mut_all_sapiens, mut_all_protbert, mut_all_ablang2, mut_all_esmc)
df <- na.omit(unique(mut_all[,c("sample", "average", "model")]))
df2 <- df %>% group_by(model) %>% summarize(Average = mean(average)) %>% ungroup()
p <- ggplot(df, aes(y = average, x = model)) +
  geom_col(data = df2, aes(y = Average, x = model), color = "black", fill = "white") +
  geom_point(aes(color = sample), size = 3) +
  scale_color_manual(values = c("Individual1" = "#99d8c9",
                                "Individual2" = "#66c2a4",
                                "Individual3" = "#41ae76",
                                "Individual4" = "#238b45",
                                "Individual5" = "#005824",
                                "Mouse1" = "#fcc5c0",
                                "Mouse2" = "#fa9fb5",
                                "Mouse3" = "#f768a1",
                                "Mouse4" = "#c51b8a",
                                "Mouse5" = "#7a0177"),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 16)) +
  ylab("Likelihood Rank") + xlab("Model") +
  scale_x_discrete(labels = c("esm" = "ESM-1b", "ablang" = "Ablang", "protbert" = "ProtBERT", "sapiens" = "Sapiens")) +
  ggtitle("Average Substitution Rank")

#Plot Figure 5D
pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figure4/AverageSubstitutionRank.pdf", width = 7, height = 5.5)
print(p)
dev.off()

#Figure E 
#Read data convserved residues
conserved_ova <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/OVA_V7/conserved_rank_table_OVA_V7_default.csv")
conserved_ova2 <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/OVA_V7/conserved_rank_table_OVA_V7_default_2.csv")
conserved_ova <- rbind(conserved_ova, conserved_ova2)
conserved_horns <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/horns2020a__VDJ_RAW/conserved_rank_table_horns2020a__VDJ_RAW_default.csv")
conserved_horns2 <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/horns2020a__VDJ_RAW/conserved_rank_table_horns2020a__VDJ_RAW_default_2.csv")
conserved_horns <- rbind(conserved_horns, conserved_horns2)
conserved_bruhn <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Bruhn/conserved_rank_table_Bruhn_default.csv")
conserved_bruhn2 <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Bruhn/conserved_rank_table_Bruhn_default_2.csv")
conserved_bruhn <- rbind(conserved_bruhn, conserved_bruhn2)
conserved_kim <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim/conserved_rank_table_Kim_default.csv")
conserved_kim2 <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim/conserved_rank_table_Kim_default_2.csv")
conserved_kim <- rbind(conserved_kim, conserved_kim2)
conserved_kim_extra <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim_extra/conserved_rank_table_Kim_extra_default.csv")
conserved_kim_extra2 <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim_extra/conserved_rank_table_Kim_extra_default_2.csv")
conserved_kim_extra <- rbind(conserved_kim_extra, conserved_kim_extra2)
conserved_mathew <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Mathew/conserved_rank_table_Mathew_default.csv")
conserved_mathew2 <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Mathew/conserved_rank_table_Mathew_default_2.csv")
conserved_mathew <- rbind(conserved_mathew, conserved_mathew2)

#change sample names
conserved_ova$sample <- case_match(conserved_ova$sample,
                                   "S1" ~ "Mouse1",
                                   "S2" ~ "Mouse2",
                                   "S3" ~ "Mouse3",
                                   "S4" ~ "Mouse4",
                                   "S5" ~ "Mouse5")
conserved_horns$sample <- case_match(conserved_horns$sample,
                                     "Influenza.vac.11.12.human.S1" ~ "Individual1",
                                     "Influenza.vac.11.12.human.S2" ~ "Human2",
                                     "Influenza.vac.11.12.human.S3" ~ "Human3",
                                     "Influenza.vac.11.12.human.S4" ~ "Human4")
conserved_horns <- conserved_horns[conserved_horns$sample == "Individual1",]
conserved_bruhn$sample <- "Individual2"
conserved_kim$sample <- case_match(conserved_kim$sample,
                                   "SRR17729703" ~ "Individual3",
                                   "SRR17729692" ~ "Individual4",
                                   "SRR17729726" ~ "Individual5")
conserved_kim <- conserved_kim[conserved_kim$sample %in% paste0("Individual",2:5),]
conserved_human <- rbind(conserved_horns, conserved_bruhn, conserved_kim)
conserved_all <- rbind(conserved_ova, conserved_human)

#separate on model
conserved_all_esmc <- conserved_all[conserved_all$model == "esmc",]
conserved_all_protbert <- conserved_all[conserved_all$model == "protbert",]
conserved_all_sapiens <- conserved_all[conserved_all$model == "sapiens",]
conserved_all_ablang2 <- conserved_all[conserved_all$model == "ablang2",]
conserved_all_ablang1 <- conserved_all[conserved_all$model == "ablang1",]
conserved_all_esm1b <- conserved_all[conserved_all$model == "esm1b",]

#read data original mutated residues
original_ova <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/OVA_V7/mutational_rank_reversed_table_OVA_V7_default.csv")
original_ova2 <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/OVA_V7/mutational_rank_reversed_table_OVA_V7_default_2.csv")
original_ova <- rbind(original_ova, original_ova2)
original_horns <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/horns2020a__VDJ_RAW/mutational_rank_reversed_table_horns2020a__VDJ_RAW_default.csv")
original_horns2 <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/horns2020a__VDJ_RAW/mutational_rank_reversed_table_horns2020a__VDJ_RAW_default_2.csv")
original_horns <- rbind(original_horns, original_horns2)
original_bruhn <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Bruhn/mutational_rank_reversed_table_Bruhn_default.csv")
original_bruhn2 <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Bruhn/mutational_rank_reversed_table_Bruhn_default_2.csv")
original_bruhn <- rbind(original_bruhn, original_bruhn2)
original_kim <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim/mutational_rank_reversed_table_Kim_default.csv")
original_kim2 <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim/mutational_rank_reversed_table_Kim_default_2.csv")
original_kim <- rbind(original_kim, original_kim2)
original_kim_extra <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim_extra/mutational_rank_reversed_table_Kim_extra_default.csv")
original_kim_extra2 <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim_extra/mutational_rank_reversed_table_Kim_extra_default_2.csv")
original_kim_extra <- rbind(original_kim_extra, original_kim_extra2)
original_mathew <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Mathew/mutational_rank_reversed_table_Mathew_default.csv")
original_mathew2 <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Mathew/mutational_rank_reversed_table_Mathew_default_2.csv")
original_mathew <- rbind(original_mathew, original_mathew2)

#change sample names
original_ova$sample <- case_match(original_ova$sample,
                                  "S1" ~ "Mouse1",
                                  "S2" ~ "Mouse2",
                                  "S3" ~ "Mouse3",
                                  "S4" ~ "Mouse4",
                                  "S5" ~ "Mouse5")
original_horns$sample <- case_match(original_horns$sample,
                                    "Influenza.vac.11.12.human.S1" ~ "Individual1",
                                    "Influenza.vac.11.12.human.S2" ~ "Human2",
                                    "Influenza.vac.11.12.human.S3" ~ "Human3",
                                    "Influenza.vac.11.12.human.S4" ~ "Human4")
original_horns <- original_horns[original_horns$sample == "Individual1",]
original_bruhn$sample <- "Individual2"
original_kim$sample <- case_match(original_kim$sample,
                                  "SRR17729703" ~ "Individual3",
                                  "SRR17729692" ~ "Individual4",
                                  "SRR17729726" ~ "Individual5")
original_kim <- original_kim[original_kim$sample %in% paste0("Individual",2:5),]
original_human <- rbind(original_horns, original_bruhn, original_kim)
original_all <- rbind(original_ova, original_human)

#separate on model
original_all_esmc <- original_all[original_all$model == "esmc",]
original_all_protbert <- original_all[original_all$model == "protbert",]
original_all_sapiens <- original_all[original_all$model == "sapiens",]
original_all_ablang2 <- original_all[original_all$model == "ablang2",]
original_all_ablang1 <- original_all[original_all$model == "ablang1",]
original_all_esm1b <- original_all[original_all$model == "esm1b",]

#Figure E - Significance boxplot of likelihood ranks
#esm
conserved_all_esmc$group <- "conserved"
original_all_esmc$group <- "mutating"
df <- rbind(conserved_all_esmc[,-8], original_all_esmc[-8])

stat.test <- df %>%
  t_test(mean_sub_rank ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test
p <- ggplot(df, aes(y = as.numeric(mean_sub_rank), x = group, color=group)) +
  geom_boxplot() +
  scale_color_manual(values = c("conserved" = "orange",
                                "mutating" = "dodgerblue3")) +
  theme_minimal() +
  theme(text = element_text(size = 18),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("RL rank") +
  geom_signif(comparisons=list(c("conserved","mutating")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ylim(0,22) +
  ggtitle("ESM-C")

#ablang2
conserved_all_ablang2$group <- "conserved"
original_all_ablang2$group <- "mutating"
df <- rbind(conserved_all_ablang2[,-8], original_all_ablang2[-8])
stat.test <- df %>%
  t_test(mean_sub_rank ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test
o <- ggplot(df, aes(y = as.numeric(mean_sub_rank), x = group, color=group)) +
  geom_boxplot() +
  scale_color_manual(values = c("conserved" = "orange",
                                "mutating" = "dodgerblue3")) +
  theme_minimal() +
  theme(text = element_text(size = 18),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("RL rank") +
  geom_signif(comparisons=list(c("conserved","mutating")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ylim(0,22) +
  ggtitle("Ablang2")

#Plot Figure 4F
pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figure4/ConservedMutating_Rank.pdf", width = 7, height = 5.5)
ggarrange(p, o, nrow = 1, ncol = 2)
dev.off()


#Figure G - Significance boxplot of likelihood of the mutating and conserved residues
conserved_all_esmc$group <- "conserved"
original_all_esmc$group <- "mutating"
df <- rbind(conserved_all_esmc[,-8], original_all_esmc[-8])
stat.test <- df %>%
  t_test(mean_sub_prob ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test
p <- ggplot(df, aes(y = as.numeric(mean_sub_prob), x = group, color=group)) +
  geom_boxplot(outliers = F) +
  scale_color_manual(values = c("conserved" = "orange",
                                "mutating" = "dodgerblue3")) +
  theme_minimal() +
  theme(text = element_text(size = 18),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("Average RL") +
  geom_signif(comparisons=list(c("conserved","mutating")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ylim(0,1.05) +
  ggtitle("ESM-C")

conserved_all_ablang2$group <- "conserved"
original_all_ablang2$group <- "mutating"
df <- rbind(conserved_all_ablang2[,-8], original_all_ablang2[-8])
stat.test <- df %>%
  t_test(mean_sub_prob ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test
o <- ggplot(df, aes(y = as.numeric(mean_sub_prob), x = group, color=group)) +
  geom_boxplot(outliers = F) +
  scale_color_manual(values = c("conserved" = "orange",
                                "mutating" = "dodgerblue3")) +
  theme_minimal() +
  theme(text = element_text(size = 18),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("Average RL") +
  geom_signif(comparisons=list(c("conserved","mutating")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ylim(0,1.05) +
  ggtitle("Ablang2")



conserved_all_esm1b$group <- "conserved"
original_all_esm1b$group <- "mutating"
df <- rbind(conserved_all_esm1b[,-8], original_all_esm1b[-8])
stat.test <- df %>%
  t_test(mean_sub_prob ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test
esm1b <- ggplot(df, aes(y = as.numeric(mean_sub_prob), x = group, color=group)) +
  geom_boxplot(outliers = F) +
  scale_color_manual(values = c("conserved" = "orange",
                                "mutating" = "dodgerblue3")) +
  theme_minimal() +
  theme(text = element_text(size = 18),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("Average RL") +
  geom_signif(comparisons=list(c("conserved","mutating")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ylim(0,1.05) +
  ggtitle("ESM-1b")

conserved_all_ablang1$group <- "conserved"
original_all_ablang1$group <- "mutating"
df <- rbind(conserved_all_ablang1[,-8], original_all_ablang1[-8])
stat.test <- df %>%
  t_test(mean_sub_prob ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test
ab1 <- ggplot(df, aes(y = as.numeric(mean_sub_prob), x = group, color=group)) +
  geom_boxplot(outliers = F) +
  scale_color_manual(values = c("conserved" = "orange",
                                "mutating" = "dodgerblue3")) +
  theme_minimal() +
  theme(text = element_text(size = 18),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("Average RL") +
  geom_signif(comparisons=list(c("conserved","mutating")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ylim(0,1.05) +
  ggtitle("Ablang1")

#Plot Figure 4H
pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figure4/ConservedMutating_RL.pdf", width = 14, height = 5.5)
ggarrange(p, esm1b, o, ab1, nrow = 1, ncol = 4)
dev.off()

#Figure 4F and G
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

source("~/OneDrive - UMC Utrecht/Documenten/Steropodon_theme.R")

ablang_prob_matrix <- read.csv("~/Documents/GitHub/PLM-likelihoods/data/OVA_V7/edges_prob_matrix/prob_matrix_edge0_ablang2_OVA_V7_default.csv", header = T)
esm_prob_matrix <- read.csv("~/Documents/GitHub/PLM-likelihoods/data/OVA_V7/edges_prob_matrix/prob_matrix_edge0_esmc_OVA_V7_default.csv", header = T)

df_esm <- esm_prob_matrix[100,]
df_esm <- df_esm[,5:24]
df_esm <- pivot_longer(df_esm, cols = colnames(df_esm), names_to = "AA", values_to = "Probability")

df_ablang <- ablang_prob_matrix[100,]
df_ablang <- pivot_longer(df_ablang, cols = colnames(df_ablang), names_to = "AA", values_to = "Probability")

a <- ggplot(df_ablang, aes(x = AA, y = Probability, fill = AA)) +
  geom_bar(stat = "identity") + 
  ylim(0,1) +
  theme_steropodon() +
  theme(text = element_text(size = 20),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ggtitle("Ablang2")

b <- ggplot(df_esm, aes(x = AA, y = Probability, fill = AA)) +
  geom_bar(stat = "identity") + 
  ylim(0,1) +
  theme_steropodon() +
  theme(text = element_text(size = 20),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ggtitle("ESM-C")

pdf("~/Documents/GitHub/PLM-likelihoods/figures/figureS11/Example_probabilities.pdf", width = 12, height = 6)
a <- ggarrange(a, b, ncol = 2)
dev.off()

#evenness index
df <- data.frame("PLM" = character(), "SEI" = numeric())
for (data in c("OVA_V7", "Bruhn", "Kim", "horns2020a__VDJ_RAW")){
  print(data)
  prob_matrix_files = list.files(paste0("~/Documents/GitHub/PLM-likelihoods/data/", data, "/edges_prob_matrix/"))
  for (model in c("ablang1", "ablang2", "sapiens", "esm1b", "esmc", "protbert")){
    file_list <- parallel::mclapply(mc.cores = num.cores, as.list(prob_matrix_files[grep(model, prob_matrix_files)]), function(file){
      prob_matrix <- read.csv(paste0("~/Documents/GitHub/PLM-likelihoods/data/", data, "/edges_prob_matrix/", file), header = T)
      if(model == "esmc"){prob_matrix <- prob_matrix[,5:24]}
      if(model == "probert"){prob_matrix <- prob_matrix[,6:25]}
      num.cores <- parallel::detectCores() -1
      prob_list <- as.list(data.frame(t(prob_matrix)))
      SEI = parallel::mclapply(mc.cores = num.cores, prob_list, function(x) {
        probs <- as.numeric(x)
        SEI <- vegan::diversity(x = probs, index = "shannon") / log(length(probs))
      })
      df_file <- data.frame("PLM" = rep(model, length(SEI)), "SEI" = unlist(SEI))
      return(df_file)
    })
    df_model <- do.call(rbind, file_list)
    df <- rbind(df, df_model)
  }
  #add esmc 600
  model <- "esmc_600"
  prob_matrix_files = list.files(paste0("~/Documents/GitHub/PLM-likelihoods/data/", data, "/edges_prob_matrix_esmc600/"))
  file_list <- parallel::mclapply(mc.cores = num.cores, as.list(prob_matrix_files), function(file){
    prob_matrix <- read.csv(paste0("~/Documents/GitHub/PLM-likelihoods/data/", data, "/edges_prob_matrix_esmc600/", file), header = T)
    prob_matrix <- prob_matrix[,5:24]
    num.cores <- parallel::detectCores() -1
    prob_list <- as.list(data.frame(t(prob_matrix)))
    SEI = parallel::mclapply(mc.cores = num.cores, prob_list, function(x) {
      probs <- as.numeric(x)
      SEI <- vegan::diversity(x = probs, index = "shannon") / log(length(probs))
    })
    df_file <- data.frame("PLM" = rep(model, length(SEI)), "SEI" = unlist(SEI))
    return(df_file)
  })
  df_model <- do.call(rbind, file_list)
  df <- rbind(df, df_model)
}
save(df, file = "~/Documents/GitHub/PLM-likelihoods/figures/figure4/SEI_df_v2.RData")

pdf("~/Documents/GitHub/PLM-likelihoods/figures/figure4/SEI_boxplot.pdf")
b <- ggplot(df, aes(x = PLM, y = SEI)) +
  geom_boxplot(outliers = F, fill= "lightblue") +
  theme_steropodon() +
  theme(text = element_text(size = 20),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Shannon Evenness Index")
dev.off()

pdf("~/Documents/GitHub/PLM-likelihoods/figures/figure4/FigureS11.pdf", width = 20, height = 6)
ggarrange(a, b, ncol = 2, nrow = 1)
dev.off()

