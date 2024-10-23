##Figure 5 - Phylogenetic analysis of clonal lineages
library(igraph)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)

#Figure 5A+B
#Read data
load("PLM_Likelihoods/data/OVA_V7/AF_OVA_V7_default_HC.RData") #af
af_ova <- af
rm(af)
load("PLM_Likelihoods/data/Horns/AF_horns2020a__VDJ_RAW_default_HC.RData")
af_horns <- af
rm(af)
load("PLM_Likelihoods/data/Bruhn/AF_Bruhn_default_HC.RData")
af_bruhn <- af
rm(af)
load("PLM_Likelihoods/data/Kim/AF_Kim_default_HC.RData")
af_kim <- af
rm(af)

#A - Lineage Trees
#Source AntibodyForests function from https://github.com/alexyermanos/AntibodyForests
source("~/Documents/GitHub/Platypus/R/AntibodyForests_plot.R")

#Plot Figure 5A
pdf("PLM_Likelihoods/figures/Figure5/LineageTree_OVA.pdf")
AntibodyForests_plot(af_ova,
                     sample = "S1",
                     clonotype = "clonotype4",
                     color.by = "IGH_evo_likelihood_esm_full_VDJ",
                     label.by = "size",
                     show.color.legend = F,
                     show.size.legend = F,
                     main.title = "Mouse 1",
                     sub.title = "Clonotype 4",
                     edge.label = "original",
                     arrow.size = 0.5)
dev.off()

#Plot Figure 5B
pdf("PLM_Likelihoods/figures/Figure5_Supplementary7/LineageTree_Kim.pdf")
AntibodyForests_plot(af_kim,
                     sample = "SRR17729703",
                     clonotype = "clonotype1",
                     color.by = "IGH_evo_likelihood_esm_full_VDJ",
                     label.by = "size",
                     show.color.legend = F,
                     show.size.legend = F,
                     main.title = "Individual 3",
                     sub.title = "Clonotype 1",
                     edge.label = "original")
dev.off()

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
      esm = unlist(lapply(node_features,function(x){unique(x[["IGH_evo_likelihood_esm_full_VDJ"]])[!is.na(unique(x[["IGH_evo_likelihood_esm_full_VDJ"]]))]}))[rownames(distance)]
      protbert = unlist(lapply(node_features,function(x){unique(x[["IGH_evo_likelihood_protbert_full_VDJ"]])[!is.na(unique(x[["IGH_evo_likelihood_protbert_full_VDJ"]]))]}))[rownames(distance)]
      ablang = unlist(lapply(node_features,function(x){unique(x[["IGH_evo_likelihood_ablang_full_VDJ"]])[!is.na(unique(x[["IGH_evo_likelihood_ablang_full_VDJ"]]))]}))[rownames(distance)]
      sapiens = unlist(lapply(node_features,function(x){unique(x[["IGH_evo_likelihood_sapiens_full_VDJ"]])[!is.na(unique(x[["IGH_evo_likelihood_sapiens_full_VDJ"]]))]}))[rownames(distance)]
      
      distance = data.frame("distance" = distance, "sample" = sample_name, "clonotype" = clonotype_name, 
                            "node" = rownames(distance), "esm" = esm, "protbert" = protbert, 
                            "ablang" = ablang, "sapiens" = sapiens)
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
#esm
cor_human <- cor.test(distance_human$germline, distance_human$esm)$estimate
cor_mouse <- cor.test(distance_ova$germline, distance_ova$esm)$estimate
p <- ggplot(distance_all, aes(x = germline, y = esm, colour = sample)) + 
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
  theme_minimal() +
  theme(text = element_text(size = 18)) +
  xlab("Distance to germline") + ylab("ESM-1b Pseudolikelihood") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#Plot Figure 4C
pdf("PLM_Likelihoods/figures/Figure4_Supplementary67/DistanceGermline_main.pdf")
print(p)
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
  theme_minimal() +
  theme(text = element_text(size = 18)) +
  xlab("Distance to germline") + ylab("ProtBERT Pseudolikelihood") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mouse, digits = 2)))
#ablang
cor_human <- cor.test(distance_human$germline, distance_human$ablang)$estimate
cor_mouse <- cor.test(distance_ova$germline, distance_ova$ablang)$estimate
b <- ggplot(distance_all, aes(x = germline, y = ablang, colour = sample)) + 
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
  theme_minimal() +
  theme(text = element_text(size = 18)) +
  xlab("Distance to germline") + ylab("Ablang Pseudolikelihood") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mouse, digits = 2)))
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
  theme_minimal() +
  theme(text = element_text(size = 18)) +
  xlab("Distance to germline") + ylab("Sapiens Pseudolikelihood") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#Plot supplementary Figure 6
pdf("PLM_Likelihoods/figures/Figure4_Supplementary67/DistanceGermline_sup.pdf", width = 8, height = 6)
ggarrange(a, b, c, ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()

#Figure 4DE
#Read data
mut_ova <- read.csv("PLM_Likelihoods/data/OVA_V7/mutational_rank_table_OVA_V7_default.csv")
mut_horns <- read.csv("PLM_Likelihoods/data/Horns/mutational_rank_table_horns2020a__VDJ_RAW_default.csv")
mut_bruhn <- read.csv("PLM_Likelihoods/data/Bruhn/mutational_rank_table_Bruhn_default.csv")
mut_kim <- read.csv("PLM_Likelihoods/data/Kim/mutational_rank_table_Kim_default.csv")
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
mut_all_esm <- mut_all[mut_all$model == "esm",]
mut_all_ablang <- mut_all[mut_all$model == "ablang",]
mut_all_protbert <- mut_all[mut_all$model == "protbert",]
mut_all_sapiens <- mut_all[mut_all$model == "sapiens",]

# C - Distribution of likelihood rank of the substitution along the edge
a <- ggplot(mut_all_esm, aes(mean_sub_rank)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 1, linewidth = 2) +
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
  theme_minimal() +
  ylim(0,250)+
  scale_x_continuous(breaks = c(1,5,10,15,20), limits = c(1,NA)) +
  theme(text = element_text(size = 16)) +
  xlab("Substitution Rank") + ylab("Number of edges") +
  ggtitle("ESM-1b")
b <- ggplot(mut_all_ablang, aes(mean_sub_rank)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 1, linewidth = 2) +
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
  theme_minimal() +
  ylim(0,250)+
  scale_x_continuous(breaks = c(1,5,10,15,20), limits = c(1,NA)) +
  theme(text = element_text(size = 16)) +
  xlab("Substitution Rank") + ylab("Number of edges") +
  ggtitle("Ablang")

#Plot Figure 5C
pdf("PLM_Likelihoods/figures/Figure4_Supplementary67/Figure4D.pdf", width = 10, height = 4.5)
print(ggarrange(a, b, ncol = 2, common.legend = T, legend = "right"))
dev.off()

#E - quantification of mean likelihood rank per sample for each model
mut_all <- mut_all[mut_all$model != "",]
mut_all_esm %>% group_by(sample) %>% mutate(average = mean(as.numeric(mean_sub_rank))) -> mut_all_esm
mut_all_ablang %>% group_by(sample) %>% mutate(average = mean(as.numeric(mean_sub_rank))) -> mut_all_ablang
mut_all_sapiens %>% group_by(sample) %>% mutate(average = mean(as.numeric(mean_sub_rank))) -> mut_all_sapiens
mut_all_protbert %>% group_by(sample) %>% mutate(average = mean(as.numeric(mean_sub_rank))) -> mut_all_protbert
mut_all <- rbind(mut_all_esm, mut_all_ablang, mut_all_sapiens, mut_all_protbert)
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
  theme_minimal() +
  theme(text = element_text(size = 16)) +
  ylab("Likelihood Rank") + xlab("Model") +
  scale_x_discrete(labels = c("esm" = "ESM-1b", "ablang" = "Ablang", "protbert" = "ProtBERT", "sapiens" = "Sapiens")) +
  ggtitle("Average Substitution Rank")

#Plot Figure 5D
pdf("PLM_Likelihoods/figures/Figure4_Supplementary67/Figure4E.pdf", width = 7, height = 5.5)
print(p)
dev.off()

#Figure E 
#Read data convserved residues
conserved_ova <- read.csv("PLM_Likelihoods/data/OVA_V7/conserved_rank_table_OVA_V7_default.csv")
conserved_horns <- read.csv("PLM_Likelihoods/data/Horns/conserved_rank_table_horns2020a__VDJ_RAW_default.csv")
conserved_bruhn <- read.csv("PLM_Likelihoods/data/Bruhn/conserved_rank_table_Bruhn_default.csv")
conserved_kim <- read.csv("PLM_Likelihoods/data/Kim/conserved_rank_table_Kim_default.csv")

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
conserved_all_esm <- conserved_all[conserved_all$model == "esm",]
conserved_all_protbert <- conserved_all[conserved_all$model == "protbert",]
conserved_all_sapiens <- conserved_all[conserved_all$model == "sapiens",]
conserved_all_ablang <- conserved_all[conserved_all$model == "ablang",]

#read data original mutated residues
original_ova <- read.csv("PLM_Likelihoods/data/OVA_V7/mutational_rank_reversed_table_OVA_V7_default.csv")
original_horns <- read.csv("PLM_Likelihoods/data/Horns/mutational_rank_reversed_table_horns2020a__VDJ_RAW_default.csv")
original_bruhn <- read.csv("PLM_Likelihoods/data/Bruhn/mutational_rank_reversed_table_Bruhn_default.csv")
original_kim <- read.csv("PLM_Likelihoods/data/Kim/mutational_rank_reversed_table_Kim_default.csv")
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
original_all_esm <- original_all[original_all$model == "esm",]
original_all_protbert <- original_all[original_all$model == "protbert",]
original_all_sapiens <- original_all[original_all$model == "sapiens",]
original_all_ablang <- original_all[original_all$model == "ablang",]

#Figure E - Significance boxplot of likelihood ranks
#esm
conserved_all_esm$group <- "conserved"
original_all_esm$group <- "mutating"
df <- rbind(conserved_all_esm[,-8], original_all_esm[-8])
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
  ylab("ESM-1b Likelihood rank") +
  geom_signif(comparisons=list(c("conserved","mutating")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ylim(0,22) +
  ggtitle("ESM-1b")

#ablang
conserved_all_ablang$group <- "conserved"
original_all_ablang$group <- "mutating"
df <- rbind(conserved_all_ablang[,-8], original_all_ablang[-8])
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
  ylab("Ablang Likelihood rank") +
  geom_signif(comparisons=list(c("conserved","mutating")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ylim(0,22) +
  ggtitle("Ablang")

#Plot Figure 4F
pdf("PLM_Likelihoods/figures/Figure4_Supplementary7/Figure4F.pdf", width = 7, height = 5.5)
ggarrange(p, o, nrow = 1, ncol = 2)
dev.off()

#Figure G - distribution of likelihood of the mutating and conserved residues along the edges
#esm original
a <- ggplot(original_all_esm, aes(mean_sub_prob)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 0.1, linewidth = 2) +
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
  theme_minimal() +
  theme(text = element_text(size = 12)) +
  xlim(0,1) + ylim(0,800) +
  xlab("Original Residue Likelihood") + ylab("Number of edges") +
  ggtitle("Original - ESM-1b")

#esm conserved
b <- ggplot(conserved_all_esm, aes(mean_sub_prob)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 0.1, linewidth = 2) +
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
  theme_minimal() +
  theme(text = element_text(size = 12)) +
  xlim(0,1) + ylim(0,800)+
  xlab("Conserved Residue Likelihood") + ylab("Number of edges") +
  ggtitle("Conserved - ESM-1b")

#ablang original
c <- ggplot(original_all_ablang, aes(mean_sub_prob)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 0.1, linewidth = 2) +
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
  theme_minimal() +
  theme(text = element_text(size = 12)) +
  xlim(0,1) + ylim(0,800) +
  xlab("Original Residue Likelihood") + ylab("Number of edges") +
  ggtitle("Original - Ablang")

#ablang conserved
d <- ggplot(conserved_all_ablang, aes(mean_sub_prob)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 0.1, linewidth = 2) +
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
  theme_minimal() +
  theme(text = element_text(size = 12)) +
  xlim(0,1) + ylim(0,800)+
  xlab("Conserved Residue Likelihood") + ylab("Number of edges") +
  ggtitle("Conserved - Ablang")

#Plot Figure 4G
pdf("PLM_Likelihoods/figures/Figure4_Supplementary67/Figure4G.pdf", width = 12, height = 4.5)
ggarrange(a, b, c, d, ncol = 4, common.legend = T, legend = "right")
dev.off()

#Figure H - Significance boxplot of likelihood of the mutating and conserved residues
conserved_all_esm$group <- "conserved"
original_all_esm$group <- "mutating"
df <- rbind(conserved_all_esm[,-8], original_all_esm[-8])
stat.test <- df %>%
  t_test(mean_sub_prob ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test
p <- ggplot(df, aes(y = as.numeric(mean_sub_prob), x = group, color=group)) +
  geom_boxplot() +
  scale_color_manual(values = c("conserved" = "orange",
                                "mutating" = "dodgerblue3")) +
  theme_minimal() +
  theme(text = element_text(size = 18),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("ESM-1b Likelihood") +
  geom_signif(comparisons=list(c("conserved","mutating")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ylim(0,1.05) +
  ggtitle("ESM-1b")

conserved_all_ablang$group <- "conserved"
original_all_ablang$group <- "mutating"
df <- rbind(conserved_all_ablang[,-8], original_all_ablang[-8])
stat.test <- df %>%
  t_test(mean_sub_prob ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test
o <- ggplot(df, aes(y = as.numeric(mean_sub_prob), x = group, color=group)) +
  geom_boxplot() +
  scale_color_manual(values = c("conserved" = "orange",
                                "mutating" = "dodgerblue3")) +
  theme_minimal() +
  theme(text = element_text(size = 18),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("Ablang Likelihood") +
  geom_signif(comparisons=list(c("conserved","mutating")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ylim(0,1.05) +
  ggtitle("Ablang")

#Plot Figure 4H
pdf("PLM_Likelihoods/figures/Figure4_Supplementary67/Figure4H.pdf", width = 7, height = 5.5)
ggarrange(p, o, nrow = 1, ncol = 2)
dev.off()



