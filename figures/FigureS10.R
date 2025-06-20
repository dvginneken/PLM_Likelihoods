##Phylogenetic analysis of clonal lineages of extra samples
library(igraph)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(AntibodyForests)

source("~/OneDrive - UMC Utrecht/Documenten/Steropodon_theme.R")

#Figure 5A+B
#Read data
load("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Mathew/AF_Mathew_default_HC.RData") #af
af_mouse <- af
rm(af)
load("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim_extra/AF_Kim_extra_default_HC.RData")
af_human <- af
rm(af)

#Read VDJ dataframes
load("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Mathew/VDJ_PLL_Mathew.RData")
vdj_mouse <- vdj_likelihood
rm(vdj_liukelihood)
load("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim_extra/VDJ_PLL_Kim_extra.RData")
vdj_human <- vdj_likelihood
rm(vdj_likelihood)

#Add pseudolikelihoods
af_mouse <- Af_add_node_feature(af_mouse, vdj_mouse, c("HC_evo_likelihood_esmc_full_VDJ",
                                                       "HC_evo_likelihood_protbert_full_VDJ",
                                                       "HC_evo_likelihood_ablang1_full_VDJ",
                                                       "HC_evo_likelihood_sapiens_full_VDJ",
                                                       "HC_evo_likelihood_esm1b_full_VDJ",
                                                 "paired_evo_likelihood_ablang2_full_VDJ"))
af_human <- Af_add_node_feature(af_human, vdj_human, c("HC_evo_likelihood_esmc_full_VDJ",
                                                       "HC_evo_likelihood_protbert_full_VDJ",
                                                       "HC_evo_likelihood_ablang1_full_VDJ",
                                                       "HC_evo_likelihood_sapiens_full_VDJ",
                                                       "HC_evo_likelihood_esm1b_full_VDJ",
                                                       "paired_evo_likelihood_ablang2_full_VDJ"))
#A - Lineage Trees
#Plot Figure 5A
Af_plot_tree(af_mouse,
             sample = "S15.D28.1.mln",
             clonotype = "clonotype6",
             color.by = "HC_evo_likelihood_esmc_full_VDJ",
             label.by = "size",
             show.size.legend = F,
             main.title = "Mouse sample: S15.D28.1.mln",
             sub.title = "Clonotype 6",
             arrow.size = 0.5,
             output.file = "/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/FigureS9/MouseTree_esmc.pdf")



#Plot Figure 5B
Af_plot_tree(af_human,
             sample = "SRR17729681",
             clonotype = "clonotype7",
             color.by = "HC_evo_likelihood_esmc_full_VDJ",
             label.by = "size",
             show.size.legend = F,
             main.title = "Human sample: SRR17729681",
             sub.title = "Clonotype 7",
             arrow.size = 0.5,
             output.file = "/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/FigureS9/HumanTree_esmc.pdf")

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
      esm1b = unlist(lapply(node_features,function(x){unique(x[["HC_evo_likelihood_esm1b_full_VDJ"]])[!is.na(unique(x[["HC_evo_likelihood_esm1b_full_VDJ"]]))]}))[rownames(distance)]
      protbert = unlist(lapply(node_features,function(x){unique(x[["HC_evo_likelihood_protbert_full_VDJ"]])[!is.na(unique(x[["HC_evo_likelihood_protbert_full_VDJ"]]))]}))[rownames(distance)]
      ablang1 = unlist(lapply(node_features,function(x){unique(x[["HC_evo_likelihood_ablang1_full_VDJ"]])[!is.na(unique(x[["HC_evo_likelihood_ablang1_full_VDJ"]]))]}))[rownames(distance)]
      sapiens = unlist(lapply(node_features,function(x){unique(x[["HC_evo_likelihood_sapiens_full_VDJ"]])[!is.na(unique(x[["HC_evo_likelihood_sapiens_full_VDJ"]]))]}))[rownames(distance)]
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
distance_mouse <- distance_germline(af_mouse)
distance_human <- distance_germline(af_human)
distance_all <- rbind(distance_mouse, distance_human)

#Plot the distance to the germline against the pseudolikelihoods
#esmc
cor_human <- cor.test(distance_human$germline, distance_human$esmc)$estimate
cor_mouse <- cor.test(distance_mouse$germline, distance_mouse$esmc)$estimate
p <- ggplot(distance_all, aes(x = germline, y = esmc, colour = sample)) + 
  geom_point(size = 1) + 
  geom_smooth(method = "lm", se=F) +
  scale_color_manual(values = c("SRR17729681" = "#99d8c9",
                                "SRR17729674" = "#66c2a4",
                                "SRR17729725" = "#41ae76",
                                "SRR17729673" = "#238b45",
                                "SRR17729714" = "#005824",
                                "S11.D14.3.spleen" = "#fcc5c0",
                                "S14.D28.1.spleen" = "#fa9fb5",
                                "S15.D28.1.mln" = "#f768a1",
                                "S17.D28.2.spleen" = "#c51b8a",
                                "S18.D28.2.mln" = "#7a0177",
                                "S20.D28.3.mln" = "#49006a"),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 18)) +
  xlab("Distance to germline") + ylab("ESM-C SP") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), "\nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#esm1b
cor_human <- cor.test(distance_human$germline, distance_human$esm1b)$estimate
cor_mouse <- cor.test(distance_mouse$germline, distance_mouse$esm1b)$estimate
p2 <- ggplot(distance_all, aes(x = germline, y = esm1b, colour = sample)) + 
  geom_point(size = 1) + 
  geom_smooth(method = "lm", se=F) +
  scale_color_manual(values = c("SRR17729681" = "#99d8c9",
                                "SRR17729674" = "#66c2a4",
                                "SRR17729725" = "#41ae76",
                                "SRR17729673" = "#238b45",
                                "SRR17729714" = "#005824",
                                "S11.D14.3.spleen" = "#fcc5c0",
                                "S14.D28.1.spleen" = "#fa9fb5",
                                "S15.D28.1.mln" = "#f768a1",
                                "S17.D28.2.spleen" = "#c51b8a",
                                "S18.D28.2.mln" = "#7a0177",
                                "S20.D28.3.mln" = "#49006a"),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 18)) +
  xlab("Distance to germline") + ylab("ESM-1b SP") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), "\nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#ablang2
cor_human <- cor.test(distance_human$germline, distance_human$ablang2)$estimate
cor_mouse <- cor.test(distance_mouse$germline, distance_mouse$ablang2)$estimate
q <- ggplot(distance_all, aes(x = germline, y = ablang2, colour = sample)) + 
  geom_point(size = 1) + 
  geom_smooth(method = "lm", se=F) +
  scale_color_manual(values = c("SRR17729681" = "#99d8c9",
                                "SRR17729674" = "#66c2a4",
                                "SRR17729725" = "#41ae76",
                                "SRR17729673" = "#238b45",
                                "SRR17729714" = "#005824",
                                "S11.D14.3.spleen" = "#fcc5c0",
                                "S14.D28.1.spleen" = "#fa9fb5",
                                "S15.D28.1.mln" = "#f768a1",
                                "S17.D28.2.spleen" = "#c51b8a",
                                "S18.D28.2.mln" = "#7a0177",
                                "S20.D28.3.mln" = "#49006a"),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 18)) +
  xlab("Distance to germline") + ylab("Ablang2 SP") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), "\nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#ablang1
cor_human <- cor.test(distance_human$germline, distance_human$ablang1)$estimate
cor_mouse <- cor.test(distance_mouse$germline, distance_mouse$ablang1)$estimate
q2 <- ggplot(distance_all, aes(x = germline, y = ablang1, colour = sample)) + 
  geom_point(size = 1) + 
  geom_smooth(method = "lm", se=F) +
  scale_color_manual(values = c("SRR17729681" = "#99d8c9",
                                "SRR17729674" = "#66c2a4",
                                "SRR17729725" = "#41ae76",
                                "SRR17729673" = "#238b45",
                                "SRR17729714" = "#005824",
                                "S11.D14.3.spleen" = "#fcc5c0",
                                "S14.D28.1.spleen" = "#fa9fb5",
                                "S15.D28.1.mln" = "#f768a1",
                                "S17.D28.2.spleen" = "#c51b8a",
                                "S18.D28.2.mln" = "#7a0177",
                                "S20.D28.3.mln" = "#49006a"),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 18)) +
  xlab("Distance to germline") + ylab("Ablang1 SP") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), "\nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))


pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figureS10/Distance.pdf", width = 20, height = 5)
ggarrange(p, p2, q, q2, ncol = 4, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()

#Read data
mut_mouse <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Mathew/mutational_rank_table_Mathew_default.csv")
mut_mouse2 <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Mathew/mutational_rank_table_Mathew_default_2.csv")
mut_mouse <- rbind(mut_mouse, mut_mouse2)
mut_human <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim_extra/mutational_rank_table_Kim_extra_default.csv")
mut_human2 <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim_extra/mutational_rank_table_Kim_extra_default_2.csv")
mut_human <- rbind(mut_human, mut_human2)
mut_human <- mut_human[mut_human$sample %in% c("SRR17729681",
                                              "SRR17729674",
                                              "SRR17729725",
                                              "SRR17729673",
                                              "SRR17729714"),]

mut_all <- rbind(mut_mouse, mut_human)

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
  scale_color_manual(values = c("SRR17729681" = "#99d8c9",
                                "SRR17729674" = "#66c2a4",
                                "SRR17729725" = "#41ae76",
                                "SRR17729673" = "#238b45",
                                "SRR17729714" = "#005824",
                                "S11.D14.3.spleen" = "#fcc5c0",
                                "S14.D28.1.spleen" = "#fa9fb5",
                                "S15.D28.1.mln" = "#f768a1",
                                "S17.D28.2.spleen" = "#c51b8a",
                                "S18.D28.2.mln" = "#7a0177",
                                "S20.D28.3.mln" = "#49006a"),
                     name = "Sample") +
  theme_steropodon() +
  ylim(0,400)+
  scale_x_continuous(breaks = c(1,5,10,15,20), limits = c(1,NA)) +
  theme(text = element_text(size = 16)) +
  xlab("Substitution Rank") + ylab("Number of edges") +
  ggtitle("ESM-C")

b <- ggplot(mut_all_esm1b, aes(mean_sub_rank)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 1, linewidth = 1) +
  scale_color_manual(values = c("SRR17729681" = "#99d8c9",
                                "SRR17729674" = "#66c2a4",
                                "SRR17729725" = "#41ae76",
                                "SRR17729673" = "#238b45",
                                "SRR17729714" = "#005824",
                                "S11.D14.3.spleen" = "#fcc5c0",
                                "S14.D28.1.spleen" = "#fa9fb5",
                                "S15.D28.1.mln" = "#f768a1",
                                "S17.D28.2.spleen" = "#c51b8a",
                                "S18.D28.2.mln" = "#7a0177",
                                "S20.D28.3.mln" = "#49006a"),
                     name = "Sample") +
  theme_steropodon() +
  ylim(0,400)+
  scale_x_continuous(breaks = c(1,5,10,15,20), limits = c(1,NA)) +
  theme(text = element_text(size = 16)) +
  xlab("Substitution Rank") + ylab("Number of edges") +
  ggtitle("ESM-1b")

c <- ggplot(mut_all_ablang2, aes(mean_sub_rank)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 1, linewidth = 1) +
  scale_color_manual(values = c("SRR17729681" = "#99d8c9",
                                "SRR17729674" = "#66c2a4",
                                "SRR17729725" = "#41ae76",
                                "SRR17729673" = "#238b45",
                                "SRR17729714" = "#005824",
                                "S11.D14.3.spleen" = "#fcc5c0",
                                "S14.D28.1.spleen" = "#fa9fb5",
                                "S15.D28.1.mln" = "#f768a1",
                                "S17.D28.2.spleen" = "#c51b8a",
                                "S18.D28.2.mln" = "#7a0177",
                                "S20.D28.3.mln" = "#49006a"),
                     name = "Sample") +
  theme_steropodon() +
  ylim(0,400)+
  scale_x_continuous(breaks = c(1,5,10,15,20), limits = c(1,NA)) +
  theme(text = element_text(size = 16)) +
  xlab("Substitution Rank") + ylab("Number of edges") +
  ggtitle("Ablang2")

d <- ggplot(mut_all_ablang1, aes(mean_sub_rank)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 1, linewidth = 1) +
  scale_color_manual(values = c("SRR17729681" = "#99d8c9",
                                "SRR17729674" = "#66c2a4",
                                "SRR17729725" = "#41ae76",
                                "SRR17729673" = "#238b45",
                                "SRR17729714" = "#005824",
                                "S11.D14.3.spleen" = "#fcc5c0",
                                "S14.D28.1.spleen" = "#fa9fb5",
                                "S15.D28.1.mln" = "#f768a1",
                                "S17.D28.2.spleen" = "#c51b8a",
                                "S18.D28.2.mln" = "#7a0177",
                                "S20.D28.3.mln" = "#49006a"),
                     name = "Sample") +
  theme_steropodon() +
  ylim(0,400)+
  scale_x_continuous(breaks = c(1,5,10,15,20), limits = c(1,NA)) +
  theme(text = element_text(size = 16)) +
  xlab("Substitution Rank") + ylab("Number of edges") +
  ggtitle("Ablang1")

#Plot Figure 5C
pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figureS10/MutationRank.pdf", width = 20, height = 5)
print(ggarrange(a, b, c, d, ncol = 4, common.legend = T, legend = "right"))
dev.off()

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
  scale_color_manual(values = c("SRR17729681" = "#99d8c9",
                                "SRR17729674" = "#66c2a4",
                                "SRR17729725" = "#41ae76",
                                "SRR17729673" = "#238b45",
                                "SRR17729714" = "#005824",
                                "S11.D14.3.spleen" = "#fcc5c0",
                                "S14.D28.1.spleen" = "#fa9fb5",
                                "S15.D28.1.mln" = "#f768a1",
                                "S17.D28.2.spleen" = "#c51b8a",
                                "S18.D28.2.mln" = "#7a0177",
                                "S20.D28.3.mln" = "#49006a"),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 16)) +
  ylab("Likelihood Rank") + xlab("Model") +
  scale_x_discrete(labels = c("esm" = "ESM-1b", "ablang" = "Ablang", "protbert" = "ProtBERT", "sapiens" = "Sapiens")) +
  ggtitle("Average Substitution Rank")

#Plot Figure 5D
pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figureS9/AverageSubstitutionRank.pdf", width = 7, height = 5.5)
print(p)
dev.off()

conserved_kim_extra <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim_extra/conserved_rank_table_Kim_extra_default.csv")
conserved_kim_extra2 <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim_extra/conserved_rank_table_Kim_extra_default_2.csv")
conserved_kim_extra <- rbind(conserved_kim_extra, conserved_kim_extra2)
conserved_mathew <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Mathew/conserved_rank_table_Mathew_default.csv")
conserved_mathew2 <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Mathew/conserved_rank_table_Mathew_default_2.csv")
conserved_mathew <- rbind(conserved_mathew, conserved_mathew2)

conserved_all <- rbind(conserved_kim_extra, conserved_mathew)

#separate on model
conserved_all_esmc <- conserved_all[conserved_all$model == "esmc",]
conserved_all_protbert <- conserved_all[conserved_all$model == "protbert",]
conserved_all_sapiens <- conserved_all[conserved_all$model == "sapiens",]
conserved_all_ablang2 <- conserved_all[conserved_all$model == "ablang2",]
conserved_all_ablang1 <- conserved_all[conserved_all$model == "ablang1",]
conserved_all_esm1b <- conserved_all[conserved_all$model == "esm1b",]


original_kim_extra <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim_extra/mutational_rank_reversed_table_Kim_extra_default.csv")
original_kim_extra2 <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim_extra/mutational_rank_reversed_table_Kim_extra_default_2.csv")
original_kim_extra <- rbind(original_kim_extra, original_kim_extra2)
original_mathew <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Mathew/mutational_rank_reversed_table_Mathew_default.csv")
original_mathew2 <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Mathew/mutational_rank_reversed_table_Mathew_default_2.csv")
original_mathew <- rbind(original_mathew, original_mathew2)

original_all <- rbind(original_kim_extra, original_mathew)

#separate on model
original_all_esmc <- original_all[original_all$model == "esmc",]
original_all_protbert <- original_all[original_all$model == "protbert",]
original_all_sapiens <- original_all[original_all$model == "sapiens",]
original_all_ablang2 <- original_all[original_all$model == "ablang2",]
original_all_ablang1 <- original_all[original_all$model == "ablang1",]
original_all_esm1b <- original_all[original_all$model == "esm1b",]


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

#esm1b
conserved_all_esm1b$group <- "conserved"
original_all_esm1b$group <- "mutating"
df <- rbind(conserved_all_esm1b[,-8], original_all_esm1b[-8])
stat.test <- df %>%
  t_test(mean_sub_prob ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test
p2 <- ggplot(df, aes(y = as.numeric(mean_sub_prob), x = group, color=group)) +
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

conserved_all_ablang1$group <- "conserved"
original_all_ablang1$group <- "mutating"
df <- rbind(conserved_all_ablang1[,-8], original_all_ablang1[-8])
stat.test <- df %>%
  t_test(mean_sub_prob ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test
o2 <- ggplot(df, aes(y = as.numeric(mean_sub_prob), x = group, color=group)) +
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
pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figureS10/ConservedMutating_RL.pdf", width = 14, height = 5.5)
ggarrange(p, p2, o, o2, nrow = 1, ncol = 4)
dev.off()

#SEI
num.cores <- parallel::detectCores() - 1
df <- data.frame("PLM" = character(), "SEI" = numeric())
for (data in c("Kim_extra", "Mathew")){
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
save(df, file = "/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figureS10/SEI_df.RData")

pdf("~/Documents/GitHub/PLM-likelihoods/figures/figureS10/SEI_boxplot.pdf", width = 14, height = 5)
p <- ggplot(df, aes(x = PLM, y = SEI)) +
  geom_boxplot(outliers = F, fill= "lightblue") +
  theme_steropodon() +
  theme(text = element_text(size = 20),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Shannon Evenness Index")
dev.off()
