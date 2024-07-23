library(ggplot2) 
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(ggpubr)

#Figure 3A
#Read data
df_ova <- read.csv("PLM_Likelihoods/data/OVA_V7/PLMCorrelation.csv",header = TRUE, sep = ",")
df_horns <- read.csv("PLM_Likelihoods/data/Horns/PLMCorrelation.csv",header = TRUE, sep = ",")
df_bruhn <- read.csv("PLM_Likelihoods/data/Bruhn/PLMCorrelation.csv",header = TRUE, sep = ",")
df_kim <- read.csv("PLM_Likelihoods/data/Kim/PLMCorrelation.csv",header = TRUE, sep = ",")

#change sample names
df_ova$sample <- case_match(df_ova$sample,
                            "S1" ~ "Mouse1",
                            "S2" ~ "Mouse2",
                            "S3" ~ "Mouse3",
                            "S4" ~ "Mouse4",
                            "S5" ~ "Mouse5")
df_horns$sample <- case_match(df_horns$sample,
                              "Influenza.vac.11.12.human.S1" ~ "Individual1",
                              "Influenza.vac.11.12.human.S2" ~ "Human2",
                              "Influenza.vac.11.12.human.S3" ~ "Human3",
                              "Influenza.vac.11.12.human.S4" ~ "Human4")
#Only keep 1 replicate
df_horns <- df_horns[df_horns$sample == "Individual1",]
df_bruhn$sample <- "Individual2"
df_kim <- df_kim[df_kim$sample %in% c("SRR17729703", "SRR17729692", "SRR17729726"),]
df_kim$sample <- case_match(df_kim$sample,
                            "SRR17729703" ~ "Individual3",
                            "SRR17729692" ~ "Individual4",
                            "SRR17729726" ~ "Individual5")
df <- rbind(df_ova, df_horns, df_bruhn, df_kim)
df <- pivot_longer(df, cols = 1:6, names_to = "PLM", values_to = "correlation")

#Plot figure 3A
pdf("PLM_Likelihoods/figures/Figure3_Supplementary3/PLMcorrelation_main.pdf", width = 8, height = 6)
ggplot(df, aes(x=PLM, y=correlation, col=factor(sample))) + 
  geom_point(size = 4) +
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
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  xlab("PLM Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~source)
dev.off()

#Supplementary Figure 3
#Read data
df_ova <- read.csv("PLM_Likelihoods/data/OVA_V7/PLMCorrelation_chains.csv",header = TRUE, sep = ",")
df_horns <- read.csv("PLM_Likelihoods/data/Horns/PLMCorrelation_chains.csv",header = TRUE, sep = ",")
df_bruhn <- read.csv("PLM_Likelihoods/PLM-likelihoods/data/Bruhn/PLMCorrelation_chains.csv",header = TRUE, sep = ",")
df_kim <- read.csv("PLM_Likelihoods/data/Kim/PLMCorrelation_chains.csv",header = TRUE, sep = ",")

#Change sample names
df_ova$sample <- case_match(df_ova$sample,
                            "S1" ~ "Mouse1",
                            "S2" ~ "Mouse2",
                            "S3" ~ "Mouse3",
                            "S4" ~ "Mouse4",
                            "S5" ~ "Mouse5")
df_horns$sample <- case_match(df_horns$sample,
                              "Influenza.vac.11.12.human.S1" ~ "Individual1",
                              "Influenza.vac.11.12.human.S2" ~ "Human2",
                              "Influenza.vac.11.12.human.S3" ~ "Human3",
                              "Influenza.vac.11.12.human.S4" ~ "Human4")
#Only keep 1 replicate
df_horns <- df_horns[df_horns$sample == "Individual1",]
df_bruhn$sample <- "Individual2"
df_kim <- df_kim[df_kim$sample %in% c("SRR17729703", "SRR17729692", "SRR17729726"),]
df_kim$sample <- case_match(df_kim$sample,
                            "SRR17729703" ~ "Individual3",
                            "SRR17729692" ~ "Individual4",
                            "SRR17729726" ~ "Individual5")
df <- rbind(df_ova, df_horns, df_bruhn, df_kim)
df <- pivot_longer(df, cols = 1:6, names_to = "PLM", values_to = "correlation")

#Plot the figures
b <- ggplot(df[df$source == "CDR3 from VDJ",], aes(x=PLM, y=correlation, col=factor(sample))) + 
  geom_point(size = 3) +
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
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  ylim(-0.25,1) +
  xlab("PLM Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~chain) + ggtitle("CDR3 from VDJ")
a <- ggplot(df[df$source == "CDR3 only",], aes(x=PLM, y=correlation, col=factor(sample))) + 
  geom_point(size = 3) +
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
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  ylim(-0.25,1) +
  xlab("PLM Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~chain) + ggtitle("CDR3 only")
c <- ggplot(df[df$source == "Full VDJ",], aes(x=PLM, y=correlation, col=factor(sample))) + 
  geom_point(size = 3) +
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
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  ylim(-0.25,1) +
  xlab("PLM Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~chain) + ggtitle("Full VDJ")

pdf("PLM_Likelihoods/figures/Figure3_Supplementary3/PLMcorrelation_sup.pdf", width = 8, height = 6)
ggarrange(a, b, c, ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()

#Figure 3B
library(ggplot2)
library(ggpubr)
library(Platypus)
library(dplyr)
library(RColorBrewer)
library(gridExtra)

#Read VDJ dataframes
load("PLM_Likelihoods/data/OVA_V7/VDJ_PLL_OVA_V7.RData")
vdj_ova <- vdj
load("PLM_Likelihoods/data/Horns/VDJ_PLL_horns2020a__VDJ_RAW.RData")
vdj_horns <- vdj
load("PLM_Likelihoods/data/Bruhn/VDJ_PLL_Bruhn.RData")
vdj_bruhn <- vdj
rm(vdj)
load("PLM_Likelihoods/data/Kim/VDJ_PLL_Kim.RData")
vdj_kim <- vdj
rm(vdj)

#Change sample names
vdj_ova$sample_id <- case_match(vdj_ova$sample_id,
                                "S1" ~ "Mouse1",
                                "S2" ~ "Mouse2",
                                "S3" ~ "Mouse3",
                                "S4" ~ "Mouse4",
                                "S5" ~ "Mouse5")
vdj_horns$sample_id <- case_match(vdj_horns$sample_id,
                                  "Influenza.vac.11.12.human.S1" ~ "Individual1",
                                  "Influenza.vac.11.12.human.S2" ~ "Human2",
                                  "Influenza.vac.11.12.human.S3" ~ "Human3",
                                  "Influenza.vac.11.12.human.S4" ~ "Human4")
#Only keep 1 replicate
vdj_horns <- vdj_horns[vdj_horns$sample_id == "Individual1",]
vdj_bruhn$sample_id <- "Individual2"
vdj_bruhn <- vdj_bruhn
vdj_kim <- vdj_kim[vdj_kim$sample_id %in% c("SRR17729703", "SRR17729692", "SRR17729726"),]
vdj_kim$sample_id <- case_match(vdj_kim$sample_id,
                                "SRR17729703" ~ "Individual3",
                                "SRR17729692" ~ "Individual4",
                                "SRR17729726" ~ "Individual5")

#Set v-gene families
vdj_ova$v_gene_family <- gsub(pattern = "-.*", replacement = "", x = vdj_ova$VDJ_vgene)
vdj_bruhn$v_gene_family <- gsub(pattern = "-.*", replacement = "", x = vdj_bruhn$VDJ_vgene)
vdj_horns$v_gene_family <- gsub(pattern = "-.*", replacement = "", x = vdj_horns$VDJ_vgene)
vdj_kim$v_gene_family <- gsub(pattern = "-.*", replacement = "", x = vdj_kim$VDJ_vgene)

#Remove NA cgene
vdj_ova <- vdj_ova[!is.na(vdj_ova$VDJ_cgene),]
vdj_kim <- vdj_kim[!is.na(vdj_kim$VDJ_cgene),]
vdj_horns <- vdj_horns[!is.na(vdj_horns$VDJ_cgene),]
vdj_bruhn <- vdj_bruhn[!is.na(vdj_bruhn$VDJ_cgene),]

#combine human samples
vdj_human <- rbind(vdj_horns, vdj_bruhn, vdj_kim)
vdj_all <- rbind(vdj_ova, vdj_human)

#Set colors
isotype_colors <- c("IgA" = "#fb6a4a",
                    "IgE" = "#B452CD",
                    "IgG" ="#74c476",
                    "IgM" = "black",
                    "IgD" = "#1874CD",
                    "NA" = "grey")

#Set isotypes
vdj_ova$isotype <- case_match(vdj_ova$isotype,
                              "IgG1" ~ "IgG",
                              "IgG2" ~ "IgG",
                              "IgG3" ~ "IgG",
                              "IGHA" ~ "IgA",
                              "IGHD" ~ "IgD",
                              "IGHE" ~ "IgE",
                              "IGHM" ~ "IgM")
vdj_ova %>% arrange(factor(isotype, levels = c('IgM', 'IgG', 'IgA', 'IgE', 'IgD'))) -> vdj_ova
vdj_human$isotype <- case_match(vdj_human$isotype,
                                "IgG1" ~ "IgG",
                                "IgG2" ~ "IgG",
                                "IgG3" ~ "IgG",
                                "IgG4" ~ "IgG",
                                "IGHA" ~ "IgA",
                                "IGHD" ~ "IgD",
                                "IGHE" ~ "IgE",
                                "IGHM" ~ "IgM",
                                "IgA1" ~ "IgA",
                                "IgA2" ~ "IgA")
vdj_human %>% arrange(factor(isotype, levels = c('IgM', 'IgG', 'IgA', 'IgE', 'IgD'))) -> vdj_human

#Plot mouse samples
b <- ggplot(vdj_ova, aes(x=IGH_evo_likelihood_esm_full_VDJ, y=IGH_evo_likelihood_protbert_full_VDJ, color=isotype)) +
  geom_point() +
  scale_color_manual(name = "Isotype", values = isotype_colors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_minimal() +
  theme(text = element_text(size = 12)) +
  xlab("ESM-1b Pseudolikelihood") + ylab("ProtBERT Pseudolikelihood") +
  ggtitle("Full VDJ")
a <- ggplot(vdj_ova, aes(x=IGH_evo_likelihood_esm_cdr3_only, y=IGH_evo_likelihood_protbert_cdr3_only, color=isotype)) +
  geom_point() +
  scale_color_manual(name = "Isotype", values = isotype_colors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_minimal() +
  theme(text = element_text(size = 12)) +
  xlab("ESM-1b Pseudolikelihood") + ylab("ProtBERT Pseudolikelihood") +
  ggtitle("CDR3 only")
d <- ggplot(vdj_ova, aes(x=IGH_evo_likelihood_ablang_full_VDJ, y=IGH_evo_likelihood_sapiens_full_VDJ, color=isotype)) +
  geom_point() +
  scale_color_manual(name = "Isotype", values = isotype_colors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_minimal() +
  theme(text = element_text(size = 12)) +
  xlab("Ablang Pseudolikelihood") + ylab("Sapiens Pseudolikelihood") +
  ggtitle("Full VDJ")
c <- ggplot(vdj_ova, aes(x=IGH_evo_likelihood_ablang_cdr3_only, y=IGH_evo_likelihood_sapiens_cdr3_only, color=isotype)) +
  geom_point() +
  scale_color_manual(name = "Isotype", values = isotype_colors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_minimal() +
  theme(text = element_text(size = 12)) +
  xlab("Ablang Pseudolikelihood") + ylab("Sapiens Pseudolikelihood") +
  ggtitle("CDR3 only")

#Plot figure3B top
pdf("PLM_Likelihoods/figures/Figure3_Supplementary3/Isotype_PLMcorrelation_mouse.pdf", width = 8, height = 6)
ggarrange(a, b, c, d, ncol = 4, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()

#Plot human samples
b <- ggplot(vdj_human, aes(x=IGH_evo_likelihood_esm_full_VDJ, y=IGH_evo_likelihood_protbert_full_VDJ, color=isotype)) +
  geom_point() +
  scale_color_manual(name = "Isotype", values = isotype_colors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_minimal() +
  theme(text = element_text(size = 12)) +
  xlab("ESM-1b Pseudolikelihood") + ylab("ProtBERT Pseudolikelihood") +
  ggtitle("Full VDJ")
a <- ggplot(vdj_human, aes(x=IGH_evo_likelihood_esm_cdr3_only, y=IGH_evo_likelihood_protbert_cdr3_only, color=isotype)) +
  geom_point() +
  scale_color_manual(name = "Isotype", values = isotype_colors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_minimal() +
  theme(text = element_text(size = 12)) +
  xlab("ESM-1b Pseudolikelihood") + ylab("ProtBERT Pseudolikelihood") +
  ggtitle("CDR3 only")
d <- ggplot(vdj_human, aes(x=IGH_evo_likelihood_ablang_full_VDJ, y=IGH_evo_likelihood_sapiens_full_VDJ, color=isotype)) +
  geom_point() +
  scale_color_manual(name = "Isotype", values = isotype_colors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_minimal() +
  theme(text = element_text(size = 12)) +
  xlab("Ablang Pseudolikelihood") + ylab("Sapiens Pseudolikelihood") +
  ggtitle("Full VDJ")
c <- ggplot(vdj_human, aes(x=IGH_evo_likelihood_ablang_cdr3_only, y=IGH_evo_likelihood_sapiens_cdr3_only, color=isotype)) +
  geom_point() +
  scale_color_manual(name = "Isotype", values = isotype_colors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_minimal() +
  theme(text = element_text(size = 12)) +
  xlab("Ablang Pseudolikelihood") + ylab("Sapiens Pseudolikelihood") +
  ggtitle("CDR3 only")

#Plot figure3B bottom
pdf("PLM_Likelihoods/figures/Figure3_Supplementary3/Isotype_PLMcorrelation_human.pdf", width = 8, height = 6)
ggarrange(a, b, c, d, ncol = 4, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()

#Figure 3C
#set colors
vgene_colors <- c("IGHV1" = "darkred",
                  "IGHV2" = "cadetblue",
                  "IGHV3" = "darkgreen",
                  "IGHV4" = "deeppink",
                  "IGHV5" = "darkgoldenrod4",
                  "IGHV6" = "snow4",
                  "IGHV7" = "gold3",
                  "IGHV8" = "hotpink",
                  "IGHV9" = "dodgerblue3",
                  "IGHV10" = "violetred4",
                  "IGHV11" = "aquamarine4",
                  "IGHV12" = "pink",
                  "IGHV13" = "lightgreen",
                  "IGHV14" = "orangered",
                  "IGHV15" = "mediumpurple")

#Plot human samples
a1 <- ggplot(vdj_human, aes(x=IGH_evo_likelihood_esm_full_VDJ, y=IGH_evo_likelihood_protbert_full_VDJ, color=v_gene_family)) +
  geom_point(size = 0.1) +
  scale_color_manual(name = "V-gene family", values = vgene_colors,
                     breaks = c("IGHV1", "IGHV2", "IGHV3","IGHV4","IGHV5","IGHV6", "IGHV7","IGHV8",
                                "IGHV9", "IGHV10", "IGHV11", "IGHV12", "IGHV13", "IGHV14", "IGHV15")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_minimal() +
  theme(text = element_text(size = 12)) +
  xlab("ESM-1b Pseudolikelihood") + ylab("ProtBERT Pseudolikelihood") +
  ggtitle("Foundational PLMs")
b1 <- ggplot(vdj_human, aes(x=IGH_evo_likelihood_ablang_full_VDJ, y=IGH_evo_likelihood_sapiens_full_VDJ, color=v_gene_family)) +
  geom_point(size = 0.1) +
  scale_color_manual(name = "V-gene family", values = vgene_colors,
                     breaks = c("IGHV1", "IGHV2", "IGHV3","IGHV4","IGHV5","IGHV6", "IGHV7","IGHV8",
                                "IGHV9", "IGHV10", "IGHV11", "IGHV12", "IGHV13", "IGHV14", "IGHV15")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_minimal() +
  theme(text = element_text(size = 12)) +
  xlab("Ablang Pseudolikelihood") + ylab("Sapiens Pseudolikelihood") +
  ggtitle("Antibody-Specific PLMs")

#Plot figure3C right
pdf("PLM_Likelihoods/figures/Figure3_Supplementary3/Vgene_PLMcorrelation_human.pdf", width = 8, height = 6)
ggarrange(a1, b1, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()

#Plot mouse samples
a2 <- ggplot(vdj_ova, aes(x=IGH_evo_likelihood_esm_full_VDJ, y=IGH_evo_likelihood_protbert_full_VDJ, color=v_gene_family)) +
  geom_point(size = 0.1) +
  scale_color_manual(name = "V-gene family", values = vgene_colors,
                     breaks = c("IGHV1", "IGHV2", "IGHV3","IGHV4","IGHV5","IGHV6", "IGHV7","IGHV8",
                                "IGHV9", "IGHV10", "IGHV11", "IGHV12", "IGHV13", "IGHV14", "IGHV15")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_minimal() +
  theme(text = element_text(size = 12)) +
  xlab("ESM-1b Pseudolikelihood") + ylab("ProtBERT Pseudolikelihood") +
  ggtitle("Foundational PLMs")
b2 <- ggplot(vdj_ova, aes(x=IGH_evo_likelihood_ablang_full_VDJ, y=IGH_evo_likelihood_sapiens_full_VDJ, color=v_gene_family)) +
  geom_point(size = 0.1) +
  scale_color_manual(name = "V-gene family", values = vgene_colors,
                     breaks = c("IGHV1", "IGHV2", "IGHV3","IGHV4","IGHV5","IGHV6", "IGHV7","IGHV8",
                                "IGHV9", "IGHV10", "IGHV11", "IGHV12", "IGHV13", "IGHV14", "IGHV15")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_minimal() +
  theme(text = element_text(size = 12)) +
  xlab("Ablang Pseudolikelihood") + ylab("Sapiens Pseudolikelihood") +
  ggtitle("Antibody-Specific PLMs")

#Plot figure3C left
pdf("PLM_Likelihoods/figures/Figure3_Supplementary3/Vgene_PLMcorrelation_mouse.pdf", width = 8, height = 6)
ggarrange(a2, b2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()
