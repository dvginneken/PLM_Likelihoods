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

## Supplementary Figure 1
source("PLM_Likelihoods/scripts/VDJ_clonal_expansion.R")
VDJ_clonal_expansion(vdj_bruhn, clones = 30, group.by = "sample_id", color.by = "VDJ_cgene", text.size=20)
VDJ_clonal_expansion(vdj_horns, clones = 30, group.by = "sample_id", color.by = "VDJ_cgene", text.size=20)
VDJ_clonal_expansion(vdj_ova, clones = 30, group.by = "sample_id", color.by = "VDJ_cgene", text.size=20)
VDJ_clonal_expansion(vdj_kim, clones = 30, group.by = "sample_id", color.by = "VDJ_cgene", text.size=20)

## Figure1B

#Bottom
#https://github.com/alexyermanos/Platypus/blob/master/R/VDJ_clonal_barplot.R
vdj_human$sample_id <- gsub(pattern = "Individual", replacement = "I", x = vdj_human$sample_id)
vdj_ova$sample_id <- gsub(pattern = "Mouse", replacement = "M", x = vdj_ova$sample_id)
out1 <- Platypus::VDJ_clonal_barplot(vdj_human, counts.to.use = "clonotype_id", group.by = "sample_id", expanded.colors = brewer.pal(n=5, "BuGn")[2:5])
out2 <- Platypus::VDJ_clonal_barplot(vdj_ova, counts.to.use = "clonotype_id", group.by = "sample_id", expanded.colors = brewer.pal(n=5, "RdPu")[2:5])
pdf("PLM_Likelihoods/figures/Figure1_Supplementary1/ClonalExpansion_barplot.pdf")
do.call("grid.arrange", c(out1, out2, ncol = 10))
dev.off()

#Top
vdj_all <- rbind(vdj_human, vdj_ova)
pdf("PLM_Likelihoods/figures/Figure1_Supplementary1/CellCount_barplot.pdf",
    height = 2)
ggplot(vdj_all, aes(x = sample_id, fill=sample_id)) +
  geom_bar(position = "dodge", stat = "count") +
  theme_minimal() +
  scale_fill_manual(values = c("I1" = "#99d8c9",
                               "I2" = "#66c2a4",
                               "I3" = "#41ae76",
                               "I4" = "#238b45",
                               "I5" = "#005824",
                               "M1" = "#fcc5c0",
                               "M2" = "#fa9fb5",
                               "M3" = "#f768a1",
                               "M4" = "#c51b8a",
                               "M5" = "#7a0177")) +
  theme(text = element_text(size = 15),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("Sample") + ylab("Number of cells")
dev.off()

## Supplementary Figure 2
vdj_all$VDJ_cgene <- gsub("IGHA1", "IGHA", vdj_all$VDJ_cgene)
vdj_all$VDJ_cgene <- gsub("IGHA2", "IGHA", vdj_all$VDJ_cgene)
vdj_all$VDJ_cgene <- gsub("IGHG1", "IGHG", vdj_all$VDJ_cgene)
vdj_all$VDJ_cgene <- gsub("IGHG2", "IGHG", vdj_all$VDJ_cgene)
vdj_all$VDJ_cgene <- gsub("IGHG2B", "IGHG", vdj_all$VDJ_cgene)
vdj_all$VDJ_cgene <- gsub("IGHG2C", "IGHG", vdj_all$VDJ_cgene)
vdj_all$VDJ_cgene <- gsub("IGHG3", "IGHG", vdj_all$VDJ_cgene)
vdj_all$VDJ_cgene <- gsub("IGHG4", "IGHG", vdj_all$VDJ_cgene)
vdj_all$VDJ_cgene <- gsub("IGHGB", "IGHG", vdj_all$VDJ_cgene)
vdj_all$VDJ_cgene <- gsub("IGHGC", "IGHG", vdj_all$VDJ_cgene)

cl_colors <- c("IGHA" = "#fb6a4a",
               "IGHE" = "#B452CD",
               "IGHG" = "#74c476",
               "IGHM" = "black",
               "IGHD" = "#1874CD",
               "NA" = "grey")
pdf("PLM_Likelihoods/figures/Figure1_Supplementary1/IsotypeDistribution_barplot.pdf",
    height = 5)
ggplot(vdj_all, aes(x=sample_id, fill=VDJ_cgene)) +
  geom_bar(position = "fill", stat = "count") +
  theme_minimal() +
  theme(text = element_text(size = 15),
        axis.ticks.x=element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = cl_colors, name = "Isotype") +
  xlab("Sample") + ylab("Percentage of cells")
dev.off()

#Souce Correlation
#Read data
df_ova <- read.csv("PLM_Likelihoods/data/OVA_V7/SourceCorrelation.csv",header = TRUE, sep = ",")
df_horns <- read.csv("PLM_Likelihoods/data/Horns/SourceCorrelation.csv",header = TRUE, sep = ",")
df_bruhn <- read.csv("PLM_Likelihoods/data/Bruhn/SourceCorrelation.csv",header = TRUE, sep = ",")
df_kim <- read.csv("PLM_Likelihoods/data/Kim/SourceCorrelation.csv",header = TRUE, sep = ",")

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
df_kim$sample <- case_match(df_kim$sample,
                            "SRR17729703" ~ "Individual3",
                            "SRR17729692" ~ "Individual4",
                            "SRR17729726" ~ "Individual5")
df_kim <- df_kim[df_kim$sample %in% paste0("Individual",2:5),]
df <- rbind(df_ova, df_horns, df_bruhn, df_kim)

#Tranform shape of dataframe
df <- pivot_longer(df, cols = 1:3, names_to = "source", values_to = "correlation")
df$source<- case_match(df$source,
                       "full_VDJ__CDR3_only" ~ "VDJ\nCDR3",
                       "full_VDJ__CDR3_from_VDJ" ~ "VDJ\nCDR3-VDJ",
                       "CDR3_only__CDR3_from_VDJ" ~ "CDR3\nCDR3-VDJ")
df_main <- df[df$model %in% c("ESM-1b", "Ablang"),]

#Figure 1 D
pdf("PLM_Likelihoods/figures/Figure2_Supplementary2/SourceCorrelation_main.pdf", width = 8, height = 6)
ggplot(df_main, aes(x=source, y=correlation, col=factor(sample))) + 
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
  theme(text = element_text(size = 18)) +
  xlab("Source Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~model)
dev.off()

#Supplementary figure 3A
pdf("PLM_Likelihoods/figures/Figure2_Supplementary2/SourceCorrelation_sup.pdf", width = 8, height = 6)
ggplot(df, aes(x=source, y=correlation, col=factor(sample))) + 
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
  theme(text = element_text(size = 14)) +
  xlab("Source Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~model, nrow = 1, ncol = 4)
dev.off()

#Supplementary figure 2B
#Read data
df_ova <- read.csv("PLM_Likelihoods/data/OVA_V7/SourceCorrelation_chains.csv",
                   header = TRUE, sep = ",")
df_horns <- read.csv("PLM_Likelihoods/data/Horns/SourceCorrelation_chains.csv",
                     header = TRUE, sep = ",")
df_bruhn <- read.csv("PLM_Likelihoods/data/Bruhn/SourceCorrelation_chains.csv",
                     header = TRUE, sep = ",")
df_kim <- read.csv("PLM_Likelihoods/data/Kim/SourceCorrelation_chains.csv",
                   header = TRUE, sep = ",")

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
df_kim$sample <- case_match(df_kim$sample,
                            "SRR17729703" ~ "Individual3",
                            "SRR17729692" ~ "Individual4",
                            "SRR17729726" ~ "Individual5")
df_kim <- df_kim[df_kim$sample %in% paste0("Individual",2:5),]

#Tranform shape of dataframe
df <- rbind(df_ova, df_horns, df_bruhn, df_kim)
df <- pivot_longer(df, cols = 1:3, names_to = "source", values_to = "correlation")
df$source<- case_match(df$source,
                       "full_VDJ__CDR3_only" ~ "VDJ\nCDR3",
                       "full_VDJ__CDR3_from_VDJ" ~ "VDJ\nCDR3-VDJ",
                       "CDR3_only__CDR3_from_VDJ" ~ "CDR3\nCDR3-VDJ")

#Plot Supplementary Figure3B
a <- ggplot(df[df$model == "ESM-1b",], aes(x=source, y=correlation, col=factor(sample))) + 
  geom_point(size = 2) +
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
  theme(text = element_text(size = 14)) +
  xlab("Source Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~chain) + ggtitle("ESM-1b")
b <- ggplot(df[df$model == "ProtBERT",], aes(x=source, y=correlation, col=factor(sample))) + 
  geom_point(size = 2) +
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
  theme(text = element_text(size = 14)) +
  xlab("Source Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~chain) + ggtitle("ProtBERT")
c <- ggplot(df[df$model == "Ablang",], aes(x=source, y=correlation, col=factor(sample))) + 
  geom_point(size = 2) +
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
  theme(text = element_text(size = 14)) +
  xlab("Source Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~chain) + ggtitle("Ablang")
d <- ggplot(df[df$model == "Sapiens",], aes(x=source, y=correlation, col=factor(sample))) + 
  geom_point(size = 2) +
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
  theme(text = element_text(size = 14)) +
  xlab("Source Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~chain) + ggtitle("Sapiens")
