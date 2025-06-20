library(ggplot2)
library(ggpubr)
library(Platypus)
library(dplyr)
library(RColorBrewer)
library(gridExtra)
library(tidyr)

source("~/OneDrive - UMC Utrecht/Documenten/Steropodon_theme.R")

#Read VDJ dataframes
load("~/Documents/GitHub/PLM-likelihoods/data/OVA_V7/VDJ_PLL_OVA_V7.RData")
vdj_ova <- vdj_likelihood
load("~/Documents/GitHub/PLM-likelihoods/data/horns2020a__VDJ_RAW/VDJ_PLL_horns2020a__VDJ_RAW.RData")
vdj_horns <- vdj_likelihood
load("~/Documents/GitHub/PLM-likelihoods/data/Bruhn/VDJ_PLL_Bruhn.RData")
vdj_bruhn <- vdj_likelihood
load("~/Documents/GitHub/PLM-likelihoods/data/Kim/VDJ_PLL_Kim.RData")
vdj_kim <- vdj_likelihood
rm(vdj_likelihood)

#Change sample names
vdj_ova$sample_id <- case_match(vdj_ova$sample_id,
                                "S1" ~ "Mouse1",
                                "S2" ~ "Mouse2",
                                "S3" ~ "Mouse3",
                                "S4" ~ "Mouse4",
                                "S5" ~ "Mouse5")
vdj_horns$sample_id <- "Individual1"
vdj_bruhn$sample_id <- "Individual2"
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
source("PLM-likelihoods/scripts/VDJ_clonal_expansion.R")
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
#Figure 1 and Supplementals 
#Read data
df_ova_hc <- read.csv("~/Documents/GitHub/PLM-likelihoods/data/OVA_V7/SourceCorrelation_HC.csv",
                   header = TRUE, sep = ",")
df_ova_lc <- read.csv("~/Documents/GitHub/PLM-likelihoods/data/OVA_V7/SourceCorrelation_LC.csv",
                   header = TRUE, sep = ",")
df_horns_hc <- read.csv("~/Documents/GitHub/PLM-likelihoods/data/horns2020a__VDJ_RAW/SourceCorrelation_HC.csv",
                     header = TRUE, sep = ",")
df_horns_lc <- read.csv("~/Documents/GitHub/PLM-likelihoods/data/horns2020a__VDJ_RAW/SourceCorrelation_LC.csv",
                     header = TRUE, sep = ",")
df_bruhn_hc <- read.csv("~/Documents/GitHub/PLM-likelihoods/data/Bruhn/SourceCorrelation_HC.csv",
                     header = TRUE, sep = ",")
df_bruhn_lc <- read.csv("~/Documents/GitHub/PLM-likelihoods/data/Bruhn/SourceCorrelation_LC.csv",
                     header = TRUE, sep = ",")
df_kim_hc <- read.csv("~/Documents/GitHub/PLM-likelihoods/data/Kim/SourceCorrelation_HC.csv",
                   header = TRUE, sep = ",")
df_kim_lc <- read.csv("~/Documents/GitHub/PLM-likelihoods/data/Kim/SourceCorrelation_LC.csv",
                   header = TRUE, sep = ",")

#Change sample names
df_ova_hc$sample <- case_match(df_ova_hc$sample,
                            "S1" ~ "Mouse1",
                            "S2" ~ "Mouse2",
                            "S3" ~ "Mouse3",
                            "S4" ~ "Mouse4",
                            "S5" ~ "Mouse5")
df_ova_lc$sample <- case_match(df_ova_lc$sample,
                            "S1" ~ "Mouse1",
                            "S2" ~ "Mouse2",
                            "S3" ~ "Mouse3",
                            "S4" ~ "Mouse4",
                            "S5" ~ "Mouse5")
df_horns_hc$sample <- "Individual1"
df_horns_lc$sample <- "Individual1"
df_bruhn_hc$sample <- "Individual2"
df_bruhn_lc$sample <- "Individual2"
df_kim_hc$sample <- case_match(df_kim_hc$sample,
                            "SRR17729703" ~ "Individual3",
                            "SRR17729692" ~ "Individual4",
                            "SRR17729726" ~ "Individual5")
df_kim_lc$sample <- case_match(df_kim_lc$sample,
                            "SRR17729703" ~ "Individual3",
                            "SRR17729692" ~ "Individual4",
                            "SRR17729726" ~ "Individual5")
df_hc <- rbind(df_ova_hc, df_horns_hc, df_bruhn_hc, df_kim_hc)
df_lc <- rbind(df_ova_lc, df_horns_lc, df_bruhn_lc, df_kim_lc)
colnames(df_hc) <- c("full_VDJ__CDR3_only", "full_VDJ__CDR3_from_VDJ", "CDR3_only__CDR3_from_VDJ", "sample", "chain", "model")
colnames(df_lc) <- c("full_VDJ__CDR3_only", "full_VDJ__CDR3_from_VDJ", "CDR3_only__CDR3_from_VDJ", "sample", "chain", "model")


#Tranform shape of dataframe
df_hc <- pivot_longer(df_hc, cols = 1:3, names_to = "source", values_to = "correlation")
df_hc$source<- case_match(df_hc$source,
                       "full_VDJ__CDR3_only" ~ "VDJ\nCDR3",
                       "full_VDJ__CDR3_from_VDJ" ~ "VDJ\nCDR3-VDJ",
                       "CDR3_only__CDR3_from_VDJ" ~ "CDR3\nCDR3-VDJ")

df_main <- df_hc[df_hc$model %in% c("ESM-C", "ESM-1b", "Ablang1", "Ablang2"),]

#Figure 1C
pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/Figure1/SourceCorrelation_main.pdf", width = 20, height = 6)
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
  theme_steropodon() +
  theme(text = element_text(size = 26)) +
  xlab("Source Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~model, nrow = 1, ncol = 4)
dev.off()

#Supplementary figure 3A
df_sup <- df_hc[df_hc$model %in% c("ProtBERT", "Sapiens"),]
pdf("PLM-likelihoods/figures/figureS1/SourceCorrelation_sup.pdf", width = 12, height = 6)
ggplot(df_sup, aes(x=source, y=correlation, col=factor(sample))) + 
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
  theme_steropodon() +
  theme(text = element_text(size = 20)) +
  xlab("Source Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~model, ncol = 2)
dev.off()


#Plot Supplementary Figure3B
#Tranform shape of dataframe
df_lc <- pivot_longer(df_lc, cols = 1:3, names_to = "source", values_to = "correlation")
df_lc$source<- case_match(df_lc$source,
                          "full_VDJ__CDR3_only" ~ "VDJ\nCDR3",
                          "full_VDJ__CDR3_from_VDJ" ~ "VDJ\nCDR3-VDJ",
                          "CDR3_only__CDR3_from_VDJ" ~ "CDR3\nCDR3-VDJ")

pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figureS3/SourceCorrelation_supB.pdf", width = 16, height = 9)
ggplot(df_lc, aes(x=source, y=correlation, col=factor(sample))) + 
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
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  xlab("Source Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~model, ncol = 3, nrow = 2)
dev.off()

#statistics
df <- rbind(df_main, df_sup)
mean(df[df$source == "VDJ\nCDR3",]$correlation)
mean(df[df$source == "CDR3\nCDR3-VDJ",]$correlation)
mean(df[df$source == "VDJ\nCDR3-VDJ",]$correlation)

df_human <- df[grep("Ind", df$sample),]
mean(df_human[df_human$source == "VDJ\nCDR3-VDJ" & df_human$model %in% c("Ablang1"),]$correlation)
mean(df_human[df_human$source == "VDJ\nCDR3-VDJ" & df_human$model %in% c("Ablang2"),]$correlation)
