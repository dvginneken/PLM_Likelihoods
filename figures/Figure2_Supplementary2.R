library(ggplot2) 
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(ggpubr)

#Figure 2B
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

#Plot fig 2b
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

#Supplementary figure 2A
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

#Plot Supplementary Figure2B
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

pdf("PLM_Likelihoods/figures/Figure2_Supplementary2/SourceCorrelation_sup_chains.pdf", width = 8, height = 6)
ggarrange(a, b, c, d, ncol = 2, nrow = 2, common.legend = TRUE, legend = "right")
dev.off()

#Figure 2C
#Read data
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
vdj_ova$sample <- vdj_ova$sample_id
vdj_horns$sample_id <- case_match(vdj_horns$sample_id,
                                  "Influenza.vac.11.12.human.S1" ~ "Individual1",
                                  "Influenza.vac.11.12.human.S2" ~ "Human2",
                                  "Influenza.vac.11.12.human.S3" ~ "Human3",
                                  "Influenza.vac.11.12.human.S4" ~ "Human4")
#Only keep 1 replicate
vdj_horns <- vdj_horns[vdj_horns$sample_id == "Individual1",]
vdj_horns$sample <- vdj_horns$sample_id
vdj_bruhn$sample_id <- "Individual2"
vdj_bruhn$sample <- vdj_bruhn$sample_id
vdj_kim <- vdj_kim[vdj_kim$sample_id %in% c("SRR17729703", "SRR17729692", "SRR17729726"),]
vdj_kim$sample_id <- case_match(vdj_kim$sample_id,
                                "SRR17729703" ~ "Individual3",
                                "SRR17729692" ~ "Individual4",
                                "SRR17729726" ~ "Individual5")
vdj_kim$sample<- vdj_kim$sample_id
#combine human samples
vdj_human <- rbind(vdj_horns, vdj_bruhn, vdj_kim)
vdj_all <- rbind(vdj_ova, vdj_human)

#Ablang
vdj_human$LC_evo_likelihood_ablang_full_VDJ <- vdj_human$IGL_evo_likelihood_ablang_full_VDJ
vdj_human[which(is.na(vdj_human$LC_evo_likelihood_ablang_full_VDJ)),"LC_evo_likelihood_ablang_full_VDJ"] <- vdj_human[which(is.na(vdj_human$LC_evo_likelihood_ablang_full_VDJ)),"IGK_evo_likelihood_ablang_full_VDJ"]
vdj_ova$LC_evo_likelihood_ablang_full_VDJ <- vdj_ova$IGL_evo_likelihood_ablang_full_VDJ
vdj_ova[which(is.na(vdj_ova$LC_evo_likelihood_ablang_full_VDJ)),"LC_evo_likelihood_ablang_full_VDJ"] <- vdj_ova[which(is.na(vdj_ova$LC_evo_likelihood_ablang_full_VDJ)),"IGK_evo_likelihood_ablang_full_VDJ"]
vdj_all <- rbind(vdj_ova, vdj_human)
cor_human <- cor.test(vdj_human$IGH_evo_likelihood_ablang_full_VDJ, vdj_human$LC_evo_likelihood_ablang_full_VDJ)$estimate
cor_mice <- cor.test(vdj_ova$IGH_evo_likelihood_ablang_full_VDJ, vdj_ova$LC_evo_likelihood_ablang_full_VDJ)$estimate
a <-ggplot(vdj_all, aes(x=IGH_evo_likelihood_ablang_full_VDJ, y=LC_evo_likelihood_ablang_full_VDJ, color=sample_id)) +
  geom_point(size = 0.7) +
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
  geom_smooth(method="lm", se = F) +
  theme_minimal() +
  theme(text = element_text(size = 12)) +
  xlab("Heavy Chain Pseudolikelihood") + ylab("Light Chain Pseudolikelihood") +
  ggtitle(paste0("Ablang \nR\u00b2 Individuals = ", round(cor_human, digits = 2), "\nR\u00b2 Mice = ", round(cor_mice, digits = 2)))

#esm
vdj_human$LC_evo_likelihood_esm_full_VDJ <- vdj_human$IGL_evo_likelihood_esm_full_VDJ
vdj_human[which(is.na(vdj_human$LC_evo_likelihood_esm_full_VDJ)),"LC_evo_likelihood_esm_full_VDJ"] <- vdj_human[which(is.na(vdj_human$LC_evo_likelihood_esm_full_VDJ)),"IGK_evo_likelihood_esm_full_VDJ"]
vdj_ova$LC_evo_likelihood_esm_full_VDJ <- vdj_ova$IGL_evo_likelihood_esm_full_VDJ
vdj_ova[which(is.na(vdj_ova$LC_evo_likelihood_esm_full_VDJ)),"LC_evo_likelihood_esm_full_VDJ"] <- vdj_ova[which(is.na(vdj_ova$LC_evo_likelihood_esm_full_VDJ)),"IGK_evo_likelihood_esm_full_VDJ"]
vdj_all <- rbind(vdj_ova, vdj_human)
cor_human <- cor.test(vdj_human$IGH_evo_likelihood_esm_full_VDJ, vdj_human$LC_evo_likelihood_esm_full_VDJ)$estimate
cor_mice <- cor.test(vdj_ova$IGH_evo_likelihood_esm_full_VDJ, vdj_ova$LC_evo_likelihood_esm_full_VDJ)$estimate
b <- ggplot(vdj_all, aes(x=IGH_evo_likelihood_esm_full_VDJ, y=LC_evo_likelihood_esm_full_VDJ, color=sample_id)) +
  geom_point(size = 0.7) +
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
  geom_smooth(method="lm", se = F) +
  theme_minimal() +
  theme(text = element_text(size = 12)) +
  xlab("Heavy Chain Pseudolikelihood") + ylab("Light Chain Pseudolikelihood") +
  ggtitle(paste0("ESM-1b \nR\u00b2 Individuals = ", round(cor_human, digits = 2), "\nR\u00b2 Mice = ", round(cor_mice, digits = 2)))

pdf("PLM_Likelihoods/figures/Figure2_Supplementary2/ChainCorrelation_main.pdf", width = 8, height = 4)
ggarrange(a, b, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()

#Supplementary figure 2C
#Sapiens
vdj_human$LC_evo_likelihood_sapiens_full_VDJ <- vdj_human$IGL_evo_likelihood_sapiens_full_VDJ
vdj_human[which(is.na(vdj_human$LC_evo_likelihood_sapiens_full_VDJ)),"LC_evo_likelihood_sapiens_full_VDJ"] <- vdj_human[which(is.na(vdj_human$LC_evo_likelihood_sapiens_full_VDJ)),"IGK_evo_likelihood_sapiens_full_VDJ"]
vdj_ova$LC_evo_likelihood_sapiens_full_VDJ <- vdj_ova$IGL_evo_likelihood_sapiens_full_VDJ
vdj_ova[which(is.na(vdj_ova$LC_evo_likelihood_sapiens_full_VDJ)),"LC_evo_likelihood_sapiens_full_VDJ"] <- vdj_ova[which(is.na(vdj_ova$LC_evo_likelihood_sapiens_full_VDJ)),"IGK_evo_likelihood_sapiens_full_VDJ"]
vdj_all <- rbind(vdj_ova, vdj_human)
cor_human <- cor.test(vdj_human$IGH_evo_likelihood_sapiens_full_VDJ, vdj_human$LC_evo_likelihood_sapiens_full_VDJ)$estimate
cor_mice <- cor.test(vdj_ova$IGH_evo_likelihood_sapiens_full_VDJ, vdj_ova$LC_evo_likelihood_sapiens_full_VDJ)$estimate
c <- ggplot(vdj_all, aes(x=IGH_evo_likelihood_sapiens_full_VDJ, y=LC_evo_likelihood_sapiens_full_VDJ, color=sample_id)) +
  geom_point(size = 1) +
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
  geom_smooth(method="lm", se = F) +
  theme_minimal() +
  theme(text = element_text(size = 12)) +
  xlab("Heavy Chain Pseudolikelihood") + ylab("Light Chain Pseudolikelihood") +
  ggtitle(paste0("Sapiens \nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mice, digits = 2)))

#protbert
vdj_human$LC_evo_likelihood_protbert_full_VDJ <- vdj_human$IGL_evo_likelihood_protbert_full_VDJ
vdj_human[which(is.na(vdj_human$LC_evo_likelihood_protbert_full_VDJ)),"LC_evo_likelihood_protbert_full_VDJ"] <- vdj_human[which(is.na(vdj_human$LC_evo_likelihood_protbert_full_VDJ)),"IGK_evo_likelihood_protbert_full_VDJ"]
vdj_ova$LC_evo_likelihood_protbert_full_VDJ <- vdj_ova$IGL_evo_likelihood_protbert_full_VDJ
vdj_ova[which(is.na(vdj_ova$LC_evo_likelihood_protbert_full_VDJ)),"LC_evo_likelihood_protbert_full_VDJ"] <- vdj_ova[which(is.na(vdj_ova$LC_evo_likelihood_protbert_full_VDJ)),"IGK_evo_likelihood_protbert_full_VDJ"]
vdj_all <- rbind(vdj_ova, vdj_human)
cor_human <- cor.test(vdj_human$IGH_evo_likelihood_protbert_full_VDJ, vdj_human$LC_evo_likelihood_protbert_full_VDJ)$estimate
cor_mice <- cor.test(vdj_ova$IGH_evo_likelihood_protbert_full_VDJ, vdj_ova$LC_evo_likelihood_protbert_full_VDJ)$estimate
d <- ggplot(vdj_all, aes(x=IGH_evo_likelihood_protbert_full_VDJ, y=LC_evo_likelihood_protbert_full_VDJ, color=sample_id)) +
  geom_point(size = 1) +
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
  geom_smooth(method="lm", se = F) +
  theme_minimal() +
  theme(text = element_text(size = 12)) +
  xlab("Heavy Chain Pseudolikelihood") + ylab("Light Chain Pseudolikelihood") +
  ggtitle(paste0("ProtBERT \nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mice, digits = 2)))

#Plot supplementary figure 2C
pdf("PLM_Likelihoods/figures/Figure2_Supplementary2/ChainCorrelation_sup.pdf", width = 8, height = 4)
ggarrange(a, b, c, d, ncol = 4, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()