library(ggplot2) 
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
library(gridExtra)

source("~/OneDrive - UMC Utrecht/Documenten/Steropodon_theme.R")


#Plot figure 3A
#Read data
df_ova_HC <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/OVA_V7/PLMCorrelation_HC.csv",header = TRUE, sep = ",")
df_horns_HC <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/horns2020a__VDJ_RAW/PLMCorrelation_HC.csv",header = TRUE, sep = ",")
df_bruhn_HC <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Bruhn/PLMCorrelation_HC.csv",header = TRUE, sep = ",")
df_kim_HC <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim/PLMCorrelation_HC.csv",header = TRUE, sep = ",")

#Change sample names
df_ova_HC$sample <- case_match(df_ova_HC$sample,
                            "S1" ~ "Mouse1",
                            "S2" ~ "Mouse2",
                            "S3" ~ "Mouse3",
                            "S4" ~ "Mouse4",
                            "S5" ~ "Mouse5")
df_horns_HC$sample <- "Individual1"
df_bruhn_HC$sample <- "Individual2"
df_kim_HC$sample <- case_match(df_kim_HC$sample,
                            "SRR17729703" ~ "Individual3",
                            "SRR17729692" ~ "Individual4",
                            "SRR17729726" ~ "Individual5")
df_HC <- rbind(df_ova_HC, df_horns_HC, df_bruhn_HC, df_kim_HC)
df_HC <- pivot_longer(df_HC, cols = 1:15, names_to = "PLM", values_to = "correlation")
df_HC$PLM <- gsub(pattern = "HC_", replacement = "", x = df_HC$PLM)

df_HC_main <- df_HC[grep("Ablang1", df_HC$PLM, invert = T),]
df_HC_main <- df_HC_main[grep("ESM1b", df_HC_main$PLM, invert = T),]

#Statistics
df_human <- df_HC[grep("Ind", df_HC$sample),]
df_human %>% group_by(source) %>% summarize(mean(correlation))
df_human %>% group_by(PLM) %>% summarize(mean(correlation))
df_human %>% summarize(mean(correlation))

df_vdj <- df_human[df_human$source == "Full VDJ",]
df_vdj %>% group_by(PLM) %>% summarize(mean(correlation))

df_ab <- df_human[df_human$PLM %in% c("Ablang1_Sapiens", "Ablang2_Sapiens", "Ablang2_Ablang1"),]
df_ab <- df_ab[df_ab$source == "Full VDJ",]
df_ab %>% group_by(PLM) %>% summarize(mean(correlation))
df_ab %>% summarize(mean(correlation))

df_mouse <- df_HC[grep("Mouse", df_HC$sample),]
df_mouse %>% group_by(source) %>% summarize(mean(correlation))
df_mouse %>% group_by(PLM) %>% summarize(mean(correlation))
df_mouse %>% summarize(mean(correlation))

df_vdj <- df_mouse[df_mouse$source == "Full VDJ",]
df_vdj %>% group_by(PLM) %>% summarize(mean(correlation))

df_ab <- df_vdj[df_vdj$PLM %in% c("ESMc_ESM1b", "ESMc_ProtBERT", "ProtBERT_ESM1b"),]
df_mouse %>% group_by(PLM) %>% summarize(mean(correlation))
df_ab %>% summarize(mean(correlation))

#main
pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figure2/PLMcorrelation_main.pdf", width = 20, height = 8)
ggplot(df_HC_main, aes(x=PLM, y=correlation, col=factor(sample))) + 
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
  theme(text = element_text(size = 26),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  xlab("PLM Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~source)
dev.off()

#supplementary
df_HC_sup <- df_HC[grep("Ablang1", df_HC$PLM),]
df_HC_sup <- rbind(df_HC_sup, df_HC[grep("ESM1b", df_HC$PLM),])

pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figureS4/PLMcorrelation_sup.pdf", width = 25, height = 6)
ggplot(df_HC_sup, aes(x=PLM, y=correlation, col=factor(sample))) + 
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
  theme(text = element_text(size = 25),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  xlab("PLM Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~source)
dev.off()

#Figure S4
#Read data
df_ova_LC <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/OVA_V7/PLMCorrelation_LC.csv",header = TRUE, sep = ",")
df_horns_LC <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/horns2020a__VDJ_RAW/PLMCorrelation_LC.csv",header = TRUE, sep = ",")
df_bruhn_LC <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Bruhn/PLMCorrelation_LC.csv",header = TRUE, sep = ",")
df_kim_LC <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim/PLMCorrelation_LC.csv",header = TRUE, sep = ",")

#Change sample names
df_ova_LC$sample <- case_match(df_ova_LC$sample,
                               "S1" ~ "Mouse1",
                               "S2" ~ "Mouse2",
                               "S3" ~ "Mouse3",
                               "S4" ~ "Mouse4",
                               "S5" ~ "Mouse5")
df_horns_LC$sample <- "Individual1"
df_bruhn_LC$sample <- "Individual2"
df_kim_LC$sample <- case_match(df_kim_LC$sample,
                               "SRR17729703" ~ "Individual3",
                               "SRR17729692" ~ "Individual4",
                               "SRR17729726" ~ "Individual5")
df_LC <- rbind(df_ova_LC, df_horns_LC, df_bruhn_LC, df_kim_LC)
df_LC <- pivot_longer(df_LC, cols = 1:15, names_to = "PLM", values_to = "correlation")
df_LC$PLM <- gsub(pattern = "LC_", replacement = "", x = df_LC$PLM)

pdf("PLM-likelihoods/figures/figureS4/PLMcorrelation_supB.pdf", width = 25, height = 6)
ggplot(df_LC, aes(x=PLM, y=correlation, col=factor(sample))) + 
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
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  xlab("PLM Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~source)
dev.off()


#Figure 2B
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

#Statistics
mean(vdj_ova[vdj_ova$v_gene_family == "IGHV1",c("HC_evo_likelihood_esmc_full_VDJ")])
mean(vdj_ova[vdj_ova$v_gene_family == "IGHV1",c("HC_evo_likelihood_esm1b_full_VDJ")])
mean(vdj_ova[vdj_ova$v_gene_family == "IGHV1",c("HC_evo_likelihood_protbert_full_VDJ")])

mean(vdj_ova[vdj_ova$v_gene_family != "IGHV1",c("HC_evo_likelihood_esmc_full_VDJ")])
mean(vdj_ova[vdj_ova$v_gene_family != "IGHV1",c("HC_evo_likelihood_esm1b_full_VDJ")])
mean(vdj_ova[vdj_ova$v_gene_family != "IGHV1",c("HC_evo_likelihood_protbert_full_VDJ")])

mean(vdj_human[vdj_human$v_gene_family %in% c("IGHV1", "IGHV4"), c("HC_evo_likelihood_esmc_full_VDJ")])
mean(vdj_human[vdj_human$v_gene_family %in% c("IGHV1", "IGHV4"),c("HC_evo_likelihood_esm1b_full_VDJ")])
mean(vdj_human[vdj_human$v_gene_family %in% c("IGHV1", "IGHV4"),c("HC_evo_likelihood_protbert_full_VDJ")])

mean(vdj_human[!vdj_human$v_gene_family %in% c("IGHV1", "IGHV4"),c("HC_evo_likelihood_esmc_full_VDJ")])
mean(vdj_human[!vdj_human$v_gene_family %in% c("IGHV1", "IGHV4"),c("HC_evo_likelihood_esm1b_full_VDJ")])
mean(vdj_human[!vdj_human$v_gene_family %in% c("IGHV1", "IGHV4"),c("HC_evo_likelihood_protbert_full_VDJ")])

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
a1 <- ggplot(vdj_human, aes(x=HC_evo_likelihood_esmc_full_VDJ, y=HC_evo_likelihood_protbert_full_VDJ, color=v_gene_family)) +
  geom_point(size = 0.1) +
  scale_color_manual(name = "V-gene family", values = vgene_colors,
                     breaks = c("IGHV1", "IGHV2", "IGHV3","IGHV4","IGHV5","IGHV6", "IGHV7","IGHV8",
                                "IGHV9", "IGHV10", "IGHV11", "IGHV12", "IGHV13", "IGHV14", "IGHV15")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_steropodon() +
  theme(text = element_text(size = 18)) +
  xlab("ESM-C SP") + ylab("ProtBERT SP")
b1 <- ggplot(vdj_human, aes(x=paired_evo_likelihood_ablang2_full_VDJ, y=HC_evo_likelihood_sapiens_full_VDJ, color=v_gene_family)) +
  geom_point(size = 0.1) +
  scale_color_manual(name = "V-gene family", values = vgene_colors,
                     breaks = c("IGHV1", "IGHV2", "IGHV3","IGHV4","IGHV5","IGHV6", "IGHV7","IGHV8",
                                "IGHV9", "IGHV10", "IGHV11", "IGHV12", "IGHV13", "IGHV14", "IGHV15")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_steropodon() +
  theme(text = element_text(size = 18)) +
  xlab("Ablang2 SP") + ylab("Sapiens SP")

#Plot figure2B right
pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figure2/Vgene_PLMcorrelation_human.pdf", width = 10, height = 4)
ggarrange(a1, b1, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()

#Plot mouse samples
a2 <- ggplot(vdj_ova, aes(x=HC_evo_likelihood_esmc_full_VDJ, y=HC_evo_likelihood_protbert_full_VDJ, color=v_gene_family)) +
  geom_point(size = 0.1) +
  scale_color_manual(name = "V-gene family", values = vgene_colors,
                     breaks = c("IGHV1", "IGHV2", "IGHV3","IGHV4","IGHV5","IGHV6", "IGHV7","IGHV8",
                                "IGHV9", "IGHV10", "IGHV11", "IGHV12", "IGHV13", "IGHV14", "IGHV15")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_steropodon() +
  theme(text = element_text(size = 18)) +
  xlab("ESM-C SP") + ylab("ProtBERT SP")
b2 <- ggplot(vdj_ova, aes(x=paired_evo_likelihood_ablang2_full_VDJ, y=HC_evo_likelihood_sapiens_full_VDJ, color=v_gene_family)) +
  geom_point(size = 0.1) +
  scale_color_manual(name = "V-gene family", values = vgene_colors,
                     breaks = c("IGHV1", "IGHV2", "IGHV3","IGHV4","IGHV5","IGHV6", "IGHV7","IGHV8",
                                "IGHV9", "IGHV10", "IGHV11", "IGHV12", "IGHV13", "IGHV14", "IGHV15")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_steropodon() +
  theme(text = element_text(size = 18)) +
  xlab("Ablang2 SP") + ylab("Sapiens SP")

#Plot figure2B left
pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figure2/Vgene_PLMcorrelation_mouse.pdf", width = 10, height = 4)
ggarrange(a2, b2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()


#Extra scatterplots
#ablang2-paired vs ablang2 heavy chain
cor_ablang <- cor.test(vdj_all$paired_evo_likelihood_ablang2_full_VDJ, vdj_all$HC_evo_likelihood_ablang2_full_VDJ)$estimate
c1 <- ggplot(vdj_all, aes(x=paired_evo_likelihood_ablang2_full_VDJ, y=HC_evo_likelihood_ablang2_full_VDJ, color=sample_id)) +
  geom_point(size = 0.1) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_steropodon() +
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
  theme(text = element_text(size = 12)) +
  xlab("Ablang2 paired SP") + ylab("Ablang2 HC SP") +
  ggtitle(paste0("R\u00b2 = ", round(cor_ablang, digits = 2)))

#ablang2-paired vs ablang2 heavy chain
cor_ablang <- cor.test(vdj_all$paired_evo_likelihood_ablang2_full_VDJ, vdj_all$LC_evo_likelihood_ablang2_full_VDJ)$estimate
c2 <- ggplot(vdj_all, aes(x=paired_evo_likelihood_ablang2_full_VDJ, y=LC_evo_likelihood_ablang2_full_VDJ, color=sample_id)) +
  geom_point(size = 0.1) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_steropodon() +
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
  theme(text = element_text(size = 12)) +
  xlab("Ablang2 paired SP") + ylab("Ablang2 LC SP") +
  ggtitle(paste0("R\u00b2 = ", round(cor_ablang, digits = 2)))

# #ESM1b vs. ESMC
# cor_esm <- cor.test(vdj_all$HC_evo_likelihood_esmc_full_VDJ, vdj_all$HC_evo_likelihood_esm1b_full_VDJ)$estimate
# c2 <- ggplot(vdj_all, aes(x=HC_evo_likelihood_esmc_full_VDJ, y=HC_evo_likelihood_esm1b_full_VDJ, color=sample_id)) +
#   geom_point(size = 0.1) +
#   guides(color = guide_legend(override.aes = list(size = 3))) +
#   theme_steropodon() +
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
#   theme(text = element_text(size = 12)) +
#   xlab("ESM-C SP") + ylab("ESM-1b SP") +
#   ggtitle(paste0("R\u00b2 = ", round(cor_esm, digits = 2)))
# 
# #Ablang2 paired vs. ESMC
# cor_esm_ablang <- cor.test(vdj_all$paired_evo_likelihood_ablang2_full_VDJ, vdj_all$HC_evo_likelihood_esmc_full_VDJ)$estimate
# c3 <- ggplot(vdj_all, aes(x=paired_evo_likelihood_ablang2_full_VDJ, y=HC_evo_likelihood_esmc_full_VDJ, color=sample_id)) +
#   geom_point(size = 0.1) +
#   guides(color = guide_legend(override.aes = list(size = 3))) +
#   theme_steropodon() +
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
#   theme(text = element_text(size = 12)) +
#   xlab("Ablang2 paired SP") + ylab("ESM-C SP") +
#   ggtitle(paste0("R\u00b2 = ", round(cor_esm_ablang, digits = 2)))
# 
# #Ablang1 vs. Ablang2 paired
# cor_ablang1_ablang2 <- cor.test(vdj_all$paired_evo_likelihood_ablang2_full_VDJ, vdj_all$HC_evo_likelihood_ablang1_full_VDJ)$estimate
# c4 <- ggplot(vdj_all, aes(x=paired_evo_likelihood_ablang2_full_VDJ, y=HC_evo_likelihood_ablang1_full_VDJ, color=sample_id)) +
#   geom_point(size = 0.1) +
#   guides(color = guide_legend(override.aes = list(size = 3))) +
#   theme_steropodon() +
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
#   theme(text = element_text(size = 12)) +
#   xlab("Ablang2 paired SP") + ylab("Ablang1 SP") +
#   ggtitle(paste0("R\u00b2 = ", round(cor_ablang1_ablang2, digits = 2)))

pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figure2/PLMcorrelation_extra.pdf", width = 10, height = 4)
ggarrange(c1, c2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()

#Correlation heavy chain / light chain
cor_df = data.frame()
for(sample in unique(vdj_all$sample_id)) {
  df <- vdj_all[vdj_all$sample_id == sample,]
  ablang1 <- cor.test(df$HC_evo_likelihood_ablang1_full_VDJ, df$LC_evo_likelihood_ablang1_full_VDJ, method = "pearson")$estimate
  ablang2 <- cor.test(df$HC_evo_likelihood_ablang2_full_VDJ, df$LC_evo_likelihood_ablang2_full_VDJ, method = "pearson")$estimate
  esmc <- cor.test(df$HC_evo_likelihood_esmc_full_VDJ, df$LC_evo_likelihood_esmc_full_VDJ, method = "pearson")$estimate
  esm1b <- cor.test(df$HC_evo_likelihood_esm1b_full_VDJ, df$LC_evo_likelihood_esm1b_full_VDJ, method = "pearson")$estimate
  protbert <- cor.test(df$HC_evo_likelihood_protbert_full_VDJ, df$LC_evo_likelihood_protbert_full_VDJ, method = "pearson")$estimate
  sapiens <- cor.test(df$HC_evo_likelihood_sapiens_full_VDJ, df$LC_evo_likelihood_sapiens_full_VDJ, method = "pearson")$estimate
  cor_df <- rbind(cor_df,
                   data.frame(sample = sample,
                              ablang1 = ablang1,
                              ablang2 = ablang2,
                              esmc = esmc,
                              esm1b = esm1b,
                              protbert = protbert,
                              sapiens = sapiens))
}
cor_df <- pivot_longer(cor_df, cols = 2:7, names_to = "PLM", values_to = "correlation")
pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figureS4/PLMcorrelation_sup_HCLC.pdf")
ggplot(cor_df, aes(x=PLM, y=correlation, col=factor(sample))) + 
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
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  xlab("PLM") + ylab("Correlation Coefficient") + ggtitle("Heavy Chain vs. Light Chain SP")
dev.off()

#Correlation v-gene distribution OAS and likelihood
oas_human <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/OAS_statistics/sequence_statistics_human.csv")
oas_mouse <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/OAS_statistics/sequence_statistics_mouse.csv")
average_likelihood_ova <- vdj_ova %>% group_by(v_gene_family) %>% 
  summarise(esmc = mean(HC_evo_likelihood_esmc_full_VDJ),
            ablang2 = mean(paired_evo_likelihood_ablang2_full_VDJ),
            ablang1 = mean(HC_evo_likelihood_ablang1_full_VDJ),
            protbert = mean(HC_evo_likelihood_protbert_full_VDJ),
            sapiens = mean(HC_evo_likelihood_sapiens_full_VDJ, na.rm = T),
            esm1b = mean(HC_evo_likelihood_esm1b_full_VDJ))
colnames(average_likelihood_ova) <- c("V_Gene_Family", "ESM-C", "Ablang2", "Ablang1", "Protbert", "Sapiens", "ESM-1b")
df_ova <- inner_join(oas_mouse[,c("V_Gene_Family", "Percentage", "Frequency")], average_likelihood_ova, by = "V_Gene_Family")
df_ova <- pivot_longer(df_ova, cols = c("ESM-C", "Ablang2", "Ablang1", "Protbert", "Sapiens", "ESM-1b"), names_to = "PLM", values_to = "Average_Pseudolikelihood")


pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figureS5/OAS_correlation_mouse.pdf", height = 6, width = 20)
ggplot(df_ova, aes(x = Percentage, y = Average_Pseudolikelihood, color = V_Gene_Family)) +
  geom_point(size = 4) +
  theme_steropodon() +
  ggtitle("Mouse") +
  xlab("Percentage in OAS") + ylab("Average SP") +
  theme(text = element_text(size = 25)) +
  scale_color_manual(name = "V-gene family", values = vgene_colors,
                     breaks = c("IGHV1", "IGHV2", "IGHV3","IGHV4","IGHV5","IGHV6", "IGHV7","IGHV8",
                                "IGHV9", "IGHV10", "IGHV11", "IGHV12", "IGHV13", "IGHV14", "IGHV15")) +
  facet_wrap(~PLM, ncol = 6)
dev.off()


df_ova <- as.data.frame(df_ova)
mean(cor.test(df_ova[df_ova$PLM == "ESM-C", "Percentage"], df_ova[df_ova$PLM == "ESM-C", "Average_Pseudolikelihood"])$estimate,
cor.test(df_ova[df_ova$PLM == "ESM-1b", "Percentage"], df_ova[df_ova$PLM == "ESM-1b", "Average_Pseudolikelihood"])$estimate,
cor.test(df_ova[df_ova$PLM == "Ablang2", "Percentage"], df_ova[df_ova$PLM == "Ablang2", "Average_Pseudolikelihood"])$estimate,
cor.test(df_ova[df_ova$PLM == "Ablang1", "Percentage"], df_ova[df_ova$PLM == "Ablang1", "Average_Pseudolikelihood"])$estimate,
cor.test(df_ova[df_ova$PLM == "Protbert", "Percentage"], df_ova[df_ova$PLM == "Protbert", "Average_Pseudolikelihood"])$estimate,
cor.test(df_ova[df_ova$PLM == "Sapiens", "Percentage"], df_ova[df_ova$PLM == "Sapiens", "Average_Pseudolikelihood"])$estimate)

average_likelihood_human <- vdj_human %>% group_by(v_gene_family) %>% 
  summarise(esmc = mean(HC_evo_likelihood_esmc_full_VDJ),
            ablang2 = mean(paired_evo_likelihood_ablang2_full_VDJ),
            ablang1 = mean(HC_evo_likelihood_ablang1_full_VDJ),
            protbert = mean(HC_evo_likelihood_protbert_full_VDJ),
            sapiens = mean(HC_evo_likelihood_sapiens_full_VDJ, na.rm = T),
            esm1b = mean(HC_evo_likelihood_esm1b_full_VDJ))
colnames(average_likelihood_human) <- c("V_Gene_Family",  "ESM-C", "Ablang2", "Ablang1", "Protbert", "Sapiens", "ESM-1b")
df_human <- inner_join(oas_human[,c("V_Gene_Family", "Percentage", "Frequency")], average_likelihood_human, by = "V_Gene_Family")
df_human <- pivot_longer(df_human, cols = c( "ESM-C", "Ablang2", "Ablang1", "Protbert", "Sapiens", "ESM-1b"), names_to = "PLM", values_to = "Average_Pseudolikelihood")


pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figureS5/OAS_correlation_human.pdf", height = 6, width = 20)
ggplot(df_human, aes(x = Percentage, y = Average_Pseudolikelihood, color = V_Gene_Family)) +
  geom_point(size = 4) +
  theme_steropodon() +
  ggtitle("Human") +
  xlab("Percentage in OAS") + ylab("Average SP") +
  theme(text = element_text(size = 25)) +
  scale_color_manual(name = "V-gene family", values = vgene_colors,
                     breaks = c("IGHV1", "IGHV2", "IGHV3","IGHV4","IGHV5","IGHV6", "IGHV7","IGHV8",
                                "IGHV9", "IGHV10", "IGHV11", "IGHV12", "IGHV13", "IGHV14", "IGHV15")) +
  facet_wrap(~PLM, ncol = 6)
dev.off()

df_human <- as.data.frame(df_human)
mean(cor.test(df_human[df_human$PLM == "ESM-C", "Percentage"], df_human[df_human$PLM == "ESM-C", "Average_Pseudolikelihood"])$estimate,
cor.test(df_human[df_human$PLM == "ESM-1b", "Percentage"], df_human[df_human$PLM == "ESM-1b", "Average_Pseudolikelihood"])$estimate,
cor.test(df_human[df_human$PLM == "Ablang2", "Percentage"], df_human[df_human$PLM == "Ablang2", "Average_Pseudolikelihood"])$estimate,
cor.test(df_human[df_human$PLM == "Ablang1", "Percentage"], df_human[df_human$PLM == "Ablang1", "Average_Pseudolikelihood"])$estimate,
cor.test(df_human[df_human$PLM == "Protbert", "Percentage"], df_human[df_human$PLM == "Protbert", "Average_Pseudolikelihood"])$estimate,
cor.test(df_human[df_human$PLM == "Sapiens", "Percentage"], df_human[df_human$PLM == "Sapiens", "Average_Pseudolikelihood"])$estimate)



