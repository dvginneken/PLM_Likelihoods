library(ggplot2) 
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(ggpubr)

#Figure 4A
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

#normalize by sample size
vdj_ova %>% group_by(sample_id) %>% summarise(sample_size = n()) -> sample_size
vdj_ova$samplesize <- sample_size$sample_size[match(vdj_ova$sample_id, sample_size$sample_id)]
vdj_ova %>% mutate(normalized_expansion = clonotype_frequency / samplesize) -> vdj_ova
vdj_human %>% group_by(sample_id) %>% summarise(sample_size = n()) -> sample_size
vdj_human$samplesize <- sample_size$sample_size[match(vdj_human$sample_id, sample_size$sample_id)]
vdj_human %>% mutate(normalized_expansion = clonotype_frequency / samplesize) -> vdj_human
vdj_all <- rbind(vdj_ova, vdj_human)

#esm
cor_human <- cor.test(vdj_human$normalized_expansion, vdj_human$IGH_evo_likelihood_esm_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$normalized_expansion, vdj_ova$IGH_evo_likelihood_esm_full_VDJ)$estimate
a <- ggplot(vdj_all, aes(x=normalized_expansion, y=IGH_evo_likelihood_esm_full_VDJ, color=sample_id)) +
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
  scale_x_continuous(trans='log10') +
  geom_smooth(method = "lm", se=F) +
  theme_minimal() +
  theme(text = element_text(size = 14)) +
  xlab("Clonal Expansion") + ylab("ESM-1b Pseudolikelihood") +
  ggtitle(paste0("ESM-1b\nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#protbert
cor_human <- cor.test(vdj_human$normalized_expansion, vdj_human$IGH_evo_likelihood_protbert_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$normalized_expansion, vdj_ova$IGH_evo_likelihood_protbert_full_VDJ)$estimate
b <- ggplot(vdj_all, aes(x=normalized_expansion, y=IGH_evo_likelihood_protbert_full_VDJ, color=sample_id)) +
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
  scale_x_continuous(trans='log10') +
  geom_smooth(method = "lm", se=F) +
  theme_minimal() +
  theme(text = element_text(size = 14)) +
  xlab("Clonal Expansion") + ylab("ProtBERT Pseudolikelihood") +
  ggtitle(paste0("ProtBERT\nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#sapiens
cor_human <- cor.test(vdj_human$normalized_expansion, vdj_human$IGH_evo_likelihood_sapiens_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$normalized_expansion, vdj_ova$IGH_evo_likelihood_sapiens_full_VDJ)$estimate
c <- ggplot(vdj_all, aes(x=normalized_expansion, y=IGH_evo_likelihood_sapiens_full_VDJ, color=sample_id)) +
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
  scale_x_continuous(trans='log10') +
  geom_smooth(method = "lm", se=F) +
  theme_minimal() +
  theme(text = element_text(size = 14)) +
  xlab("Clonal Expansion") + ylab("Sapiens Pseudolikelihood") +
  ggtitle(paste0("Sapiens\nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#ablang
cor_human <- cor.test(vdj_human$normalized_expansion, vdj_human$IGH_evo_likelihood_ablang_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$normalized_expansion, vdj_ova$IGH_evo_likelihood_ablang_full_VDJ)$estimate
d <- ggplot(vdj_all, aes(x=normalized_expansion, y=IGH_evo_likelihood_ablang_full_VDJ, color=sample_id)) +
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
  scale_x_continuous(trans='log10') +
  geom_smooth(method = "lm", se=F) +
  theme_minimal() +
  theme(text = element_text(size = 14)) +
  xlab("Clonal Expansion") + ylab("Ablang Pseudolikelihood") +
  ggtitle(paste0("Ablang\nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#Plot figure4A
pdf("PLM_Likelihoods/figures/Figure4AB_Supplementary4/Expansion_main.pdf", width = 8, height = 6)
ggarrange(a, d, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()

#Plot supplementary figure4A
pdf("PLM_Likelihoods/figures/Figure4AB_Supplementary4/Expansion_sup.pdf", width = 8, height = 6)
ggarrange(b, c, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()

#Figure 4B
#esm
cor_human <- cor.test(vdj_human$SHM_count, vdj_human$IGH_evo_likelihood_esm_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$SHM_count, vdj_ova$IGH_evo_likelihood_esm_full_VDJ)$estimate
a <- ggplot(vdj_all, aes(x=SHM_count, y=IGH_evo_likelihood_esm_full_VDJ, color=sample_id)) +
  geom_point(size =1) +
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
  geom_smooth(method = "lm", se=F) +
  theme_minimal() +
  theme(text = element_text(size = 14)) +
  xlab("SHM count") + ylab("ESM-1b Pseudolikelihood") +
  ggtitle(paste0("ESM-1b\nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#protbert
cor_human <- cor.test(vdj_human$SHM_count, vdj_human$IGH_evo_likelihood_protbert_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$SHM_count, vdj_ova$IGH_evo_likelihood_protbert_full_VDJ)$estimate
b <- ggplot(vdj_all, aes(x=SHM_count, y=IGH_evo_likelihood_protbert_full_VDJ, color=sample_id)) +
  geom_point(size =1) +
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
  geom_smooth(method = "lm", se=F) +
  theme_minimal() +
  theme(text = element_text(size = 14)) +
  xlab("SHM count") + ylab("ProtBERT Pseudolikelihood") +
  ggtitle(paste0("ProtBERT\nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#Sapiens
cor_human <- cor.test(vdj_human$SHM_count, vdj_human$IGH_evo_likelihood_sapiens_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$SHM_count, vdj_ova$IGH_evo_likelihood_sapiens_full_VDJ)$estimate
c <- ggplot(vdj_all, aes(x=SHM_count, y=IGH_evo_likelihood_sapiens_full_VDJ, color=sample_id)) +
  geom_point(size =1) +
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
  geom_smooth(method = "lm", se=F) +
  theme_minimal() +
  theme(text = element_text(size = 14)) +
  xlab("SHM count") + ylab("Sapiens Pseudolikelihood") +
  ggtitle(paste0("Sapiens\nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#ablang
cor_human <- cor.test(vdj_human$SHM_count, vdj_human$IGH_evo_likelihood_ablang_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$SHM_count, vdj_ova$IGH_evo_likelihood_ablang_full_VDJ)$estimate
d <- ggplot(vdj_all, aes(x=SHM_count, y=IGH_evo_likelihood_ablang_full_VDJ, color=sample_id)) +
  geom_point(size =1) +
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
  geom_smooth(method = "lm", se=F) +
  theme_minimal() +
  theme(text = element_text(size = 14)) +
  xlab("SHM count") + ylab("Ablang Pseudolikelihood") +
  ggtitle(paste0("Ablang\nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#Plot figure4B
pdf("PLM_Likelihoods/figures/Figure4AB_Supplementary4/SHM_main.pdf", width = 8, height = 6)
ggarrange(a, d, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()

#Plot supplementary figure4B
pdf("PLM_Likelihoods/figures/Figure4AB_Supplementary4/SHM_sup.pdf", width = 8, height = 6)
ggarrange(b, c, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()
