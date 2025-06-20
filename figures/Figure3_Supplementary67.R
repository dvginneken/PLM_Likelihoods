library(ggplot2)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(tidyr)

source("~/OneDrive - UMC Utrecht/Documenten/Steropodon_theme.R")
setwd("/Users/dginneke/Documents/GitHub")
#Figure 3A
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

#Set colors
isotype_colors <- c("IGHA" = "#fb6a4a",
                    "IGHE" = "#B452CD",
                    "IGHG" ="#74c476",
                    "IGHM" = "black",
                    "IGHD" = "#1874CD",
                    "NA" = "grey")

#Set isotypes
vdj_ova$isotype <- case_match(vdj_ova$isotype,
                              "IgG1" ~ "IGHG",
                              "IgG2" ~ "IGHG",
                              "IgG3" ~ "IGHG",
                              "IGHA" ~ "IGHA",
                              "IGHD" ~ "IGHD",
                              "IGHE" ~ "IGHE",
                              "IGHM" ~ "IGHM")
vdj_ova %>% arrange(factor(isotype, levels = c('IGHM', 'IGHG', 'IGHA', 'IGHE', 'IGHD'))) -> vdj_ova

vdj_human$isotype <- case_match(vdj_human$isotype,
                                "IgG1" ~ "IGHG",
                                "IgG2" ~ "IGHG",
                                "IgG3" ~ "IGHG",
                                "IgG4" ~ "IGHG",
                                "IGHA" ~ "IGHA",
                                "IGHD" ~ "IGHD",
                                "IGHE" ~ "IGHE",
                                "IGHM" ~ "IGHM",
                                "IgA1" ~ "IGHA",
                                "IgA2" ~ "IGHA")
vdj_human %>% arrange(factor(isotype, levels = c('IGHM', 'IGHG', 'IGHA', 'IGHE', 'IGHD'))) -> vdj_human

vdj_ova <- vdj_ova[vdj_ova$isotype %in% c("IGHG", "IGHA", "IGHM"),]
vdj_human <- vdj_human[vdj_human$isotype %in% c("IGHG", "IGHA", "IGHM"),]

vdj_all <- rbind(vdj_ova, vdj_human)

ova_esm <- ggplot(vdj_ova, aes(x=isotype, y=HC_evo_likelihood_esmc_full_VDJ, fill = isotype)) +
  geom_violin() +
  scale_fill_manual(values = isotype_colors) +
  ggsignif::geom_signif(comparisons=list(c("IGHA", "IGHG"), c("IGHG","IGHM"), c("IGHA", "IGHM")), 
                        step_increase = 0.1, test = "t.test", map_signif_level = T) +
  theme_steropodon() +
  theme(legend.position = "none", text = element_text(size = 18)) +
  xlab("Isotype") + ylab("ESM-C SP") + ggtitle("Mice samples")


human_esm <- ggplot(vdj_human, aes(x=isotype, y=HC_evo_likelihood_esmc_full_VDJ, fill = isotype)) +
  geom_violin() +
  scale_fill_manual(values = isotype_colors) +
  ggsignif::geom_signif(comparisons=list(c("IGHA", "IGHG"), c("IGHG","IGHM"), c("IGHA", "IGHM")), 
                        step_increase = 0.1, test = "t.test", map_signif_level = T) +
  theme_steropodon() +
  theme(legend.position = "none", text = element_text(size = 18)) +
  xlab("Isotype") + ylab("ESM-C SP") + ggtitle("Human samples")



ova_ablang <- ggplot(vdj_ova, aes(x=isotype, y=paired_evo_likelihood_ablang2_full_VDJ, fill = isotype)) +
  geom_violin() +
  scale_fill_manual(values = isotype_colors) +
  ggsignif::geom_signif(comparisons=list(c("IGHA", "IGHG"), c("IGHG","IGHM"), c("IGHA", "IGHM")), 
                        step_increase = 0.1, test = "t.test", map_signif_level = T) +  
  theme_steropodon() + 
  theme(legend.position = "none", text = element_text(size = 18)) +
  xlab("Isotype") + ylab("Ablang2 SP") + ggtitle("Mice samples")

human_ablang <- ggplot(vdj_human, aes(x=isotype, y=paired_evo_likelihood_ablang2_full_VDJ, fill = isotype)) +
  geom_violin() +
  scale_fill_manual(values = isotype_colors) +
  ggsignif::geom_signif(comparisons=list(c("IGHA", "IGHG"), c("IGHG","IGHM"), c("IGHA", "IGHM")), 
                        step_increase = 0.1, test = "t.test", map_signif_level = T) +
  theme_steropodon() +
  theme(legend.position = "none", text = element_text(size = 18)) +
  xlab("Isotype") + ylab("Ablang2 SP") + ggtitle("Human samples")

pdf("PLM-likelihoods/figures/figure3/Figure3A.pdf", width = 15, height = 5)
ggarrange(ova_esm, human_esm, ova_ablang, human_ablang, ncol = 4)
dev.off()

#Correlation isotype distribution OAS and likelihood
oas_human <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/OAS_statistics/isotype_statistics_human.csv")
oas_mouse <- read.csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/OAS_statistics/isotype_statistics_mouse.csv")

average_likelihood_ova <- vdj_ova %>% group_by(isotype) %>% 
  summarise(esmc = mean(HC_evo_likelihood_esmc_full_VDJ),
            ablang2 = mean(paired_evo_likelihood_ablang2_full_VDJ),
            ablang1 = mean(HC_evo_likelihood_ablang1_full_VDJ),
            protbert = mean(HC_evo_likelihood_protbert_full_VDJ),
            sapiens = mean(HC_evo_likelihood_sapiens_full_VDJ, na.rm = T),
            esm1b = mean(HC_evo_likelihood_esm1b_full_VDJ))
colnames(average_likelihood_ova) <- c("Isotype", "ESM-C", "Ablang2", "Ablang1", "Protbert", "Sapiens", "ESM-1b")
df_ova <- inner_join(oas_mouse[,c("Isotype", "Percentage_Unique_Sequences")], average_likelihood_ova, by = "Isotype")
df_ova <- pivot_longer(df_ova, cols = c("ESM-C", "Ablang2", "Ablang1", "Protbert", "Sapiens", "ESM-1b"), names_to = "PLM", values_to = "Average_Pseudolikelihood")


pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figureS6/OAS_correlation_mouse.pdf", height = 6, width = 20)
ggplot(df_ova, aes(x = Percentage_Unique_Sequences, y = Average_Pseudolikelihood, color = Isotype)) +
  geom_point(size = 4) +
  theme_steropodon() +
  ggtitle("Mouse") +
  xlab("Percentage in OAS") + ylab("Average SP") +
  theme(text = element_text(size = 25)) +
  scale_color_manual(values = isotype_colors) +
  facet_wrap(~PLM, ncol = 6)
dev.off()

df_ova <- as.data.frame(df_ova)
cor.test(df_ova[df_ova$PLM == "ESM-C", "Percentage_Unique_Sequences"], df_ova[df_ova$PLM == "ESM-C", "Average_Pseudolikelihood"])$estimate
cor.test(df_ova[df_ova$PLM == "ESM-1b", "Percentage_Unique_Sequences"], df_ova[df_ova$PLM == "ESM-1b", "Average_Pseudolikelihood"])$estimate
cor.test(df_ova[df_ova$PLM == "Ablang2", "Percentage_Unique_Sequences"], df_ova[df_ova$PLM == "Ablang2", "Average_Pseudolikelihood"])$estimate
cor.test(df_ova[df_ova$PLM == "Ablang1", "Percentage_Unique_Sequences"], df_ova[df_ova$PLM == "Ablang1", "Average_Pseudolikelihood"])$estimate
cor.test(df_ova[df_ova$PLM == "Protbert", "Percentage_Unique_Sequences"], df_ova[df_ova$PLM == "Protbert", "Average_Pseudolikelihood"])$estimate
cor.test(df_ova[df_ova$PLM == "Sapiens", "Percentage_Unique_Sequences"], df_ova[df_ova$PLM == "Sapiens", "Average_Pseudolikelihood"])$estimate


average_likelihood_human <- vdj_human %>% group_by(isotype) %>% 
  summarise(esmc = mean(HC_evo_likelihood_esmc_full_VDJ),
            ablang2 = mean(paired_evo_likelihood_ablang2_full_VDJ),
            ablang1 = mean(HC_evo_likelihood_ablang1_full_VDJ),
            protbert = mean(HC_evo_likelihood_protbert_full_VDJ),
            sapiens = mean(HC_evo_likelihood_sapiens_full_VDJ, na.rm = T),
            esm1b = mean(HC_evo_likelihood_esm1b_full_VDJ))
colnames(average_likelihood_human) <- c("Isotype", "ESM-C", "Ablang2", "Ablang1", "Protbert", "Sapiens", "ESM-1b")
df_human <- inner_join(oas_human[,c("Isotype", "Percentage_Unique_Sequences")], average_likelihood_human, by = "Isotype")
df_human <- pivot_longer(df_human, cols = c("ESM-C", "Ablang2", "Ablang1", "Protbert", "Sapiens", "ESM-1b"), names_to = "PLM", values_to = "Average_Pseudolikelihood")


pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figureS6/OAS_correlation_human.pdf", height = 6, width = 20)
ggplot(df_human, aes(x = Percentage_Unique_Sequences, y = Average_Pseudolikelihood, color = Isotype)) +
  geom_point(size = 4) +
  theme_steropodon() +
  ggtitle("Human") +
  xlab("Percentage in OAS") + ylab("Average SP") +
  theme(text = element_text(size = 25)) +
  scale_color_manual(values = isotype_colors) +
  facet_wrap(~PLM, ncol = 6)
dev.off()

df_human <- as.data.frame(df_human)
cor.test(df_human[df_human$PLM == "ESM-C", "Percentage_Unique_Sequences"], df_human[df_human$PLM == "ESM-C", "Average_Pseudolikelihood"])$estimate
cor.test(df_human[df_human$PLM == "ESM-1b", "Percentage_Unique_Sequences"], df_human[df_human$PLM == "ESM-1b", "Average_Pseudolikelihood"])$estimate
cor.test(df_human[df_human$PLM == "Ablang2", "Percentage_Unique_Sequences"], df_human[df_human$PLM == "Ablang2", "Average_Pseudolikelihood"])$estimate
cor.test(df_human[df_human$PLM == "Ablang1", "Percentage_Unique_Sequences"], df_human[df_human$PLM == "Ablang1", "Average_Pseudolikelihood"])$estimate
cor.test(df_human[df_human$PLM == "Protbert", "Percentage_Unique_Sequences"], df_human[df_human$PLM == "Protbert", "Average_Pseudolikelihood"])$estimate
cor.test(df_human[df_human$PLM == "Sapiens", "Percentage_Unique_Sequences"], df_human[df_human$PLM == "Sapiens", "Average_Pseudolikelihood"])$estimate

average_likelihood_all <- vdj_all %>% group_by(isotype) %>% 
  summarise(esmc = mean(HC_evo_likelihood_esmc_full_VDJ),
            ablang2 = mean(paired_evo_likelihood_ablang2_full_VDJ),
            ablang1 = mean(HC_evo_likelihood_ablang1_full_VDJ),
            protbert = mean(HC_evo_likelihood_protbert_full_VDJ),
            sapiens = mean(HC_evo_likelihood_sapiens_full_VDJ, na.rm = T),
            esm1b = mean(HC_evo_likelihood_esm1b_full_VDJ))
colnames(average_likelihood_all) <- c("Isotype", "ESM-C", "Ablang2", "Ablang1", "Protbert", "Sapiens", "ESM-1b")

#statistics
#mean IgM
mean(as.numeric(average_likelihood_all[3,2:7]))
#mean IgA and IgG
mean(as.numeric(unlist(average_likelihood_all[1:2,2:7])))

#mean IgM
mean(as.numeric(average_likelihood_ova[3,c(2,3,5)]))
#mean IgA and IgG
mean(as.numeric(unlist(average_likelihood_ova[1:2,c(2,3,5)])))


#Supplementary Figure 5
ova_protbert <- ggplot(vdj_ova, aes(x=isotype, y=HC_evo_likelihood_protbert_full_VDJ, fill = isotype)) +
  geom_violin() +
  scale_fill_manual(values = isotype_colors) +
  ggsignif::geom_signif(comparisons=list(c("IGHA", "IGHG"), c("IGHG","IGHM"), c("IGHA", "IGHM")), 
                        step_increase = 0.1, test = "t.test", map_signif_level = T) +
  theme_minimal() +
  theme(legend.position = "none", text = element_text(size = 18)) +
  xlab("Isotype") + ylab("ProtBERT SP") 

human_protbert <- ggplot(vdj_human, aes(x=isotype, y=HC_evo_likelihood_protbert_full_VDJ, fill = isotype)) +
  geom_violin() +
  scale_fill_manual(values = isotype_colors) +
  ggsignif::geom_signif(comparisons=list(c("IGHA", "IGHG"), c("IGHG","IGHM"), c("IGHA", "IGHM")), 
                        step_increase = 0.1, test = "t.test", map_signif_level = T) +
  theme_minimal() + 
  theme(legend.position = "none", text = element_text(size = 18)) +
  xlab("Isotype") + ylab("ProtBERT SP") 


ova_sapiens <- ggplot(vdj_ova, aes(x=isotype, y=HC_evo_likelihood_sapiens_full_VDJ, fill = isotype)) +
  geom_violin() +
  scale_fill_manual(values = isotype_colors) +
  ggsignif::geom_signif(comparisons=list(c("IGHA", "IGHG"), c("IGHG","IGHM"), c("IGHA", "IGHM")), 
                        step_increase = 0.1, test = "t.test", map_signif_level = T) +  
  theme_minimal() + 
  theme(legend.position = "none", text = element_text(size = 18)) +
  xlab("Isotype") + ylab("Sapiens SP") 

human_sapiens <- ggplot(vdj_human, aes(x=isotype, y=HC_evo_likelihood_sapiens_full_VDJ, fill = isotype)) +
  geom_violin() +
  scale_fill_manual(values = isotype_colors) +
  ggsignif::geom_signif(comparisons=list(c("IGHA", "IGHG"), c("IGHG","IGHM"), c("IGHA", "IGHM")), 
                        step_increase = 0.1, test = "t.test", map_signif_level = T) +
  theme_minimal() +
  theme(legend.position = "none", text = element_text(size = 18)) +
  xlab("Isotype") + ylab("Sapiens SP") 


ova_esm1b <- ggplot(vdj_ova, aes(x=isotype, y=HC_evo_likelihood_esm1b_full_VDJ, fill = isotype)) +
  geom_violin() +
  scale_fill_manual(values = isotype_colors) +
  ggsignif::geom_signif(comparisons=list(c("IGHA", "IGHG"), c("IGHG","IGHM"), c("IGHA", "IGHM")), 
                        step_increase = 0.1, test = "t.test", map_signif_level = T) +  
  theme_minimal() + 
  theme(legend.position = "none", text = element_text(size = 18)) +
  xlab("Isotype") + ylab("ESM-1b SP") 

human_esm1b <- ggplot(vdj_human, aes(x=isotype, y=HC_evo_likelihood_esm1b_full_VDJ, fill = isotype)) +
  geom_violin() +
  scale_fill_manual(values = isotype_colors) +
  ggsignif::geom_signif(comparisons=list(c("IGHA", "IGHG"), c("IGHG","IGHM"), c("IGHA", "IGHM")), 
                        step_increase = 0.1, test = "t.test", map_signif_level = T) +
  theme_minimal() + 
  theme(legend.position = "none", text = element_text(size = 18)) +
  xlab("Isotype") + ylab("ESM-1b SP")

ova_ablang1 <- ggplot(vdj_ova, aes(x=isotype, y=HC_evo_likelihood_ablang1_full_VDJ, fill = isotype)) +
  geom_violin() +
  scale_fill_manual(values = isotype_colors) +
  ggsignif::geom_signif(comparisons=list(c("IGHA", "IGHG"), c("IGHG","IGHM"), c("IGHA", "IGHM")), 
                        step_increase = 0.1, test = "t.test", map_signif_level = T) +  
  theme_minimal() + 
  theme(legend.position = "none", text = element_text(size = 18)) +
  xlab("Isotype") + ylab("Ablang1 SP")

human_ablang1 <- ggplot(vdj_human, aes(x=isotype, y=HC_evo_likelihood_ablang1_full_VDJ, fill = isotype)) +
  geom_violin() +
  scale_fill_manual(values = isotype_colors) +
  ggsignif::geom_signif(comparisons=list(c("IGHA", "IGHG"), c("IGHG","IGHM"), c("IGHA", "IGHM")), 
                        step_increase = 0.1, test = "t.test", map_signif_level = T) +
  theme_minimal() + 
  theme(legend.position = "none", text = element_text(size = 18)) +
  xlab("Isotype") + ylab("Ablang1 SP")

pdf("PLM-likelihoods/figures/figureS6/Isotype_HC_Human.pdf", width = 15, height = 5)
ggarrange(human_ablang1, human_sapiens, human_esm1b, human_protbert, ncol = 4)
dev.off()

pdf("PLM-likelihoods/figures/figureS6/Isotype_HC_Mouse.pdf", width = 15, height = 5)
ggarrange(ova_ablang1, ova_sapiens, ova_esm1b, ova_protbert, ncol = 4)
dev.off()


#Clonal Expansion
#normalize by sample size
vdj_ova %>% group_by(sample_id) %>% summarise(sample_size = n()) -> sample_size
vdj_ova$samplesize <- sample_size$sample_size[match(vdj_ova$sample_id, sample_size$sample_id)]
vdj_ova %>% mutate(normalized_expansion = clonotype_frequency / samplesize) -> vdj_ova
vdj_human %>% group_by(sample_id) %>% summarise(sample_size = n()) -> sample_size
vdj_human$samplesize <- sample_size$sample_size[match(vdj_human$sample_id, sample_size$sample_id)]
vdj_human %>% mutate(normalized_expansion = clonotype_frequency / samplesize) -> vdj_human
vdj_all <- rbind(vdj_ova, vdj_human)

#esm
cor_human <- cor.test(vdj_human$normalized_expansion, vdj_human$HC_evo_likelihood_esmc_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$normalized_expansion, vdj_ova$HC_evo_likelihood_esmc_full_VDJ)$estimate
a <- ggplot(vdj_all, aes(x=normalized_expansion, y=HC_evo_likelihood_esmc_full_VDJ, color=sample_id)) +
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
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  xlab("Clonal Expansion") + ylab("ESM-C SP") +
  ggtitle(paste0("ESM-C\nR\u00b2 Individuals = ", round(cor_human, digits = 2), "\nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#ablang
cor_human <- cor.test(vdj_human$normalized_expansion, vdj_human$paired_evo_likelihood_ablang2_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$normalized_expansion, vdj_ova$paired_evo_likelihood_ablang2_full_VDJ)$estimate
d <- ggplot(vdj_all, aes(x=normalized_expansion, y=paired_evo_likelihood_ablang2_full_VDJ, color=sample_id)) +
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
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  xlab("Clonal Expansion") + ylab("Ablang2 SP") +
  ggtitle(paste0("Ablang2\nR\u00b2 Individuals = ", round(cor_human, digits = 2), "\nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))


pdf("PLM-likelihoods/figures/figure3/Figure3B.pdf", width = 12, height = 6)
ggarrange(a, d, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()

#SHM
cor_human <- cor.test(vdj_human$SHM_count, vdj_human$HC_evo_likelihood_esmc_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$SHM_count, vdj_ova$HC_evo_likelihood_esmc_full_VDJ)$estimate
a <- ggplot(vdj_all, aes(x=SHM_count, y=HC_evo_likelihood_esmc_full_VDJ, color=sample_id)) +
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
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  xlab("SHM count") + ylab("ESM-C SP") +
  ggtitle(paste0("ESM-C\nR\u00b2 Individuals = ", round(cor_human, digits = 2), "\nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#ablang
cor_human <- cor.test(vdj_human$SHM_count, vdj_human$paired_evo_likelihood_ablang2_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$SHM_count, vdj_ova$paired_evo_likelihood_ablang2_full_VDJ)$estimate
d <- ggplot(vdj_all, aes(x=SHM_count, y=paired_evo_likelihood_ablang2_full_VDJ, color=sample_id)) +
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
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  xlab("SHM count") + ylab("Ablang2 SP") +
  ggtitle(paste0("Ablang2\nR\u00b2 Individuals = ", round(cor_human, digits = 2), "\nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#Plot figure4B
pdf("PLM-likelihoods/figures/figure3/Figure3C.pdf", width = 12, height = 6)
ggarrange(a, d, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()

#Supplement 5
#Heavy chain
#protbert
cor_human <- cor.test(vdj_human$normalized_expansion, vdj_human$HC_evo_likelihood_protbert_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$normalized_expansion, vdj_ova$HC_evo_likelihood_protbert_full_VDJ)$estimate
a <- ggplot(vdj_all, aes(x=normalized_expansion, y=HC_evo_likelihood_protbert_full_VDJ, color=sample_id)) +
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
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  xlab("Clonal Expansion") + ylab("ProtBERT SP") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), "\nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#sapiens
cor_human <- cor.test(vdj_human$normalized_expansion, vdj_human$HC_evo_likelihood_sapiens_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$normalized_expansion, vdj_ova$HC_evo_likelihood_sapiens_full_VDJ)$estimate
d <- ggplot(vdj_all, aes(x=normalized_expansion, y=HC_evo_likelihood_sapiens_full_VDJ, color=sample_id)) +
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
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  xlab("Clonal Expansion") + ylab("Sapiens SP") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), "\nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#ESM1b
cor_human <- cor.test(vdj_human$normalized_expansion, vdj_human$HC_evo_likelihood_esm1b_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$normalized_expansion, vdj_ova$HC_evo_likelihood_esm1b_full_VDJ)$estimate
e <- ggplot(vdj_all, aes(x=normalized_expansion, y=HC_evo_likelihood_esm1b_full_VDJ, color=sample_id)) +
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
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  xlab("Clonal Expansion") + ylab("ESM-1b SP") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), "\nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#Ablang1
cor_human <- cor.test(vdj_human$normalized_expansion, vdj_human$HC_evo_likelihood_ablang1_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$normalized_expansion, vdj_ova$HC_evo_likelihood_ablang1_full_VDJ)$estimate
f <- ggplot(vdj_all, aes(x=normalized_expansion, y=HC_evo_likelihood_ablang1_full_VDJ, color=sample_id)) +
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
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  xlab("Clonal Expansion") + ylab("Ablang1 SP") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), "\nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))


#Plot figure4A
pdf("PLM-likelihoods/figures/figureS5/ClonalExpansion_HC.pdf", width = 18, height = 5)
ggarrange(e, f, a, d, ncol = 4, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()


#Figure S5
#protbert
cor_human <- cor.test(vdj_human$SHM_count, vdj_human$HC_evo_likelihood_protbert_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$SHM_count, vdj_ova$HC_evo_likelihood_protbert_full_VDJ)$estimate
a <- ggplot(vdj_all, aes(x=SHM_count, y=HC_evo_likelihood_protbert_full_VDJ, color=sample_id)) +
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
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  xlab("SHM count") + ylab("ProtBERT SP") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), "\nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#sapiens
cor_human <- cor.test(vdj_human$SHM_count, vdj_human$HC_evo_likelihood_sapiens_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$SHM_count, vdj_ova$HC_evo_likelihood_sapiens_full_VDJ)$estimate
d <- ggplot(vdj_all, aes(x=SHM_count, y=HC_evo_likelihood_sapiens_full_VDJ, color=sample_id)) +
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
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  xlab("SHM count") + ylab("Sapiens SP") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), "\nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#ESM1b
cor_human <- cor.test(vdj_human$SHM_count, vdj_human$HC_evo_likelihood_esm1b_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$SHM_count, vdj_ova$HC_evo_likelihood_esm1b_full_VDJ)$estimate
e <- ggplot(vdj_all, aes(x=SHM_count, y=HC_evo_likelihood_esm1b_full_VDJ, color=sample_id)) +
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
  geom_smooth(method = "lm", se=F) +
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  xlab("SHM count") + ylab("ESM-1b SP") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), "\nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#Ablang1
cor_human <- cor.test(vdj_human$SHM_count, vdj_human$HC_evo_likelihood_ablang1_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$SHM_count, vdj_ova$HC_evo_likelihood_ablang1_full_VDJ)$estimate
f <- ggplot(vdj_all, aes(x=SHM_count, y=HC_evo_likelihood_ablang1_full_VDJ, color=sample_id)) +
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
  geom_smooth(method = "lm", se=F) +
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  xlab("SHM count") + ylab("Ablang1 SP") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), "\nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))


pdf("PLM-likelihoods/figures/figureS5/SHM_HC.pdf", width = 18, height = 5)
ggarrange(e, f, a, d, ncol = 4, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()
