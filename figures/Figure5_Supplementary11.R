library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)
library(AntibodyForests)
library(igraph)

source("~/OneDrive - UMC Utrecht/Documenten/Steropodon_theme.R")
#Figure S11
##OVA
#ESMC
df<-read.delim("~/Documents/GitHub/PLM-likelihoods/data/OVA_V7/Specificity/OVA_elisa_evo_likelihood_esmc.csv",sep=",", header = T)
df$ELISA <- case_match(df$Bind..ELISA.signal.0.2.,
                       "yes" ~ "Binder",
                       "no" ~ "Non-binder")
df$ELISA <- as.factor(df$ELISA)
stat.test <- df %>%
  t_test(evo_likelihood ~ ELISA) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

mouse_esmc <- ggplot(df, aes(y = evo_likelihood, x = ELISA, color=ELISA)) +
  geom_boxplot() +
  scale_color_manual(values = c("Binder" = "#7a0177",
                                "Non-binder" = "#f768a1")) +
  theme_steropodon() +
  theme(text = element_text(size = 14),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("SP") +
  ylim(-1.5,0) +
  geom_signif(comparisons=list(c("Binder","Non-binder")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ggtitle("ESM-C")

#Ablang2
df<-read.delim("~/Documents/GitHub/PLM-likelihoods/data/OVA_V7/Specificity/OVA_elisa_evo_likelihood_ablang2.csv",sep=",", header = T)
df$ELISA <- case_match(df$Bind..ELISA.signal.0.2.,
                       "yes" ~ "Binder",
                       "no" ~ "Non-binder")
df$ELISA <- as.factor(df$ELISA)
stat.test <- df %>%
  t_test(pair_evo_likelihood ~ ELISA) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

mouse_ablang2 <- ggplot(df, aes(y = pair_evo_likelihood, x = ELISA, color=ELISA)) +
  geom_boxplot() +
  scale_color_manual(values = c("Binder" = "#7a0177",
                                "Non-binder" = "#f768a1")) +
  theme_steropodon() +
  theme(text = element_text(size = 14),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("SP") +
  ylim(-1.5,0) +
  geom_signif(comparisons=list(c("Binder","Non-binder")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ggtitle("Ablang2")

#ESM1b
df<-read.delim("~/Documents/GitHub/PLM-likelihoods/data/OVA_V7/Specificity/OVA_elisa_evo_likelihood_esm1b.csv",sep=",", header = T)
df$ELISA <- case_match(df$Bind..ELISA.signal.0.2.,
                       "yes" ~ "Binder",
                       "no" ~ "Non-binder")
df$ELISA <- as.factor(df$ELISA)
stat.test <- df %>%
  t_test(evo_likelihood ~ ELISA) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

mouse_esm1b <- ggplot(df, aes(y = evo_likelihood, x = ELISA, color=ELISA)) +
  geom_boxplot() +
  scale_color_manual(values = c("Binder" = "#7a0177",
                                "Non-binder" = "#f768a1")) +
  theme_steropodon() +
  theme(text = element_text(size = 14),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("SP") +
  ylim(-1.5,0) +
  geom_signif(comparisons=list(c("Binder","Non-binder")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ggtitle("ESM-1b")

#ProtBERT
df<-read.delim("~/Documents/GitHub/PLM-likelihoods/data/OVA_V7/Specificity/OVA_elisa_evo_likelihood_protbert.csv",sep=",", header = T)
df$ELISA <- case_match(df$Bind..ELISA.signal.0.2.,
                       "yes" ~ "Binder",
                       "no" ~ "Non-binder")
df$ELISA <- as.factor(df$ELISA)
stat.test <- df %>%
  t_test(evo_likelihood ~ ELISA) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test
mouse_protbert <- ggplot(df, aes(y = evo_likelihood, x = ELISA, color=ELISA)) +
  geom_boxplot() +
  scale_color_manual(values = c("Binder" = "#7a0177",
                                "Non-binder" = "#f768a1")) +
  theme_steropodon() +
  theme(text = element_text(size = 14),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("SP") +
  ylim(-1.5,0) +
  geom_signif(comparisons=list(c("Binder","Non-binder")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ggtitle("ProtBERT")

#ablang
df<-read.delim("~/Documents/GitHub/PLM-likelihoods/data/OVA_V7/Specificity/OVA_elisa_evo_likelihood_ablang.csv",sep=",", header = T)
df$ELISA <- case_match(df$Bind..ELISA.signal.0.2.,
                       "yes" ~ "Binder",
                       "no" ~ "Non-binder")
df$ELISA <- as.factor(df$ELISA)
stat.test <- df %>%
  t_test(evo_likelihood ~ ELISA) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

mouse_ablang <- ggplot(df, aes(y = evo_likelihood, x = ELISA, color=ELISA)) +
  geom_boxplot() +
  scale_color_manual(values = c("Binder" = "#7a0177",
                                "Non-binder" = "#f768a1")) +
  theme_steropodon() +
  theme(text = element_text(size = 14),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("SP") +
  ylim(-1.5,0) +
  geom_signif(comparisons=list(c("Binder","Non-binder")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ggtitle("Ablang1")

#sapiens
df<-read.delim("~/Documents/GitHub/PLM-likelihoods/data/OVA_V7/Specificity/OVA_elisa_evo_likelihood_sapiens.csv",sep=",", header = T)
df$ELISA <- case_match(df$Bind..ELISA.signal.0.2.,
                       "yes" ~ "Binder",
                       "no" ~ "Non-binder")
df$ELISA <- as.factor(df$ELISA)
stat.test <- df %>%
  t_test(evo_likelihood ~ ELISA) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

mouse_sapiens <- ggplot(df, aes(y = evo_likelihood, x = ELISA, color=ELISA)) +
  geom_boxplot() +
  scale_color_manual(values = c("Binder" = "#7a0177",
                                "Non-binder" = "#f768a1")) +
  theme_steropodon() +
  theme(text = element_text(size = 14),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("SP") +
  ylim(-1.5,0) +
  geom_signif(comparisons=list(c("Binder","Non-binder")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ggtitle("Sapiens")

pdf("~/Documents/GitHub/PLM-likelihoods/figures/figureS11/ELISA_mice_sup.pdf", height= 6, width = 20)
ggarrange(mouse_esmc, mouse_ablang2, mouse_esm1b, mouse_protbert, mouse_ablang, mouse_sapiens, ncol = 6, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()

#Kim
#ESMC
df_esmc<-read.delim("~/Documents/GitHub/PLM-likelihoods/data/Kim/Specificity/Kim_elisa_evo_likelihood_esmc.csv",sep=",", header = T)
df_esmc$ELISA <- case_match(df_esmc$elisa,
                       "True" ~ "Binder",
                       "False" ~ "Non-binder")
df_esmc$ELISA <- as.factor(df_esmc$ELISA)
stat.test <- df_esmc %>%
  t_test(evo_likelihood ~ ELISA) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

human_esmc <- ggplot(df_esmc, aes(y = evo_likelihood, x = ELISA, color=ELISA)) +
  geom_boxplot() +
  scale_color_manual(values = c("Binder" = "#005824",
                                "Non-binder" = "#41ae76")) +
  theme_steropodon() +
  theme(text = element_text(size = 14),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("SP") +
  ylim(-1.5,0) +
  geom_signif(comparisons=list(c("Binder","Non-binder")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ggtitle("ESM-C")

#Ablang2
df_ablang2<-read.delim("~/Documents/GitHub/PLM-likelihoods/data/Kim/Specificity/Kim_elisa_evo_likelihood_ablang2.csv",sep=",", header = T)
df_ablang2$ELISA <- case_match(df_ablang2$elisa,
                       "True" ~ "Binder",
                       "False" ~ "Non-binder")
df_ablang2$ELISA <- as.factor(df_ablang2$ELISA)
stat.test <- df_ablang2 %>%
  t_test(pair_evo_likelihood ~ ELISA) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

human_ablang2 <- ggplot(df_ablang2, aes(y = pair_evo_likelihood, x = ELISA, color=ELISA)) +
  geom_boxplot() +
  scale_color_manual(values = c("Binder" = "#005824",
                                "Non-binder" = "#41ae76")) +
  theme_steropodon() +
  theme(text = element_text(size = 14),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("SP") +
  ylim(-1.5,0) +
  geom_signif(comparisons=list(c("Binder","Non-binder")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ggtitle("Ablang2")

#ESM1b
df_esm1b<-read.delim("~/Documents/GitHub/PLM-likelihoods/data/Kim/Specificity/Kim_elisa_evo_likelihood_esm1b.csv",sep=",", header = T)
df_esm1b$ELISA <- case_match(df_esm1b$elisa,
                       "True" ~ "Binder",
                       "False" ~ "Non-binder")
df_esm1b$ELISA <- as.factor(df_esm1b$ELISA)
stat.test <- df_esm1b %>%
  t_test(evo_likelihood ~ ELISA) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

human_esm1b <- ggplot(df_esm1b, aes(y = evo_likelihood, x = ELISA, color=ELISA)) +
  geom_boxplot() +
  scale_color_manual(values = c("Binder" = "#005824",
                                "Non-binder" = "#41ae76")) +
  theme_steropodon() +
  theme(text = element_text(size = 14),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("SP") +
  ylim(-1.5,0) +
  geom_signif(comparisons=list(c("Binder","Non-binder")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ggtitle("ESM-1b")

#ProtBERT
df_protbert<-read.delim("~/Documents/GitHub/PLM-likelihoods/data/Kim/Specificity/Kim_elisa_evo_likelihood_protbert.csv",sep=",", header = T)
df_protbert$ELISA <- case_match(df_protbert$elisa,
                       "True" ~ "Binder",
                       "False" ~ "Non-binder")
df_protbert$ELISA <- as.factor(df_protbert$ELISA)
stat.test <- df_protbert %>%
  t_test(evo_likelihood ~ ELISA) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test
human_protbert <- ggplot(df_protbert, aes(y = evo_likelihood, x = ELISA, color=ELISA)) +
  geom_boxplot() +
  scale_color_manual(values = c("Binder" = "#005824",
                                "Non-binder" = "#41ae76")) +
  theme_steropodon() +
  theme(text = element_text(size = 14),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("SP") +
  ylim(-1.5,0) +
  geom_signif(comparisons=list(c("Binder","Non-binder")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ggtitle("ProtBERT")

#ablang
df_ablang<-read.delim("~/Documents/GitHub/PLM-likelihoods/data/Kim/Specificity/Kim_elisa_evo_likelihood_ablang.csv",sep=",", header = T)
df_ablang$ELISA <- case_match(df_ablang$elisa,
                       "True" ~ "Binder",
                       "False" ~ "Non-binder")
df_ablang$ELISA <- as.factor(df_ablang$ELISA)
stat.test <- df_ablang %>%
  t_test(evo_likelihood ~ ELISA) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

human_ablang <- ggplot(df_ablang, aes(y = evo_likelihood, x = ELISA, color=ELISA)) +
  geom_boxplot() +
  scale_color_manual(values = c("Binder" = "#005824",
                                "Non-binder" = "#41ae76")) +
  theme_steropodon() +
  theme(text = element_text(size = 14),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("SP") +
  ylim(-1.5,0) +
  geom_signif(comparisons=list(c("Binder","Non-binder")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ggtitle("Ablang1")

#sapiens
df_sapiens<-read.delim("~/Documents/GitHub/PLM-likelihoods/data/Kim/Specificity/Kim_elisa_evo_likelihood_sapiens.csv",sep=",", header = T)
df_sapiens$ELISA <- case_match(df_sapiens$elisa,
                       "True" ~ "Binder",
                       "False" ~ "Non-binder")
df_sapiens$ELISA <- as.factor(df_sapiens$ELISA)
stat.test <- df_sapiens %>%
  t_test(evo_likelihood ~ ELISA) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

human_sapiens <- ggplot(df_sapiens, aes(y = evo_likelihood, x = ELISA, color=ELISA)) +
  geom_boxplot() +
  scale_color_manual(values = c("Binder" = "#005824",
                                "Non-binder" = "#41ae76")) +
  theme_steropodon() +
  theme(text = element_text(size = 14),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("SP") +
  ylim(-1.5,0) +
  geom_signif(comparisons=list(c("Binder","Non-binder")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ggtitle("Sapiens")

#Plot Supplementary Figure 8A
pdf("~/Documents/GitHub/PLM-likelihoods/figures/figureS11/ELISA_human_sup.pdf", height= 6, width = 20)
ggarrange(human_esmc, human_ablang2, human_esm1b, human_protbert, human_ablang, human_sapiens, ncol = 6, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()

#Statistics
human_binders <- c(df_ablang[df_ablang$isa == "True",]$evo_likelihood,
                   df_esmc[df_esmc$elisa == "True",]$evo_likelihood,
                   df_ablang2[df_ablang2$elisa == "True",]$pair_evo_likelihood,
                   df_esm1b[df_esm1b$elisa == "True",]$evo_likelihood,
                   df_protbert[df_protbert$elisa == "True",]$evo_likelihood,
                   df_sapiens[df_sapiens$elisa == "True",]$evo_likelihood)
human_nonbinder <- c(df_ablang[df_ablang$isa == "False",]$evo_likelihood,
                   df_esmc[df_esmc$elisa == "False",]$evo_likelihood,
                   df_ablang2[df_ablang2$elisa == "False",]$pair_evo_likelihood,
                   df_esm1b[df_esm1b$elisa == "False",]$evo_likelihood,
                   df_protbert[df_protbert$elisa == "False",]$evo_likelihood,
                   df_sapiens[df_sapiens$elisa == "False",]$evo_likelihood)

mean(human_binders)
mean(human_nonbinder)

#Figure 5A
##OVA
#ESMC
df<-read.delim("~/Documents/GitHub/PLM-likelihoods/data/OVA_V7/Specificity/OVA_elisa_evo_likelihood_esmc.csv",sep=",", header = T)
df<-df[df$octet.affinity..nM. != "",]
df<-df[df$octet.affinity..nM. != "nd",]
df$octet.affinity..nM. <- gsub(",",".",df$octet.affinity..nM.)
df$Mouse_clone_HC <- case_match(df$Mouse_clone_HC,
                                1 ~ "Mouse1",
                                2 ~ "Mouse2",
                                3 ~ "Mouse3",
                                4 ~ "Mouse4",
                                5 ~ "Mouse5")

cor_mouse <- cor.test(df$evo_likelihood, as.numeric(df$octet.affinity..nM.), method = "spearman")$estimate

#Plot Figure 5A
pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figure5/Affinity_polyclonal_OVA_esmc.pdf")
ggplot(df, aes(x = evo_likelihood, y = as.numeric(octet.affinity..nM.), color = as.factor(Mouse_clone_HC))) +
  geom_point(size=2) +
  scale_color_manual(values = c("Mouse1" = "#fcc5c0",
                                "Mouse2" = "#f768a1",
                                "Mouse5" = "#7a0177"),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  ylab("Affinity (Kd)") +
  xlab("ESM-C SP") +
  ggtitle(paste0("Mouse - Polyclonal\nR\u00b2 = ", round(cor_mouse, digits = 3)))
dev.off()

#Ablang2
df<-read.delim("~/Documents/GitHub/PLM-likelihoods/data/OVA_V7/Specificity/OVA_elisa_evo_likelihood_ablang2.csv",sep=",", header = T)
df<-df[df$octet.affinity..nM. != "",]
df<-df[df$octet.affinity..nM. != "nd",]
df$octet.affinity..nM. <- gsub(",",".",df$octet.affinity..nM.)
df$Mouse_clone_HC <- case_match(df$Mouse_clone_HC,
                                1 ~ "Mouse1",
                                2 ~ "Mouse2",
                                3 ~ "Mouse3",
                                4 ~ "Mouse4",
                                5 ~ "Mouse5")

cor_mouse <- cor.test(df$pair_evo_likelihood, as.numeric(df$octet.affinity..nM.), method = "spearman")$estimate

#Plot Figure 5A
pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figure5/Affinity_polyclonal_OVA_ablang2.pdf")
ggplot(df, aes(x = pair_evo_likelihood, y = as.numeric(octet.affinity..nM.), color = as.factor(Mouse_clone_HC))) +
  geom_point(size=2) +
  scale_color_manual(values = c("Mouse1" = "#fcc5c0",
                                "Mouse2" = "#f768a1",
                                "Mouse5" = "#7a0177"),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  ylab("Affinity (Kd)") +
  xlab("Ablang2 SP") +
  ggtitle(paste0("Mouse - Polyclonal\nR\u00b2 = ", round(cor_mouse, digits = 3)))
dev.off()


#Supplementary Figure 11B
#ProtBERT
df<-read.delim("~/Documents/GitHub/PLM-likelihoods/data/OVA_V7/Specificity/OVA_elisa_evo_likelihood_protbert.csv",sep=",", header = T)
df<-df[df$octet.affinity..nM. != "",]
df<-df[df$octet.affinity..nM. != "nd",]
df$octet.affinity..nM. <- gsub(",",".",df$octet.affinity..nM.)
df$Mouse_clone_HC <- case_match(df$Mouse_clone_HC,
                                1 ~ "Mouse1",
                                2 ~ "Mouse2",
                                3 ~ "Mouse3",
                                4 ~ "Mouse4",
                                5 ~ "Mouse5")

cor_mouse <- cor.test(df$evo_likelihood, as.numeric(df$octet.affinity..nM.), method = "spearman")$estimate
protbert <- ggplot(df, aes(x = evo_likelihood, y = as.numeric(octet.affinity..nM.), color = as.factor(Mouse_clone_HC))) +
  geom_point(size=2) +
  scale_color_manual(values = c("Mouse1" = "#fcc5c0",
                                "Mouse2" = "#f768a1",
                                "Mouse5" = "#7a0177"),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  ylab("Affinity (Kd)") +
  xlab("SP") +
  ggtitle(paste0("ProtBERT\nR\u00b2 = ", round(cor_mouse, digits = 3)))

#Ablang
df<-read.delim("~/Documents/GitHub/PLM-likelihoods/data/OVA_V7/Specificity/OVA_elisa_evo_likelihood_ablang.csv",sep=",", header = T)
df<-df[df$Bind..ELISA.signal.0.2. == "yes",]
df<-df[df$octet.affinity..nM. != "",]
df<-df[df$octet.affinity..nM. != "nd",]
df$octet.affinity..nM. <- gsub(",",".",df$octet.affinity..nM.)
df$Mouse_clone_HC <- case_match(df$Mouse_clone_HC,
                                1 ~ "Mouse1",
                                2 ~ "Mouse2",
                                3 ~ "Mouse3",
                                4 ~ "Mouse4",
                                5 ~ "Mouse5")

cor_mouse <- cor.test(df$evo_likelihood, as.numeric(df$octet.affinity..nM.), method = "spearman")$estimate
ablang <- ggplot(df, aes(x = evo_likelihood, y = as.numeric(octet.affinity..nM.), color = as.factor(Mouse_clone_HC))) +
  geom_point(size=2) +
  scale_color_manual(values = c("Mouse1" = "#fcc5c0",
                                "Mouse2" = "#f768a1",
                                "Mouse5" = "#7a0177"),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  ylab("Affinity (Kd)") +
  xlab("SP") +
  ggtitle(paste0("Ablang1\nR\u00b2 = ", round(cor_mouse, digits = 3)))

#Sapiens
df<-read.delim("~/Documents/GitHub/PLM-likelihoods/data/OVA_V7/Specificity/OVA_elisa_evo_likelihood_sapiens.csv",sep=",", header = T)
df<-df[df$Bind..ELISA.signal.0.2. == "yes",]
df<-df[df$octet.affinity..nM. != "",]
df<-df[df$octet.affinity..nM. != "nd",]
df$octet.affinity..nM. <- gsub(",",".",df$octet.affinity..nM.)
df$Mouse_clone_HC <- case_match(df$Mouse_clone_HC,
                                1 ~ "Mouse1",
                                2 ~ "Mouse2",
                                3 ~ "Mouse3",
                                4 ~ "Mouse4",
                                5 ~ "Mouse5")

cor_mouse <- cor.test(df$evo_likelihood, as.numeric(df$octet.affinity..nM.), method = "spearman")$estimate
sapiens <- ggplot(df, aes(x = evo_likelihood, y = as.numeric(octet.affinity..nM.), color = as.factor(Mouse_clone_HC))) +
  geom_point(size=2) +
  scale_color_manual(values = c("Mouse1" = "#fcc5c0",
                                "Mouse2" = "#f768a1",
                                "Mouse5" = "#7a0177"),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  ylab("Affinity (Kd)") +
  xlab("SP") +
  ggtitle(paste0("Sapiens\nR\u00b2 = ", round(cor_mouse, digits = 3)))

#esm1b
df<-read.delim("~/Documents/GitHub/PLM-likelihoods/data/OVA_V7/Specificity/OVA_elisa_evo_likelihood_esm1b.csv",sep=",", header = T)
df<-df[df$octet.affinity..nM. != "",]
df<-df[df$octet.affinity..nM. != "nd",]
df$octet.affinity..nM. <- gsub(",",".",df$octet.affinity..nM.)
df$Mouse_clone_HC <- case_match(df$Mouse_clone_HC,
                                1 ~ "Mouse1",
                                2 ~ "Mouse2",
                                3 ~ "Mouse3",
                                4 ~ "Mouse4",
                                5 ~ "Mouse5")
cor_mouse <- cor.test(df$evo_likelihood, as.numeric(df$octet.affinity..nM.), method = "spearman")$estimate
esm1b <- ggplot(df, aes(x = evo_likelihood, y = as.numeric(octet.affinity..nM.), color = as.factor(Mouse_clone_HC))) +
  geom_point(size=2) +
  scale_color_manual(values = c("Mouse1" = "#fcc5c0",
                                "Mouse2" = "#f768a1",
                                "Mouse5" = "#7a0177"),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  ylab("Affinity (Kd)") +
  xlab("SP") +
  ggtitle(paste0("ESM-1b\nR\u00b2 = ", round(cor_mouse, digits = 3)))

#Plot Supplementary Figure 11
pdf("~/Documents/GitHub/PLM-likelihoods/figures/figureS11/Affinity_ELISA_ova_sup.pdf", height= 6, width = 20)
ggarrange(esm1b, protbert, ablang, sapiens, ncol = 4, common.legend = T, legend = "right")
dev.off()

#Figure S11D - Polyclonal Affinity human
#ESMC
df_kim <- read.csv("~/Documents/GitHub/PLM-likelihoods/data/Kim/Affinity/Kim_affinity_evo_likelihood_esmc.csv", header = T)
df_kim$donor <- gsub("368-", "", df_kim$donor)
cor <- cor.test(df_kim$evo_likelihood, df_kim$K_D_nM, method = "spearman")$estimate

#Plot Figure 5B
pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figure5/Affinity_polyclonal_human_esmc.pdf")
ggplot(df_kim, aes(x=evo_likelihood, y=K_D_nM, color=donor)) +
  geom_point(size=2) +
  scale_color_manual(values = c('#ccece6','#99d8c9','#66c2a4','#41ae76','#238b45','#006d2c','#00441b'),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  ylab("Affinity (Kd)") + xlab("ESM-C SP") +
  ggtitle(paste0("Human - Polyclonal\nR\u00b2 = ", round(cor, digits = 3)))
dev.off()

#Ablang2
df_kim <- read.csv("~/Documents/GitHub/PLM-likelihoods/data/Kim/Affinity/Kim_affinity_evo_likelihood_ablang2.csv", header = T)
df_kim$donor <- gsub("368-", "", df_kim$donor)
cor <- cor.test(df_kim$pair_evo_likelihood, df_kim$K_D_nM, method = "spearman")$estimate

#Plot Figure 5B
pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figure5/Affinity_polyclonal_human_ablang2.pdf")
ggplot(df_kim, aes(x=pair_evo_likelihood, y=K_D_nM, color=donor)) +
  geom_point(size=2) +
  scale_color_manual(values = c('#ccece6','#99d8c9','#66c2a4','#41ae76','#238b45','#006d2c','#00441b'),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  ylab("Affinity (Kd)") + xlab("Ablang2 SP") +
  ggtitle(paste0("Human - Polyclonal\nR\u00b2 = ", round(cor, digits = 3)))
dev.off()

#supplementals
#probert
df_kim <- read.csv("~/Documents/GitHub/PLM-likelihoods/data/Kim/Affinity/Kim_affinity_evo_likelihood_protbert.csv", header = T)
df_kim$donor <- gsub("368", "Participant", df_kim$donor)
cor <- cor.test(df_kim$evo_likelihood, df_kim$K_D_nM, method = "spearman")$estimate
protbert <- ggplot(df_kim, aes(x=evo_likelihood, y=K_D_nM, color=donor)) +
  geom_point(size=2) +
  scale_color_manual(values = c('#ccece6','#99d8c9','#66c2a4','#41ae76','#238b45','#006d2c','#00441b'),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  ylab("Affinity (Kd)") + xlab("ProtBERT Pseudolikelihood") +
  ggtitle(paste0("ProtBERT\nR\u00b2 = ", round(cor, digits = 3)))

#ablang
df_kim <- read.csv("~/Documents/GitHub/PLM-likelihoods/data/Kim/Affinity/Kim_affinity_evo_likelihood_ablang.csv", header = T)
df_kim$donor <- gsub("368", "Participant", df_kim$donor)
cor <- cor.test(df_kim$evo_likelihood, df_kim$K_D_nM, method = "spearman")$estimate
ablang <- ggplot(df_kim, aes(x=evo_likelihood, y=K_D_nM, color=donor)) +
  geom_point(size=2) +
  scale_color_manual(values = c('#ccece6','#99d8c9','#66c2a4','#41ae76','#238b45','#006d2c','#00441b'),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  ylab("Affinity (Kd)") + xlab("Ablang Ppseudolikelihood") +
  ggtitle(paste0("Ablang\nR\u00b2 = ", round(cor, digits = 3)))

#sapiens
df_kim <- read.csv("~/Documents/GitHub/PLM-likelihoods/data/Kim/Affinity/Kim_affinity_evo_likelihood_sapiens.csv", header = T)
df_kim$donor <- gsub("368", "Participant", df_kim$donor)
cor <- cor.test(df_kim$evo_likelihood, df_kim$K_D_nM, method = "spearman")$estimate
sapiens <- ggplot(df_kim, aes(x=evo_likelihood, y=K_D_nM, color=donor)) +
  geom_point(size=2) +
  scale_color_manual(values = c('#ccece6','#99d8c9','#66c2a4','#41ae76','#238b45','#006d2c','#00441b'),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  ylab("Affinity (Kd)") + xlab("Sapiens Pseudolikelihood") +
  ggtitle(paste0("Sapiens\nR\u00b2 = ", round(cor, digits = 3)))

#esm1b
df_kim <- read.csv("~/Documents/GitHub/PLM-likelihoods/data/Kim/Affinity/Kim_affinity_evo_likelihood_esm1b.csv", header = T)
df_kim$donor <- gsub("368", "Participant", df_kim$donor)
cor <- cor.test(df_kim$evo_likelihood, df_kim$K_D_nM, method = "spearman")$estimate
esm1b <- ggplot(df_kim, aes(x=evo_likelihood, y=K_D_nM, color=donor)) +
  geom_point(size=2) +
  scale_color_manual(values = c('#ccece6','#99d8c9','#66c2a4','#41ae76','#238b45','#006d2c','#00441b'),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  ylab("Affinity (Kd)") + xlab("ESM-1b Pseudolikelihood") +
  ggtitle(paste0("ESM-1b\nR\u00b2 = ", round(cor, digits = 3)))

#Plot Supplementary Figure 11D
pdf("~/Documents/GitHub/PLM-likelihoods/figures/figureS11/Affinity_ELISA_human_sup.pdf", height= 6, width = 20)
ggarrange(esm1b, protbert, ablang, sapiens, ncol = 4, common.legend = T, legend = "right")
dev.off()

#Figure 5 - Variant Tree
load("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/OVA_V7/Affinity/AF_default_HCLC.RData")

#Plot Figure 5D
#ESMC
affinity_esmc <- read_csv("~/Documents/GitHub/PLM-likelihoods/data/OVA_V7/Affinity/OVA_affinity_evo_likelihood_esmc.csv")
af <- Af_add_node_feature(af, affinity_esmc, c("evo_likelihood", "octet.affinity"))
Af_plot_tree(af,
                     sample = "S1",
                     clonotype = "clonotype1",
                     color.by = "evo_likelihood",
                     label.by = "size",
                     node.label.size = 0.95,
                     show.size.legend = F,
                     color.legend.title = "ESM-C SP",
             output.file = "/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figure5/VariantTree_esmc.pdf")

#Ablang2
affinity_ablang2 <- read_csv("~/Documents/GitHub/PLM-likelihoods/data/OVA_V7/Affinity/OVA_affinity_evo_likelihood_ablang2.csv")
af <- Af_add_node_feature(af, affinity_ablang2, c("pair_evo_likelihood", "HC_evo_likelihood"))
Af_plot_tree(af,
             sample = "S1",
             clonotype = "clonotype1",
             color.by = "pair_evo_likelihood",
             label.by = "size",
             node.label.size = 0.95,
             show.size.legend = F,
             color.legend.title = "Ablang2 SP",
             output.file = "/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figure5/VariantTree_ablang2.pdf")


#Plot Figure 5E
pdf("PLM_Likelihoods/figures/Figure5_Supplementary8/LineageTree_affinity.pdf")
AntibodyForests_plot(af,
                     sample = "S1",
                     clonotype = "clonotype1",
                     color.by = "octet.affinity",
                     node.color.gradient = c("white", "red"),
                     label.by = "size",
                     node.label.size = 0.95,
                     show.size.legend = F,
                     color.legend.title = "Affinity (Kd)")
dev.off()

#Figure 5F - Affinity vs pseudolikelihood with distance to germline
#ESMC
df <- read.csv("~/Documents/GitHub/PLM-likelihoods/data/OVA_V7/Affinity/OVA_affinity_evo_likelihood_esmc.csv", header = T)
cor <- cor.test(df$evo_likelihood, df$octet.affinity, method = "spearman")$estimate
tree = af[["S1"]][["clonotype1"]][["igraph"]]
nodes = igraph::V(tree)[names(igraph::V(tree)) != "germline"]
node_features = af[["S1"]][["clonotype1"]][["nodes"]][names(af[["S1"]][["clonotype1"]][["nodes"]]) != "germline"]
#Get the total length of shortest paths between each node and the germline
distance <- igraph::distances(tree, v = "germline", to = nodes, algorithm = "dijkstra",
                              weights = edge_attr(tree)$edge.length)
distance = t(as.data.frame(distance))
esmc = unlist(lapply(node_features,function(x){unique(x[["evo_likelihood"]])[!is.na(unique(x[["evo_likelihood"]]))]}))[rownames(distance)]
barcode = unlist(lapply(node_features,function(x){unique(x[["barcodes"]])[!is.na(unique(x[["barcodes"]]))]}))[rownames(distance)]
distance = data.frame("distance" = distance, "node" = rownames(distance), "barcode" = barcode)
df <- left_join(df, distance, by = "barcode")

#Plot Figure 5F
pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figure5/VariantTree_correlation_esmc.pdf")
ggplot(df, aes(x=evo_likelihood, y=octet.affinity, color=germline)) +
  geom_point(size=2) +
  scale_colour_gradient(low = "black", high = "orange",
                        name = "Distance to\ngermline") +
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  ylab("Affinity (Kd)") + xlab("ESM-C SP") +
  ggtitle(paste0("Mouse - Monoclonal\nR\u00b2 = ", round(cor, digits = 3)))
dev.off()

#Ablang2
df <- read.csv("~/Documents/GitHub/PLM-likelihoods/data/OVA_V7/Affinity/OVA_affinity_evo_likelihood_ablang2.csv", header = T)
cor <- cor.test(df$pair_evo_likelihood, df$octet.affinity, method = "spearman")$estimate
tree = af[["S1"]][["clonotype1"]][["igraph"]]
nodes = igraph::V(tree)[names(igraph::V(tree)) != "germline"]
node_features = af[["S1"]][["clonotype1"]][["nodes"]][names(af[["S1"]][["clonotype1"]][["nodes"]]) != "germline"]
#Get the total length of shortest paths between each node and the germline
distance <- igraph::distances(tree, v = "germline", to = nodes, algorithm = "dijkstra",
                              weights = edge_attr(tree)$edge.length)
distance = t(as.data.frame(distance))
ablang2 = unlist(lapply(node_features,function(x){unique(x[["pair_evo_likelihood"]])[!is.na(unique(x[["pair_evo_likelihood"]]))]}))[rownames(distance)]
barcode = unlist(lapply(node_features,function(x){unique(x[["barcodes"]])[!is.na(unique(x[["barcodes"]]))]}))[rownames(distance)]
distance = data.frame("distance" = distance, "node" = rownames(distance), "barcode" = barcode)
df <- left_join(df, distance, by = "barcode")

#Plot Figure 5F
pdf("/Users/dginneke/Documents/GitHub/PLM-likelihoods/figures/figure5/VariantTree_correlation_ablang2.pdf")
ggplot(df, aes(x=pair_evo_likelihood, y=octet.affinity, color=germline)) +
  geom_point(size=2) +
  scale_colour_gradient(low = "black", high = "orange",
                        name = "Distance to\ngermline") +
  theme_steropodon() +
  theme(text = element_text(size = 25)) +
  ylab("Affinity (Kd)") + xlab("Ablang2 SP") +
  ggtitle(paste0("Mouse - Monoclonal\nR\u00b2 = ", round(cor, digits = 3)))
dev.off()
