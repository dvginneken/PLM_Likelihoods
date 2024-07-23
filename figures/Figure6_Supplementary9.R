library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)

#Figure 6
##OVA
#ESM
df<-read.delim("PLM_Likelihoods/data/OVA_V7/evo-likelihoods/evo_likelihood_esm.csv",sep=",", header = T)
df$ELISA <- case_match(df$Bind..ELISA.signal.0.2.,
                       "yes" ~ "Binder",
                       "no" ~ "Non-binder")
df$ELISA <- as.factor(df$ELISA)
stat.test <- df %>%
  t_test(evo_likelihood ~ ELISA) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

#Figure 6A (left)
mouse_esm <- ggplot(df, aes(y = evo_likelihood, x = ELISA, color=ELISA)) +
  geom_boxplot() +
  scale_color_manual(values = c("Binder" = "#7a0177",
                                "Non-binder" = "#f768a1")) +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("ESM-1b Pseudolikelihood") +
  ylim(-0.75,0) +
  geom_signif(comparisons=list(c("Binder","Non-binder")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ggtitle("Mouse")

##Kim
#ESM
df<-read.delim("PLM_Likelihoods/data/Kim/evo-likelihoods/evo_likelihood_esm.csv",sep=",", header = T)
df$ELISA <- case_match(df$elisa,
                       "True" ~ "Binder",
                       "False" ~ "Non-binder")
df$ELISA <- as.factor(df$ELISA)
stat.test <- df %>%
  t_test(evo_likelihood ~ ELISA) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

human_esm <- ggplot(df, aes(y = evo_likelihood, x = ELISA, color=ELISA)) +
  geom_boxplot() +
  scale_color_manual(values = c("Binder" = "#005824",
                                "Non-binder" = "#41ae76")) +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("ESM-1b Pseudolikelihood") +
  ylim(-1.2,0) +
  geom_signif(comparisons=list(c("Binder","Non-binder")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ggtitle("Human")

#Plot Figure 6A
pdf("PLM_Likelihoods/figures/Figure6_Supplementary9/Binders_main.pdf")
ggarrange(mouse_esm, human_esm, ncol =2)
dev.off()

#Supplementary Figure 9A
#ProtBERT
df<-read.delim("PLM_Likelihoods/data/OVA_V7/evo-likelihoods/evo_likelihood_protbert.csv",sep=",", header = T)
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
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("ProtBERT Pseudolikelihood") +
  ylim(-0.75,0) +
  geom_signif(comparisons=list(c("Binder","Non-binder")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ggtitle("ProtBERT")

#ablang
df<-read.delim("PLM_Likelihoods/data/OVA_V7/evo-likelihoods/evo_likelihood_ablang.csv",sep=",", header = T)
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
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("Ablang Pseudolikelihood") +
  ylim(-0.75,0) +
  geom_signif(comparisons=list(c("Binder","Non-binder")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ggtitle("Ablang")

#sapiens
df<-read.delim("PLM_Likelihoods/data/OVA_V7/evo-likelihoods/evo_likelihood_sapiens.csv",sep=",", header = T)
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
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("Sapiens Pseudolikelihood") +
  ylim(-0.75,0) +
  geom_signif(comparisons=list(c("Binder","Non-binder")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ggtitle("Sapiens")

#Plot Supplementary Figure 9A
pdf("PLM_Likelihoods/figures/Figure6_Supplementary9/ELISA_mice_sup.pdf")
ggarrange(mouse_protbert, mouse_ablang, mouse_sapiens, ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()

#Supplementary Figure 9C
#Kim
#ProtBERT
df<-read.delim("PLM_Likelihoods/data/Kim/evo-likelihoods/evo_likelihood_protbert.csv",sep=",", header = T)
df$ELISA <- case_match(df$elisa,
                       "True" ~ "Binder",
                       "False" ~ "Non-binder")
df$ELISA <- as.factor(df$ELISA)
stat.test <- df %>%
  t_test(evo_likelihood ~ ELISA) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

human_protbert <- ggplot(df, aes(y = evo_likelihood, x = ELISA, color=ELISA)) +
  geom_boxplot() +
  scale_color_manual(values = c("Binder" = "#005824",
                                "Non-binder" = "#41ae76")) +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("ProtBERT Pseudolikelihood") +
  ylim(-1.2,0) +
  geom_signif(comparisons=list(c("Binder","Non-binder")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ggtitle("ProtBERT")

#ablang
df<-read.delim("PLM_Likelihoods/data/Kim/evo-likelihoods/evo_likelihood_ablang.csv",sep=",", header = T)
df$ELISA <- case_match(df$elisa,
                       "True" ~ "Binder",
                       "False" ~ "Non-binder")
df$ELISA <- as.factor(df$ELISA)
stat.test <- df %>%
  t_test(evo_likelihood ~ ELISA) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

human_ablang <- ggplot(df, aes(y = evo_likelihood, x = ELISA, color=ELISA)) +
  geom_boxplot() +
  scale_color_manual(values = c("Binder" = "#005824",
                                "Non-binder" = "#41ae76")) +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("Ablang Pseudolikelihood") +
  ylim(-1.2,0) +
  geom_signif(comparisons=list(c("Binder","Non-binder")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ggtitle("Ablang")

#ablang
df<-read.delim("PLM_Likelihoods/data/Kim/evo-likelihoods/evo_likelihood_sapiens.csv",sep=",", header = T)
df$ELISA <- case_match(df$elisa,
                       "True" ~ "Binder",
                       "False" ~ "Non-binder")
df$ELISA <- as.factor(df$ELISA)
stat.test <- df %>%
  t_test(evo_likelihood ~ ELISA) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

human_sapiens <- ggplot(df, aes(y = evo_likelihood, x = ELISA, color=ELISA)) +
  geom_boxplot() +
  scale_color_manual(values = c("Binder" = "#005824",
                                "Non-binder" = "#41ae76")) +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("Sapiens Pseudolikelihood") +
  ylim(-1.2,0) +
  geom_signif(comparisons=list(c("Binder","Non-binder")),
              map_signif_level = TRUE,
              annotations=stat.test$p.adj.signif,
              color = "black",
              size = 1,
              textsize = 5) +
  ggtitle("Sapiens")

#Plot Supplementary Figure 9C
pdf("PLM_Likelihoods/figures/Figure6_Supplementary9/ELISA_human_sup.pdf")
ggarrange(human_protbert, human_ablang, human_sapiens, ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()

#Figure 6B
##OVA
#ESM
df<-read.delim("PLM_Likelihoods/data/OVA_V7/evo-likelihoods/evo_likelihood_esm.csv",sep=",", header = T)
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

#Plot Figure 6B
pdf("PLM_Likelihoods/figures/Figure6_Supplementary9/Affinity_polyclonal_OVA_main.pdf")
ggplot(df, aes(x = evo_likelihood, y = as.numeric(octet.affinity..nM.), color = as.factor(Mouse_clone_HC))) +
  geom_point(size=2) +
  scale_color_manual(values = c("Mouse1" = "#fcc5c0",
                                "Mouse2" = "#f768a1",
                                "Mouse5" = "#7a0177"),
                     name = "Sample") +
  theme_minimal() +
  theme(text = element_text(size = 14)) +
  ylab("Affinity (Kd)") +
  xlab("ESM-1b Pseudolikelihood") +
  ggtitle(paste0("Mouse - Polyclonal, R\u00b2 = ", round(cor_mouse, digits = 3)))
dev.off()

#Supplementary Figure 9B
#ProtBERT
df<-read.delim("PLM_Likelihoods/data/OVA_V7/evo-likelihoods/evo_likelihood_protbert.csv",sep=",", header = T)
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
protbert <- ggplot(df, aes(x = evo_likelihood, y = as.numeric(octet.affinity..nM.), color = as.factor(Mouse_clone_HC))) +
  geom_point(size=2) +
  scale_color_manual(values = c("Mouse1" = "#fcc5c0",
                                "Mouse2" = "#f768a1",
                                "Mouse5" = "#7a0177"),
                     name = "Sample") +
  theme_minimal() +
  theme(text = element_text(size = 14)) +
  ylab("Affinity (Kd)") +
  xlab("ProtBERTPseudolikelihood") +
  ggtitle(paste0("ProtBERT, R\u00b2 = ", round(cor_mouse, digits = 3)))

#Ablang
df<-read.delim("PLM_Likelihoods/data/OVA_V7/evo-likelihoods/evo_likelihood_ablang.csv",sep=",", header = T)
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
  theme_minimal() +
  theme(text = element_text(size = 14)) +
  ylab("Affinity (Kd)") +
  xlab("Ablang Pseudolikelihood") +
  ggtitle(paste0("Ablang, R\u00b2 = ", round(cor_mouse, digits = 3)))

#Sapiens
df<-read.delim("PLM_Likelihoods/data/OVA_V7/evo-likelihoods/evo_likelihood_sapiens.csv",sep=",", header = T)
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
  theme_minimal() +
  theme(text = element_text(size = 14)) +
  ylab("Affinity (Kd)") +
  xlab("Sapiens Pseudolikelihood") +
  ggtitle(paste0("Sapiens, R\u00b2 = ", round(cor_mouse, digits = 3)))

#Plot Supplementary Figure 9B
pdf("PLM_Likelihoods/figures/Figure6_Supplementary9/Affinity_polyclonal_OVA_sup.pdf")
ggarrange(protbert, ablang, sapiens, ncol = 3, common.legend = T, legend = "right")
dev.off()

#Figure 6C - Polyclonal Affinity human
df_kim <- read.csv("PLM_Likelihoods/data/Kim/evo-likelihoods/evo_likelihood_esm.csv", header = T)
df_kim$donor <- gsub("368", "Participant", df_kim$donor)
cor <- cor.test(df_kim$evo_likelihood, df_kim$K_D_nM, method = "spearman")$estimate

#Plot Figure 6B
pdf("PLM_Likelihoods/figures/Figure6_Supplementary9/Affinity_polyclonal_human_main.pdf")
ggplot(df_kim, aes(x=evo_likelihood, y=K_D_nM, color=donor)) +
  geom_point(size=2) +
  scale_color_manual(values = c('#ccece6','#99d8c9','#66c2a4','#41ae76','#238b45','#006d2c','#00441b'),
                     name = "Sample") +
  theme_minimal() +
  theme(text = element_text(size = 14)) +
  ylab("Affinity (Kd)") + xlab("ESM-1b pseudolikelihood") +
  ggtitle(paste0("Human - Polyclonal, R\u00b2 = ", round(cor, digits = 3)))
dev.off()

#supplementals
#probert
df_kim <- read.csv("PLM_Likelihoods/data/Kim/evo-likelihoods/evo_likelihood_protbert.csv", header = T)
df_kim$donor <- gsub("368", "Participant", df_kim$donor)
cor <- cor.test(df_kim$evo_likelihood, df_kim$K_D_nM, method = "spearman")$estimate
protbert <- ggplot(df_kim, aes(x=evo_likelihood, y=K_D_nM, color=donor)) +
  geom_point(size=2) +
  scale_color_manual(values = c('#ccece6','#99d8c9','#66c2a4','#41ae76','#238b45','#006d2c','#00441b'),
                     name = "Sample") +
  theme_minimal() +
  theme(text = element_text(size = 14)) +
  ylab("Affinity (Kd)") + xlab("ProtBERT Pseudolikelihood") +
  ggtitle(paste0("ProtBERT, R\u00b2 = ", round(cor, digits = 3)))

#ablang
df_kim <- read.csv("PLM_Likelihoods/data/Kim/evo-likelihoods/evo_likelihood_ablang.csv", header = T)
df_kim$donor <- gsub("368", "Participant", df_kim$donor)
cor <- cor.test(df_kim$evo_likelihood, df_kim$K_D_nM, method = "spearman")$estimate
ablang <- ggplot(df_kim, aes(x=evo_likelihood, y=K_D_nM, color=donor)) +
  geom_point(size=2) +
  scale_color_manual(values = c('#ccece6','#99d8c9','#66c2a4','#41ae76','#238b45','#006d2c','#00441b'),
                     name = "Sample") +
  theme_minimal() +
  theme(text = element_text(size = 14)) +
  ylab("Affinity (Kd)") + xlab("Ablang Ppseudolikelihood") +
  ggtitle(paste0("Ablang, R\u00b2 = ", round(cor, digits = 3)))

#sapiens
df_kim <- read.csv("PLM_Likelihoods/data/Kim/evo-likelihoods/evo_likelihood_sapiens.csv", header = T)
df_kim$donor <- gsub("368", "Participant", df_kim$donor)
cor <- cor.test(df_kim$evo_likelihood, df_kim$K_D_nM, method = "spearman")$estimate
sapiens <- ggplot(df_kim, aes(x=evo_likelihood, y=K_D_nM, color=donor)) +
  geom_point(size=2) +
  scale_color_manual(values = c('#ccece6','#99d8c9','#66c2a4','#41ae76','#238b45','#006d2c','#00441b'),
                     name = "Sample") +
  theme_minimal() +
  theme(text = element_text(size = 14)) +
  ylab("Affinity (Kd)") + xlab("Sapiens Pseudolikelihood") +
  ggtitle(paste0("Sapiens, R\u00b2 = ", round(cor, digits = 3)))

#Plot Supplementary Figure 9D
pdf("PLM_Likelihoods/figures/Figure6_Supplementary9/Affinity_polyclonal_human_sup.pdf")
ggarrange(protbert, ablang, sapiens, ncol = 3, common.legend = T, legend = "right")
dev.off()

#Figure 6 - Variant Tree
load("PLM_Likelihoods/data/OVA_V7/AF_VariantTree.RData") #af

#Source AntibodyForests function from https://github.com/alexyermanos/AntibodyForests
source("~/Documents/GitHub/Platypus/R/AntibodyForests_plot.R")

#Plot Figure 6D
pdf("PLM_Likelihoods/figures/Figure6_Supplementary9/LineageTree_pll.pdf")
AntibodyForests_plot(af,
                     sample = "S1",
                     clonotype = "clonotype1",
                     color.by = "evo_likelihood",
                     label.by = "size",
                     node.label.size = 0.95,
                     show.size.legend = F,
                     color.legend.title = "Pseudolikelihood")
dev.off()

#Plot Figure 6E
pdf("PLM_Likelihoods/figures/Figure6_Supplementary9/LineageTree_affinity.pdf")
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

#Figure 6F - Affinity vs pseudolikelihood with distance to germline
df <- read.csv("PLM_Likelihoods/data/OVA_V7/VariantTree/evo_likelihoods/evo_likelihood_esm.csv", header = T)
cor <- cor.test(df$evo_likelihood, df$octet.affinity, method = "spearman")$estimate
tree = af[["S1"]][["clonotype1"]][["igraph"]]
nodes = igraph::V(tree)[names(igraph::V(tree)) != "germline"]
node_features = af[["S1"]][["clonotype1"]][["nodes"]][names(af[["S1"]][["clonotype1"]][["nodes"]]) != "germline"]
#Get the total length of shortest paths between each node and the germline
distance <- igraph::distances(tree, v = "germline", to = nodes, algorithm = "dijkstra",
                              weights = edge_attr(tree)$edge.length)
distance = t(as.data.frame(distance))
esm = unlist(lapply(node_features,function(x){unique(x[["evo_likelihood"]])[!is.na(unique(x[["evo_likelihood"]]))]}))[rownames(distance)]
barcode = unlist(lapply(node_features,function(x){unique(x[["barcodes"]])[!is.na(unique(x[["barcodes"]]))]}))[rownames(distance)]
distance = data.frame("distance" = distance, "node" = rownames(distance), "barcode" = barcode)
df <- left_join(df, distance, by = "barcode")

#Plot Figure 6F
pdf("PLM_Likelihoods/figures/Figure6_Supplementary9/Figure6F.pdf")
ggplot(df, aes(x=evo_likelihood, y=octet.affinity, color=germline)) +
  geom_point(size=2) +
  scale_colour_gradient(low = "black", high = "orange",
                        name = "Distance to\ngermline") +
  theme_steropodon() +
  theme(text = element_text(size = 14)) +
  ylab("Affinity (Kd)") + xlab("ESM-1b pseudolikelihood") +
  ggtitle(paste0("Mouse - Monoclonal, R\u00b2 = ", round(cor, digits = 3)))
dev.off()