library(igraph)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)

#Supplementary Figure 7
#Load data
load("PLM_Likelihoods/data/OVA_V7/AF_OVA_V7_mp_HC.RData") #af
af_ova <- af
rm(af)
load("PLM_Likelihoods/data/Horns/AF_horns2020a__VDJ_RAW_mp_HC.RData")
af_horns <- af
rm(af)
load("PLM_Likelihoods/data/Bruhn/AF_Bruhn_mp_HC.RData")
af_bruhn <- af
rm(af)
load("PLM_Likelihoods/data/Kim/AF_Kim_mp_HC.RData")
af_kim <- af
rm(af)

#A - Lineage Trees
#Source AntibodyForests function from https://github.com/alexyermanos/AntibodyForests
source("~/Documents/GitHub/Platypus/R/AntibodyForests_plot.R")

#Plot Supplementary Figure 7 A (left)
pdf("PLM_Likelihoods/figures/Figure4_Supplementary67/LineageTree_OVA.pdf")
AntibodyForests_plot(af_ova,
                     sample = "S1",
                     clonotype = "clonotype4",
                     color.by = "IGH_evo_likelihood_esm_full_VDJ",
                     label.by = "size",
                     show.color.legend = F,
                     show.size.legend = F,
                     main.title = "Mouse 1",
                     sub.title = "Clonotype 4",
                     arrow.size = 0.5,
                     show.inner.nodes = T)
dev.off()

#Plot Supplementary Figure 7 A (right)
pdf("PLM_Likelihoods/figures/Figure4_Supplementary67/LineageTree_Kim.pdf")
AntibodyForests_plot(af_kim,
                     sample = "SRR17729703",
                     clonotype = "clonotype1",
                     color.by = "IGH_evo_likelihood_esm_full_VDJ",
                     label.by = "size",
                     show.color.legend = F,
                     show.size.legend = F,
                     main.title = "Individual 3",
                     sub.title = "Clonotype 1",
                     arrow.size = 0.5,
                     show.inner.nodes = T)
dev.off()

#B - Distance to germline
#Function to get a dataframe with the distance to the germline for each node and the pseudolikelihoods
distance_germline <- function(input){
  #Go over each tree in the AntibodyForests object and create a metric dataframe
  metric_df <- lapply(seq_along(input), function(sample){
    lapply(seq_along(input[[sample]]), function(clonotype){
      sample_name <- names(input)[sample]
      clonotype_name <- names(input[[sample_name]])[clonotype]
      
      if (length(input[[sample_name]][[clonotype_name]][["nodes"]]) < 3){
        tree = input[[sample_name]][[clonotype_name]][["igraph"]]
      }else{
        tree = input[[sample_name]][[clonotype_name]][["igraph.with.inner.nodes"]]}
      nodes = igraph::V(tree)[names(igraph::V(tree)) != "germline"]
      node_features = input[[sample_name]][[clonotype_name]][["nodes"]][names(input[[sample_name]][[clonotype_name]][["nodes"]]) != "germline"]
      
      #Get the total length of shortest paths between each node and the germline
      distance <- igraph::distances(tree, v = "germline", to = nodes, algorithm = "dijkstra",
                                    weights = edge_attr(tree)$edge.length)
      distance = t(as.data.frame(distance))
      
      
      #Evo likelihood
      esm = unlist(lapply(node_features,function(x){unique(x[["IGH_evo_likelihood_esm_full_VDJ"]])[!is.na(unique(x[["IGH_evo_likelihood_esm_full_VDJ"]]))]}))[rownames(distance)]
      protbert = unlist(lapply(node_features,function(x){unique(x[["IGH_evo_likelihood_protbert_full_VDJ"]])[!is.na(unique(x[["IGH_evo_likelihood_protbert_full_VDJ"]]))]}))[rownames(distance)]
      ablang = unlist(lapply(node_features,function(x){unique(x[["IGH_evo_likelihood_ablang_full_VDJ"]])[!is.na(unique(x[["IGH_evo_likelihood_ablang_full_VDJ"]]))]}))[rownames(distance)]
      sapiens = unlist(lapply(node_features,function(x){unique(x[["IGH_evo_likelihood_sapiens_full_VDJ"]])[!is.na(unique(x[["IGH_evo_likelihood_sapiens_full_VDJ"]]))]}))[rownames(distance)]
      
      
      distance = data.frame("distance" = distance, "sample" = sample_name, "clonotype" = clonotype_name, 
                            "node" = rownames(distance), "esm" = esm, "protbert" = protbert, 
                            "ablang" = ablang, "sapiens" = sapiens)
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
distance_ova <- distance_germline(af_ova)
distance_horns <- distance_germline(af_horns)
distance_bruhn <- distance_germline(af_bruhn)
distance_kim <- distance_germline(af_kim)

#change sample names
distance_ova$sample <- case_match(distance_ova$sample,
                                  "S1" ~ "Mouse1",
                                  "S2" ~ "Mouse2",
                                  "S3" ~ "Mouse3",
                                  "S4" ~ "Mouse4",
                                  "S5" ~ "Mouse5")

distance_horns$sample <- case_match(distance_horns$sample,
                                    "Influenza.vac.11.12.human.S1" ~ "Individual1",
                                    "Influenza.vac.11.12.human.S2" ~ "Human2",
                                    "Influenza.vac.11.12.human.S3" ~ "Human3",
                                    "Influenza.vac.11.12.human.S4" ~ "Human4")
#Only keep 1 replicate
distance_horns <- distance_horns[distance_horns$sample == "Individual1",]

distance_bruhn$sample <- "Individual2"

distance_kim$sample <- case_match(distance_kim$sample,
                                  "SRR17729703" ~ "Individual3",
                                  "SRR17729692" ~ "Individual4",
                                  "SRR17729726" ~ "Individual5")
#Only keep specific samples
distance_kim <- distance_kim[distance_kim$sample %in% paste0("Individual",3:5),]

distance_all <- rbind(distance_ova, distance_horns, distance_bruhn, distance_kim)
distance_human <- rbind(distance_horns, distance_bruhn, distance_kim)

#Plot the distance to the germline against the pseudolikelihoods
#esm
cor_human <- cor.test(distance_human$germline, distance_human$esm)$estimate
cor_mouse <- cor.test(distance_ova$germline, distance_ova$esm)$estimate
d <- ggplot(distance_all, aes(x = germline, y = esm, colour = sample)) + 
  geom_point(size = 1) + 
  geom_smooth(method = "lm", se=F) +
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
  theme(text = element_text(size = 12)) +
  xlab("Distance to germline") + ylab("ESM-1b Pseudolikelihood") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), ", \nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#protbert
cor_human <- cor.test(distance_human$germline, distance_human$protbert)$estimate
cor_mouse <- cor.test(distance_ova$germline, distance_ova$protbert)$estimate
a <- ggplot(distance_all, aes(x = germline, y = protbert, colour = sample)) + 
  geom_point(size = 1) + 
  geom_smooth(method = "lm", se=F) +
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
  theme(text = element_text(size = 12)) +
  xlab("Distance to germline") + ylab("ProtBERT Pseudolikelihood") +
  ggtitle(paste0("ProtBERT\nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", \nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#sapiens
cor_human <- cor.test(distance_human$germline, distance_human$sapiens)$estimate
cor_mouse <- cor.test(distance_ova$germline, distance_ova$sapiens)$estimate
c <- ggplot(distance_all, aes(x = germline, y = sapiens, colour = sample)) + 
  geom_point(size = 1) + 
  geom_smooth(method = "lm", se=F) +
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
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  xlab("Distance to germline") + ylab("Sapiens Pseudolikelihood") +
  ggtitle(paste0("Sapiens\nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", \nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#ablang
cor_human <- cor.test(distance_human$germline, distance_human$ablang)$estimate
cor_mouse <- cor.test(distance_ova$germline, distance_ova$ablang)$estimate
b <- ggplot(distance_all, aes(x = germline, y = ablang, colour = sample)) + 
  geom_point(size = 1) + 
  geom_smooth(method = "lm", se=F) +
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
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  xlab("Distance to germline") + ylab("Ablang Pseudolikelihood") +
  ggtitle(paste0("Ablang\nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", \nR\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#Plot Supplementary Figure 7B
pdf("PLM_Likelihoods/figures/Figure4_Supplementary67/DistanceGermline_sup.pdf")
ggarrange(d, a, b, c, ncol = 2, nrow = 2, common.legend = TRUE, legend = "right")
dev.off()

