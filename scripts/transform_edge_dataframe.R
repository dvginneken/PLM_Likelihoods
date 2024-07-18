#Transform the antibodyforest object into a csv with edge information
library(dplyr)
library(igraph)

args = commandArgs(trailingOnly=TRUE)
dataset <- args[1]
method <- args[2]

load(paste0("../data/",dataset,"/AF_",dataset,"_",method,"_HC.RData")) #af

df_per_clone <- function(sample, clonotype){
  #Igraph object
  tree <- af[[sample]][[clonotype]][["igraph"]]

  #Get edgelist
  edges <- as_edgelist(tree, names = T)
  edges <- as.data.frame(edges)
  
  #Remove germline from the edge list
  edges <- edges[edges$V1 != "germline" & edges$V1 != "germline",]
  colnames(edges) <- c("node1", "node2")
  
  #If there are not enough edges, return NA
  if (nrow(edges) > 0){
    #Get the sequences
    nodes <- af[[sample]][[clonotype]][["nodes"]]
    seq_1 <- c()
    seq_2 <- c()
    for (row in 1:nrow(edges)){
      node1 <- edges[row, "node1"]
      node2 <- edges[row, "node2"]
      seq1 <- nodes[[node1]][["VDJ_sequence_aa_trimmed"]]
      seq2 <- nodes[[node2]][["VDJ_sequence_aa_trimmed"]]
      seq_1 <- c(seq_1, seq1)
      seq_2 <- c(seq_2, seq2)
    }
    df <- cbind(edges, seq_1, seq_2, sample, clonotype)
  }else{
    df <- data.frame(node1 = NA, node2 = NA, seq_1 = NA, seq_2 = NA, sample = sample, clonotype = clonotype)
  }
  return(df)
}

output_df <- data.frame()
for (sample in names(af)){
  for (clonotype in names(af[[sample]])){
    df <- df_per_clone(sample, clonotype)
    output_df <- rbind(output_df, df)
  }
}

output_df <- na.omit(output_df)
write.csv(output_df, file = paste0("../data/",dataset,"/edges_",dataset,"_",method,"_HC.csv"), row.names = F)

