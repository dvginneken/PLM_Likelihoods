SHM_per_alignment <- function(seq_1, seq_2){
  
  hamming_dist <- stringdist(seq_1,seq_2, method = "hamming")
  gaps_1 <- str_count(seq_1,"-")
  gaps_2 <- str_count(seq_2,"-")
  
  
  return(hamming_dist - gaps_1 - gaps_2)
  
} 

SHM_calculator <- function(VDJ_obj){
  
  ref_nt <- gsub("-","",VDJ_obj$VDJ_germline_nt_trimmed)
  alignments <- pairwiseAlignment(VDJ_obj$VDJ_sequence_nt_trimmed, ref_nt)
  
  
  SHM_values <- numeric(length(alignments))
  for (i in 1:length(alignments)){
    
    seq_1 <- as.character(alignments[i]@pattern)
    seq_2 <- as.character(alignments[i]@subject)
    
    if (startsWith(seq_1,"NANA") | startsWith(seq_2,"NANA")){
      SHM_values[i] <- NA
    }
    
    shm <- SHM_per_alignment(seq_1,seq_2)
    SHM_values[i] <- shm
    
  }
  VDJ_obj$SHM_count <- SHM_values
  return(VDJ_obj)
  
}

SHM_calculator_aa <- function(VDJ_obj){
  
  ref_aa <- gsub("-","",VDJ_obj$VDJ_ref.aa)
  alignments <- pairwiseAlignment(VDJ_obj$full_VDJ_aa, ref_aa)
  
  
  SHM_values <- numeric(length(alignments))
  for (i in 1:length(alignments)){
    
    seq_1 <- as.character(alignments[i]@pattern)
    seq_2 <- as.character(alignments[i]@subject)
    
    if (startsWith(seq_1,"NANA") | startsWith(seq_2,"NANA")){
      SHM_values[i] <- NA
    }
    
    shm <- SHM_per_alignment(seq_1,seq_2)
    SHM_values[i] <- shm
    
  }
  VDJ_obj$SHM_count_aa <- SHM_values
  return(VDJ_obj)
  
}

IMGT_SHM_calculator <- function(vdj_obj, dataset = "OVA_mouse", level = "aa"){
  
  IMGT_germlines <- read.table(paste0("../data/",dataset,"/vgene_germline_sequences.csv"),header = TRUE, sep = ",")
  
  vdj_obj$IMGT_SHM_count <- NA
  vdj_obj$v_gene_sequence <- apply(vdj_obj[,c("HC_fwr1","HC_cdr1","HC_fwr2","HC_cdr2","HC_fwr3")], 1, paste, collapse = "")
  
  
  
  for (i in 1:nrow(IMGT_germlines)){
    germline_sequence = IMGT_germlines[i,"VDJ_ref.aa"]
    germline_sequence = gsub('.','', germline_sequence, fixed =TRUE)
    v_gene = IMGT_germlines[i,"v_gene"]
    
    
    alignments <- pairwiseAlignment(vdj_obj[vdj_obj["v_gene"] == v_gene,]$v_gene_sequence, germline_sequence)
    
    
    SHM_values <- numeric(length(alignments))
    
    for (i in 1:length(alignments)){
      
      seq_1 <- as.character(alignments[i]@pattern)
      seq_2 <- as.character(alignments[i]@subject)
      
      if (startsWith(seq_1,"NANA") | startsWith(seq_2,"NANA")){
        SHM_values[i] <- NA
      }
      
      shm <- SHM_per_alignment(seq_1,seq_2)
      SHM_values[i] <- shm
      
    }
    
    vdj_obj[vdj_obj["v_gene"] == v_gene,"IMGT_SHM_count"] <- SHM_values
  }
  return(vdj_obj)
}

distance_calculator <- function(VDJ_obj, distance = "levenshtein"){
  
  ref_nt <- str_replace(VDJ_obj$VDJ_ref.nt,"-","")
  full_VDJ <- VDJ_obj$full_VDJ
  
  distance_values <- stringdist::stringdist(ref_nt,full_VDJ, method = distance)
  
  for (i in 1:length(full_VDJ)){
    
    seq_1 <- full_VDJ[i]
    
    if (startsWith(seq_1,"NANA")){
      distance_values[i] <- NA
    }
    
  }
  
  VDJ_obj[,c(paste0(distance,"_dist"))] <- distance_values
  
  return(VDJ_obj)
  
}