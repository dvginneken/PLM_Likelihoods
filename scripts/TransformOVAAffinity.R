#Transform affinity measurements from the variant tree from Neumeier et al. (https://www.pnas.org/doi/full/10.1073/pnas.2113766119) 
# into a VDJ dataframe

df<-read.delim("../data/OVA_V7/Variant_tree.csv", sep = ",", header = T)

#Extract germline
germline <- df[18,]
df <- df[-18,]

#Add germline sequences to the dataframe
df$VDJ_germline_nt <- germline$HC.nt.seq
df$VJ_germline_nt <- germline$LC.nt.seq

#Translate sequences
translate_DNA<- function(sequence){
  
  #Translate a nucleotide sequence into an amino acid sequence
  #Arguments:
  #- sequence: nucleotide sequence to be translated
  
  if (sequence == ""){
    return(NA)
  }
  
  #Genetic code
  genetic_code <- list(
    "TTT"="F", "TTC"="F", "TTA"="L", "TTG"="L",
    "TCT"="S", "TCC"="S", "TCA"="S", "TCG"="S",
    "TAT"="Y", "TAC"="Y", "TAA"="*", "TAG"="*",
    "TGT"="C", "TGC"="C", "TGA"="*", "TGG"="W",
    "CTT"="L", "CTC"="L", "CTA"="L", "CTG"="L",
    "CCT"="P", "CCC"="P", "CCA"="P", "CCG"="P",
    "CAT"="H", "CAC"="H", "CAA"="Q", "CAG"="Q",
    "CGT"="R", "CGC"="R", "CGA"="R", "CGG"="R",
    "ATT"="I", "ATC"="I", "ATA"="I", "ATG"="M",
    "ACT"="T", "ACC"="T", "ACA"="T", "ACG"="T",
    "AAT"="N", "AAC"="N", "AAA"="K", "AAG"="K",
    "AGT"="S", "AGC"="S", "AGA"="R", "AGG"="R",
    "GTT"="V", "GTC"="V", "GTA"="V", "GTG"="V",
    "GCT"="A", "GCC"="A", "GCA"="A", "GCG"="A",
    "GAT"="D", "GAC"="D", "GAA"="E", "GAG"="E",
    "GGT"="G", "GGC"="G", "GGA"="G", "GGG"="G"
  )
  #Split the sequence into codons
  codons <- strsplit(sequence, "(?<=.{3})", perl=TRUE)[[1]]
  
  #Translate the codons
  for (codon_id in 1:length(codons)){
    #Remove codons that are not complete
    if(nchar(codons[codon_id]) < 3){
      codons[codon_id] = ""
    }
    #Codons that contain "-" are replaced with "-"
    else if (grepl("-", codons[codon_id], fixed = TRUE)){
      codons[codon_id] = "-"
    }
    #Codons that contain "N" are replaced with "-"
    else if (grepl("N", codons[codon_id], fixed = TRUE)){
      codons[codon_id] = "-"
    }
    #Translate codons according to the genetic code
    else{
      codons[codon_id] = genetic_code[[codons[codon_id]]]
    }
  }
  
  #Paste the codons together
  sequence <- paste(codons, collapse="")
  
  #Return the sequence
  return(sequence)
}
df$VDJ_germline_aa <- apply(df, 1, function(x) translate_DNA(x["VDJ_germline_nt"]))
df$VJ_germline_aa <- apply(df, 1, function(x) translate_DNA(x["VJ_germline_nt"]))
df$VDJ_sequence_aa <- apply(df, 1, function(x) translate_DNA(x["HC.nt.seq"]))
df$VJ_sequence_aa <- apply(df, 1, function(x) translate_DNA(x["LC.nt.seq"]))
df$full_sequence <- df$VDJ_sequence_aa

#Extend rows by clonal expansion
replicated_df <- df[rep(row.names(df), times = df$cells), ]

#Save as VDJ dataframe for pseudolikelihood calculations
write.csv(replicated_df, file = "../data/OVA_V7/VDJ_Variant_tree.csv", row.names = F)
