#Join the BCR heavy chain sequences with affinity measurements and transform nucleotide sequences into amino acids

library("readxl")
library(dplyr)

df_kim<-read.delim("../data/Kim/WU368_kim_et_al_nature_2022_bcr_heavy.tsv",sep="\t", header = T, row.names = NULL)
affinity_kim <- read_xlsx("../data/Kim/2022-12-05_ed_table_6.xlsx")
affinity_kim$sequence_id <- affinity_kim$h_id
affinity_kim$sequence_id %in% df_kim$sequence_id
new_df <- left_join(affinity_kim, df_kim, by = "sequence_id")

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
new_df$sequence_alignment<-gsub("\\.","-",new_df$sequence_alignment)
new_df$full_sequence <- apply(new_df, 1, function(x) translate_DNA(x["sequence_alignment"]))
new_df$full_sequence <-gsub("-","",new_df$full_sequence)
write.csv(new_df, file = "../data/Kim/Affinity_dataframe_HC.csv", row.names = FALSE)