#Transform processed BCR heavy chain data from https://zenodo.org/records/5895181 to input csv for the evolikelihood analysis
df<-read.delim("../data/Kim/WU368_kim_et_al_nature_2022_bcr_heavy.tsv",sep="\t")

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

#only keep sequences with ELISA results
df_elisa <- df[!is.na(df$elisa),]

#Translate nt to aa for the VDJ sequence
df_elisa$sequence_alignment<-gsub("\\.","-",df_elisa$sequence_alignment)
df_elisa$full_sequence <- apply(df_elisa, 1, function(x) translate_DNA(x["sequence_alignment"]))

#Prepare dataframe for further analysis
df_elisa$barcode <- df_elisa$cell_id
df_elisa$contig_id <- df_elisa$sequence_id
df_elisa$chain <- df_elisa $locus
df_elisa$v_gene <- df_elisa$v_call
df_elisa$d_gene <- df_elisa$d_call
df_elisa$j_gene <- df_elisa$j_call
df_elisa$c_gene <- df_elisa$isotype
df_elisa$clonotype_id <- df_elisa$clone_id
df_elisa$sample_id <- df_elisa$donor

#remove gaps
df_elisa$full_sequence <-gsub("-","",df_elisa$full_sequence)

#save to csv
write.csv(df_elisa, file = "../data/Kim/WU368_kim_bcr_heavy_elisa.csv", row.names = FALSE)
