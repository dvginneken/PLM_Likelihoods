library(readr)
library(dplyr)

#Read data
WU368_kim_bcr_heavy_elisa <- read_csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim/WU368_kim_bcr_heavy_elisa.csv")
WU368_kim_bcr_light_elisa <- read_csv("/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim/WU368_kim_bcr_light_elisa.csv")


#Join heavy + light chain
kim_elisa <- left_join(
  WU368_kim_bcr_heavy_elisa[,c("cell_id", "full_sequence", "elisa", "sample_id")],
  WU368_kim_bcr_light_elisa[,c("cell_id", "sequence_alignment")],
  by = "cell_id"
)
colnames(kim_elisa) <- c("cell_id", "sequence_HC", "elisa", "sample_id", "sequence_LC")

#Remove sequences with N's
kim_elisa <- kim_elisa[grep("N", kim_elisa$sequence_LC, invert = T),]

#Remove dots from the light chain sequences
kim_elisa$sequence_LC <- gsub("\\.", "", kim_elisa$sequence_LC)

#Translate light chain nucleotide sequences to amino acid sequences
translate_DNA<- function(sequence){
  
  #Translate a nucleotide sequence into an amino acid sequence
  #Arguments:
  #- sequence: nucleotide sequence to be translated
  
  if (sequence == ""){
    return("")
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
  
kim_elisa$sequence_LC_aa <- apply(kim_elisa, 1, function(x) translate_DNA(x["sequence_LC"]))

#Manually change out-of-frame sequences
kim_elisa[grep("\\*", kim_elisa$sequence_LC_aa),] -> df

remove_first <- c("368-02a_s40@ACGAGCCAGACCACGA-1", "368-02a_s41@ACATCAGTCGCTGATA-1", "368-02a_s42@GACGTGCGTTAAGAAC-1",
"368-02a_s37@GAACATCTCTCGCATC-1", "368-02a_s40@TATCTCATCGTTACAG-1", "368-02a_s37@ATCCACCAGGTAGCTG-1",
"368-02a_s40@GGCGACTGTCTTGTCC-1", "368-13_s33@CGTGAGCTCCGGCACA-1", "368-20_s27@GACGCGTTCTTTACGT-1",
"368-02a_s38@CATCAGAGTGGCTCCA-1", "368-01a_s15@TTGTAGGCATCGTCGG-1", "368-02a_s38@TCGGGACCAAGCCCAC-1",
"368-04_s23@AACTCTTAGATAGTCA-1", "368-04_s24@CGACCTTAGTATCGAA-1", "368-10_s44@TTCGAAGAGACCACGA-1",
"368-22_s18@GAGTCCGTCCACGACG-1")

for (cell_id in remove_first){
  kim_elisa[kim_elisa$cell_id == cell_id, "sequence_LC"] <- substr(kim_elisa[kim_elisa$cell_id == cell_id, "sequence_LC"],
                                                                 start = 2, stop = nchar(kim_elisa[kim_elisa$cell_id == cell_id, "sequence_LC"]))
}

remove_first2 <- c("368-02a_s38@CACCTTGAGGATGCGT-1", "368.07.28.GC.3D07", "368-07_s31@CCTACCATCAGCTCTC-1",
                   "368-07_s30@GACGTGCAGAGACTTA-1", "368-07_s32@AATCGGTCAGGCTGAA-1", "368-10_s45@GGGCATCTCAGTTCGA-1",
                   "368-22_s18@GAGTCCGTCCACGACG-1", "368-22_s18@GAGTCCGTCCACGACG-1", "368-13_s35@TACTTGTCACCTGGTG-1",
                   "368-20_s28@GGAACTTAGGGTGTGT-1", "368-22_s17@TACTCGCTCAAACCAC-1")

for (cell_id in remove_first2){
  kim_elisa[kim_elisa$cell_id == cell_id, "sequence_LC"] <- substr(kim_elisa[kim_elisa$cell_id == cell_id, "sequence_LC"],
                                                                 start = 3, stop = nchar(kim_elisa[kim_elisa$cell_id == cell_id, "sequence_LC"]))
}

#Translate again
kim_elisa$sequence_LC_aa <- apply(kim_elisa, 1, function(x) translate_DNA(x["sequence_LC"]))

write.csv(kim_elisa, "/Users/dginneke/Documents/GitHub/PLM-likelihoods/data/Kim/Kim_elisa.csv", row.names = F)
