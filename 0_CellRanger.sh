#!/bin/bash
#SBATCH --time=2-00:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=CellRanger
#SBATCH --error=log/CellRanger.error
#SBATCH --output=log/CellRanger.out

#Download fastq files via SRA-toolkit for the Bruhn and Kim dataset and run the CellRanger (v8.0.0) V(D)J pipeline.

#Bruhn data
cd /PLM_likelihoods/data/Bruhn/VDJ
prefetch SRR24140776 --output-directory fastq/
fastq-dump --origfmt -I --split-files fastq/SRR24140776/SRR24140776.sra --outdir fastq/SRR24140776/
mv fastq/SRR24140776/SRR24140776_1.fastq fastq/SRR24140776/SRR24140776_S1_L001_R1_001.fastq
mv fastq/SRR24140776/SRR24140776_2.fastq fastq/SRR24140776/SRR24140776_S1_L001_R2_001.fastq
gzip fastq/SRR24140776/SRR24140776_S1_L001_R1_001.fastq
gzip fastq/SRR24140776/SRR24140776_S1_L001_R2_001.fastq
cellranger vdj --id=SRR24140776\
         --reference=/hpc/dla_lti/dvanginneken/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0 \
         --fastqs=fastq/SRR24140776 \
         --sample=SRR24140776
echo "Done SRR24140776"

#Kim data
cd /PLM_likelihoods/data/Kim/VDJ
identifiers=("SRR17729726" "SRR17729703" "SRR17729692")
for identifiers in "${srr[@]}"
do
    prefetch $srr --max-size 100G --output-directory fastq/
    fastq-dump --origfmt -I --split-files fastq/${srr}/${srr}.sra --outdir fastq/${srr}
    mv fastq/${srr}/${srr}_1.fastq fastq/${srr}/${srr}_S1_L001_R1_001.fastq
    mv fastq/${srr}/${srr}_2.fastq fastq/${srr}/${srr}_S1_L001_R2_001.fastq
    gzip fastq/${srr}/${srr}_S1_L001_R1_001.fastq
    gzip fastq/${srr}/${srr}_S1_L001_R2_001.fastq
    cellranger vdj --id=${srr}\
         --reference=/hpc/dla_lti/dvanginneken/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0 \
         --fastqs=fastq/${srr} \
         --sample=${srr}
    echo "Done ${srr}"
done