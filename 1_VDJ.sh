#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=VDJ
#SBATCH --error=log/CreateVDJ.error
#SBATCH --output=log/CreateVDJ.out

#Create VDJ-dataframe from CellRanger output using the Platypus package

cd scripts
datasets=("OVA_V7" "horns2020a__VDJ_RAW" "Kim" "Bruhn" "Kim_extra" "Mathew")
for data in "${datasets[@]}"
do
    Rscript CreateVDJ.R $data
done
