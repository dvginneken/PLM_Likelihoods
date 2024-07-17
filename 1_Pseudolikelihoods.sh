#!/bin/bash
#SBATCH --job-name=pll_calculator
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=100gb
#SBATCH --time=1-00:00:00
#SBATCH --output=log/pll_calculator.out
#SBATCH --error=log/pll_calculator.error

#Calculate the pseudolikelihoods for all datasets for the three modes: cdr3_only, full_VDJ, cdr3_from_VDJ
#Concatenate the results and calculate the correlations (for heavy and light chains seperately and combined)

cd scripts
modes=("cdr3_only" "full_VDJ" "cdr3_from_VDJ")
datasets=("OVA_V7" "horns2020a__VDJ_RAW" "Bruhn" "Kim")
for data in "${datasets[@]}"
do
    for mode in "${modes[@]}"
    do
        python3 pll_calculator.py --dataset=$data --mode=$mode
    done
    Rscript pll_concatenate.R $data
    Rscript SourceCorrelation_chains.R $data
    Rscript PLMCorrelation_chains.R $data
    Rscript SourceCorrelation.R $data
    Rscript PLMCorrelation.R $data
done