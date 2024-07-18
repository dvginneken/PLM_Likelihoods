#!/bin/bash
#SBATCH --job-name=evovelocity
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=100gb
#SBATCH --time=1-00:00:00
#SBATCH --output=log/evovelocity.out
#SBATCH --error=log/evovelocity.error

#Create the adata objects for the evolocity analysis for all sequences and for specific v-gene subfamilies

cd scripts

#All mouse samples
python3 evo_velo_adata_maker_OVA.py --dataset="OVA_V7" --model="esm" --input_source="full_VDJ"
#Specific v-gene subfamilies + IMGT germlines
genes=("IGHV3-1" "IGHV5-6" "IGHV9-3")
for gene in "${genes[@]}"
do
    python3 evo_velo_adata_maker_OVA_subfamily.py --model="esm" --v_gene=$gene
done

#All human samples
python3 evo_velo_adata_maker_human.py --model="esm" --input_source="full_VDJ"
#Specific v-gene subfamilies + IMGT germlines
genes=("IGHV3-23" "IGHV3-33" "IGHV4-59" "IGHV1-69D")
for gene in "${genes[@]}"
do
    python3 evo_velo_adata_maker_human_subfamily.py --model="esm" --v_gene=$gene
done