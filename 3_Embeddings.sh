#!/bin/bash
#SBATCH --job-name=emb_calculator
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=100gb
#SBATCH --time=2-00:00:00
#SBATCH --output=log/emb_calculator.out
#SBATCH --error=log/emb_calculator.error

#Calculate the embeddings for all sequences and for specific IMGT germline v-gene sequences

cd scripts

#Calculate the embeddings of all sequences in the datasets
modes=("cdr3_only" "full_VDJ" "cdr3_from_VDJ" "v_gene_only")
datasets=("OVA_V7" "horns2020a__VDJ_RAW" "Bruhn" "Kim")
for data in "${datasets[@]}"
do
    for mode in "${modes[@]}"
    do
        python3 embeddings_calculator.py --dataset=$data --mode=$mode
        echo "done $data $mode"
    done
done

#Calculate the embeddings of specific IMGT germline v-gene sequences
python3 embeddings_calculator.py --dataset="OVA_V7" --mode="general" --file_path="/hpc/dla_lti/dvanginneken/PLM-likelihoods/data/OVA_IMGT_germlines.csv"
python3 embeddings_calculator.py --dataset="Human" --mode="general" --file_path="/hpc/dla_lti/dvanginneken/PLM-likelihoods/data/Human_IMGT_germlines.csv"