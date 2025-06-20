#!/bin/bash
#SBATCH --time=5-00:00:00
#SBATCH --mem=32G
#SBATCH --partition=gpu
#SBATCH --job-name=LineageTrees
#SBATCH --error=log/LineageTrees.error
#SBATCH --output=log/LineageTrees.out

#Construct lineage trees for all clonotypes in the datasets with an mst-like algorithm (default) and with a maximum parsimony algorithm

cd scripts
datasets=("OVA_V7" "horns2020a__VDJ_RAW" "Bruhn" "Kim" "Kim_extra" "Mathew")
methods=("phylo.network.default" "phylo.tree.mp")
for data in "${datasets[@]}"
  do
      for method in "${methods[@]}"
        do
          Rscript LineageTrees.R $data $method
        done
  done
echo "finished"