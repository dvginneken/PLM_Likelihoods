#!/bin/bash
#SBATCH --time=5-00:00:00
#SBATCH --mem=32G
#SBATCH --partition=gpu
#SBATCH --job-name=LineageTrees
#SBATCH --error=log/LineageTrees.error
#SBATCH --output=log/LineageTrees.out

cd scripts
datasets=("OVA_V7" "horns2020a__VDJ_RAW" "Bruhn" "Kim")
methods=("phylo.network.default" "phylo.tree.mp")
for data in "${datasets[@]}"
  do
      for method in "${methods[@]}"
        do
          Rscript LineageTrees.R $data $method
        done
  done
echo "finished"