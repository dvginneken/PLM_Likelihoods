#!/bin/bash
#SBATCH --job-name=edges
#SBATCH --nodes=1
#SBATCH --mem=24gb
#SBATCH --time=1-00:00:00
#SBATCH --output=log/edges.out
#SBATCH --error=log/edges.error

#Calculate the per-residue likelihoods of the mutations along the edges of the lineage trees

cd scripts
datasets=("OVA_V7" "horns2020a__VDJ_RAW" "Kim" "Bruhn")
for data in "${datasets[@]}"
  do
      #Transform AntibodyForests-object into csv
      Rscript transform_edge_dataframe.R $data "default"

      #Calculate the per-residue likelihoods of the mutations
      python3 clonal_evolution_prob.py -d $data -m "default"
    
      #Calculate the per-residue likelihoods of the original residues of the mutations
      python3 clonal_evolution_prob_reversed.py -d $data -m "default"

      #Calculate the per-residue likelihoods of the conserved residues along the edges
      python3 clonal_evolution_prob_conserved.py -d $data -m "default"

  done
echo "finished"