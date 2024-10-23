#!/bin/bash
#SBATCH --job-name=affinity
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=100gb
#SBATCH --time=1-00:00:00
#SBATCH --output=log/affinity.out
#SBATCH --error=log/affinity.error


cd scripts

#ELISA results from Kim et al. (https://www.nature.com/articles/s41586-022-04527-1)
Rscript TransformKimELISA.R
python3 pll_calculator.py --dataset="Kim" --mode="general" --file_path="data/Kim/WU368_kim_bcr_heavy_elisa.csv"

#ELISA results from Neumeier et al. (https://www.pnas.org/doi/full/10.1073/pnas.2113766119)
python3 pll_calculator.py --dataset="OVA_V7" --mode="general" --file_path="data/OVA_V7/Binder_properties_cleaned.csv"

#Affinity results from Kim et al. (https://www.nature.com/articles/s41586-022-04527-1)
Rscript TransformKimAffinity.R
python3 pll_calculator.py --dataset="Kim" --mode="general" --file_path="data/Kim/Affinity_dataframe_HC.csv"

#Affinity results variant tree from Neumeier et al. (https://www.pnas.org/doi/full/10.1073/pnas.2113766119)
Rscript TransformOVAAffinity.R
python3 pll_calculator.py --dataset="OVA_V7" --mode="general" --file_path="/hpc/dla_lti/dvanginneken/PLM-likelihoods/data/OVA_V7/VDJ_Variant_tree.csv"
#Create lineage tree
Rscript VariantTree_default.R