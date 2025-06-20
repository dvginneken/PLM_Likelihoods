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
python3 pll_calculator.py --dataset="Kim" --mode="general" -hc="sequence_HC" -lc="sequence_LC_aa" \
    --file_path="data/Kim/Kim_elisa.csv" --output_dir="data/Kim/" --file_name="Kim_elisa"

#ELISA results from Neumeier et al. (https://www.pnas.org/doi/full/10.1073/pnas.2113766119)
Rscript TransformOVAELISA.R
python3 pll_calculator.py --dataset="OVA_V7" --mode="general" --file_path="data/OVA_V7/Binder_properties_cleaned.csv"

#Affinity results from Kim et al. (https://www.nature.com/articles/s41586-022-04527-1)
Rscript TransformKimAffinity.R
python3 pll_calculator.py --dataset="Kim" --mode="general"  -hc="full_sequence" -lc="LC_sequence" \
    --file_path="data/Kim/Affinity_dataframe.csv" --output_dir="data/Kim/" --file_name="Kim_affinity"

#Affinity results variant tree from Neumeier et al. (https://www.pnas.org/doi/full/10.1073/pnas.2113766119)
Rscript TransformOVAAffinity.R
python3 pll_calculator.py --dataset="OVA_V7" --mode="general"  -hc="VDJ_sequence_aa" -lc="VJ_sequence_aa" \
    --file_path="data/OVA_V7/Affinity/VDJ_Variant_tree.csv" --output_dir="data/OVA_V7/Affinity/" --file_name="OVA_affinity"

#Create lineage tree
Rscript VariantTree_default.R