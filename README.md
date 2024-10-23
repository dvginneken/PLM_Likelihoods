# PLM Likelihoods
<img width="292" alt="image" src="https://github.com/user-attachments/assets/3eed26ce-0674-47d7-8229-8aed5b0b55aa">

This repository contains the code to reproduce the results from the manuscript "Protein Language Models capture features of in vivo B cell selection and evolution".
In this study, we explore foundational and antibody-specific PLMs on partially labeled human and mouse single-cell BCR data in the context of PLM-based pseudolikelihoods.

### Set up environment
1. Install conda environment
   `conda env create -f environment.yml`
2. Install SRA-toolkit https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit
3. Install Cell Ranger 8.0.0 https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-in  
   Install the reference data https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0.tar.gz
3. Install the R-package Platypus https://github.com/alexyermanos/Platypus

### Input files
Four publicly available datasets were used for this analysis, all CellRanger output files should be saved in PLM_Likelihoods/data/[dataset_name]/VDJ/[cellranger_output]. 
1. The five mice samples from Neumeier et al. 2023 can be downloaded from the Platypus Database (https://alexyermanos.github.io/Platypus/PlatypusDB_index.html) with project ID “neumeier2021b”
2. Individual1, from Horns et al. 2020 can be downloaded from the Platypus Database (https://alexyermanos.github.io/Platypus/PlatypusDB_index.html) with project ID “horns2020a”
   only the sample “Influenza.vac.11.12.human.S1” was used here
3. Individual2, from Bruhn et al. 2024 was downloaded from the NCBI Sequence Read Archive and processed with CellRanger VDJ pipeline (v8.0.0) with the script 0_CellRanger.sh
4. Individual3-5, from Kim et al. 2022 was downloaded from the NCBI Sequence Read Archive and processed with CellRanger VDJ pipeline (v8.0.0) with the script 0_CellRanger.sh

### Run code
1. Run the bash scripts in order 0 to 5 for all calculations.
2. To reproduce the figures from the manuscript, use the code in the folder "figures".
