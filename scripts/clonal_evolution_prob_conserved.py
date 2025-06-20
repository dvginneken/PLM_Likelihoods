import pandas as pd
import numpy as np
import os
import sys
import argparse
import torch

sys.path.append("../src")

from ESMC_model import ESMc
from ablang2_model import Ablang2
from ablang_model import Ablang
from ESM1b_model import ESM1b
from sapiens_model import Sapiens
from protbert import ProtBert

#Set arguments
parser = argparse.ArgumentParser()
parser.add_argument('-d','--dataset')
parser.add_argument('-m','--method')   
args = parser.parse_args()
dataset = args.dataset
method = args.method

#Define models
model_names = ["ablang1","protbert","sapiens","esm1b","esmc","ablang2"]
model_classes = [Ablang,ProtBert, Sapiens,ESM1b,ESMc, Ablang2]

#Read edge file
edges_file = pd.read_csv(os.path.join("..","data",dataset,f"edges_{dataset}_{method}_HC.csv"))

#Output file
save_path = os.path.join("..","data",dataset,f"conserved_rank_table_{dataset}_{method}.csv")

#initiate output table
output_table = pd.DataFrame()

#Loop over each model
for j,model in enumerate(model_names):
    
    torch.cuda.empty_cache()
    
    #initiate model classes (heavy chain)
    if model == "ablang1":
        model_init = Ablang(chain="heavy")
    if model == "sapiens":
        model_init = Sapiens(chain_type="H")
    else:
        model_init = model_classes[j]()

    #Initiate output table for this model
    output_table_model = []



    #Loop over each edge in the edge file
    for i in range(len(edges_file)):
        
        #Extract sequences from edge file
        seq_1 = edges_file.iloc[i,2]
        seq_2 = edges_file.iloc[i,3]

        #Check if sequences are the same length
        if len(seq_1) != len(seq_2):
            continue
        
        #Get the position that do not mutate
        conserved_positions = [k for k in range(len(seq_1)) if seq_1[k] == seq_2[k]]
        mutating_positions = [k for k in range(len(seq_1)) if seq_1[k] != seq_2[k]]
        
        try:
             # If model is ablang2 add "|" to the sequence (HC and LC are separated with "|")
            if model == "ablang2":
                seq_1 = seq_1 + "|"

            #Calculate probability matrix (row is position in sequence, column is amino acid)
            prob_matrix = model_init.calc_probability_matrix(seq_1)
        
            #initiate rank list
            substitute_ranks = []

            #initiate probability list
            substitute_probabilities = []
            
            #Loop over the non-mutating positions
            for pos in conserved_positions:
                #Get the aa probabilities for this position
                likelihood_values = pd.Series(prob_matrix.iloc[pos,:])

                #Get the rank for each aa (rank 1 is highest probability)
                ranks = likelihood_values.rank(ascending=False)

                #Get the rank of the conserved residues over this edge
                mut_rank = ranks[seq_1[pos]]
                substitute_ranks.append(mut_rank)

                #Get the probability of the conserved residue over this edge
                probability = likelihood_values[seq_1[pos]]
                substitute_probabilities.append(probability)

            #get the average rank and probability
            mean_substitute_rank = np.average(substitute_ranks)   
            mean_substitute_probability = np.average(substitute_probabilities)
        except:
            #Sapiens issue sequence length fix:
            mean_substitute_rank = None
            mean_substitute_probability = None

        output_table_model.append({"model":model,"n_subs":len(mutating_positions), 
                                   "mean_sub_rank":mean_substitute_rank, 
                                   "mean_sub_prob":mean_substitute_probability})

    #Convert dictionary list to dataframe
    output_table_model = pd.DataFrame.from_dict(output_table_model)
    #Concatenate the columns of the edge file with the models' output table
    output_table_model = pd.concat([edges_file, output_table_model], axis=1)
    #Concatenate the rows of the model's output table to the final output table
    output_table = pd.concat([output_table, output_table_model])

    del(model_init)
          
output_table = pd.DataFrame(output_table)   
output_table.to_csv(save_path, index=False)
