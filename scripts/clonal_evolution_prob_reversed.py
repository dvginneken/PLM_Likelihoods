import pandas as pd
import numpy as np
import os
import sys
import argparse
import torch

sys.path.append("../src")

from ablang_model import Ablang
from ESM_model import ESM
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
model_names = ["ablang","protbert","sapiens","esm"]
model_classes = [Ablang,ProtBert, Sapiens,ESM]

#Read edge file
edges_file = pd.read_csv(os.path.join("..","data",dataset,f"edges_{dataset}_{method}_HC.csv"))

#Output file
save_path = os.path.join("..","data",dataset,f"mutational_rank_reversed_table_{dataset}_{method}.csv")

#initiate output table
output_table = pd.DataFrame()

#Loop over each model
for j,model in enumerate(model_names):
    
    torch.cuda.empty_cache()
    
    #initiate model classes (heavy chain)
    if model == "ablang":
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
        
        #Get the position of mutations
        diff_positions = [k for k in range(len(seq_1)) if seq_1[k] != seq_2[k]]
        
        try:
            #Calculate probability matrix (row is position in sequence, column is amino acid)
            prob_matrix = model_init.calc_probability_matrix(seq_1)
        
            #initiate rank list
            substitute_ranks = []

            #initiate probability list
            substitute_probabilities = []
            
            #Loop over the mutational positions
            for pos in diff_positions:
                #Get the aa probabilities for this position
                likelihood_values = pd.Series(prob_matrix.iloc[pos,:])

                #Get the rank for each aa (rank 1 is highest probability)
                ranks = likelihood_values.rank(ascending=False)

                #Get the rank of the original residue over this edge
                mut_rank = ranks[seq_1[pos]]
                substitute_ranks.append(mut_rank)

                #Get the probability of the original residue over this edge
                probability = likelihood_values[seq_1[pos]]
                substitute_probabilities.append(probability)

            #If there are multiple mutations in a sequence, get the average rank    
            mean_substitute_rank = np.average(substitute_ranks)

            #If there are multiple mutations in a sequence, get the average rank    
            mean_substitute_probability = np.average(substitute_probabilities)
        except:
            #Sapiens issue sequence length fix:
            mean_substitute_rank = None
            mean_substitute_probability = None

        output_table_model.append({"model":model,"n_subs":len(diff_positions), 
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