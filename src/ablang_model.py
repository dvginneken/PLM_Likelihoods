import ablang
import numpy as np
import pandas as pd
import pickle as pkl
import scipy
from tqdm import tqdm
import os
import sys
import torch

sys.path.append("../scripts")

from utils import get_pseudo_likelihood
class Ablang():

    """
    Class for the protein Model Ablang
    """

    def __init__(self, chain = "heavy",file_name = ".", method = "seqcoding"):
        """
        Creates the instance of the language model instance; either light or heavy
        
        method: `str`
        Which token to use to extract embeddings of the last layer
        
        file_name: `str`
        The name of the folder to store the embeddings
        """
        self.device = "cuda:0" if torch.cuda.is_available() else "cpu"
        self.model = ablang.pretrained(chain,device=self.device)
        #dont update the weights
        self.model.freeze()

        self.file = file_name
        self.mode = method
    


    def fit_transform(self, sequences, starts, ends):
        """
        Fits the model and outputs the embeddings.
        
        parameters
        ----------

        sequences: `list` 
        List with sequences to be transformed
        ------

        None, saved the embeddings in the embeddings.csv
        """
        output = self.model(sequences, mode=self.mode)
        if self.mode == "seqcoding":
            #The embeddings are made my averaging across all residues    
            # pd.DataFrame(output).to_csv("outfiles/"+self.file+"/embeddings.csv")
            return pd.DataFrame(output,columns=[f"dim_{i}" for i in range(output.shape[1])])

    def calc_evo_likelihood_matrix_per_position(self, sequences:list):

        probs = []
        for sequence in tqdm(sequences):
            logits = self.model(sequence, mode="likelihood")[0]
            prob = scipy.special.softmax(logits,axis = 1)
            df = pd.DataFrame(prob, columns = list(self.model.tokenizer.vocab_to_aa.values())[4:])
            #removing CLS and SEP
            df = df.iloc[1:-1,:]
            df = df.reindex(sorted(df.columns), axis=1)
            probs.append(df)

        likelihoods = get_pseudo_likelihood(probs, sequences) 
        pkl.dump([probs,likelihoods],open("outfiles/"+self.file+"/probabilities_pseudo.pkl","wb"))

    def calc_pseudo_likelihood_sequence(self, sequences: list, starts, ends):
        pll_all_sequences = []
        for j,sequence in enumerate(tqdm(sequences)):
            try:
                amino_acids = list(sequence)
                logits = self.model(sequence, mode="likelihood")[0]
                prob = scipy.special.softmax(logits,axis = 1)
                df = pd.DataFrame(prob, columns = list(self.model.tokenizer.vocab_to_aa.values())[4:])
                df = df.iloc[1:-1,:]
                df = df.reindex(sorted(df.columns), axis=1)


                per_position_ll = []
                for i in range(starts[j],ends[j]):
                    aa_i = amino_acids[i]
                    if aa_i == "-" or aa_i == "*":
                        continue
                    ll_i = np.log(df.iloc[i,:][aa_i])
                    per_position_ll.append(ll_i)
            

                pll_seq = np.average(per_position_ll)
                pll_all_sequences.append(pll_seq)
            except:
                pll_all_sequences.append(None)

        return pll_all_sequences

    def calc_probability_matrix(self, sequence:str):
        logits = self.model(sequence, mode="likelihood")[0]
        prob = scipy.special.softmax(logits,axis = 1)
        df = pd.DataFrame(prob, columns = list(self.model.tokenizer.vocab_to_aa.values())[4:])
        df = df.iloc[1:-1,:]
        df = df.reindex(sorted(df.columns), axis=1)

        return df