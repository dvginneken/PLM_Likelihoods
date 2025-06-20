import pandas as pd
import numpy as np
import os
import sys
import argparse

sys.path.append("../src")

parser = argparse.ArgumentParser()
parser.add_argument('-d','--dataset')   
parser.add_argument('-m','--mode') #PLM input source
parser.add_argument('-o','--output_dir') # output directory
parser.add_argument('-hc', '--heavy_chain', default="full_sequence") # column name in the input file with heavy chain sequences
parser.add_argument('-lc', '--light_chain', default="full_sequence") # column name in the input file with light chain sequences
parser.add_argument('-fp','--file_path') #path to the file with sequences
parser.add_argument('-f','--file_name', default="output") #name of the output file
args = parser.parse_args()
dataset = args.dataset
mode = args.mode
output_dir = args.output_dir
file_name = args.file_name

suffixes = ["protbert","esm1b","esmc","ablang","ablang2","sapiens"]
if mode == "general":
    save_path = os.path.join(output_dir,"evo_likelihoods")

    if not os.path.exists(save_path):
        os.mkdir(save_path)

    for model_suffix in suffixes:
        repertoire_file  = pd.read_csv(args.file_path)
        print(f"Calculating evo likelihoods for {model_suffix} model")
        if model_suffix != "ablang2":
            if model_suffix == "ablang":
                from ablang_model import Ablang
                model = Ablang(chain = "heavy")
            elif model_suffix == "sapiens":
                from sapiens_model import Sapiens
                model = Sapiens(chain_type="H")
            elif model_suffix == "protbert":
                from protbert import ProtBert
                model = ProtBert()
            elif model_suffix == "esm1b":
                from ESM1b_model import ESM1b
                model = ESM1b()
            elif model_suffix == "esmc":
                from ESMC_model import ESMc
                model = ESMc()
            #only heavy chain
            sequences_column_HC = args.heavy_chain
            starts_hc = [0]*repertoire_file.shape[0]
            ends_hc = repertoire_file[sequences_column_HC].apply(len)
            #calculate likelihoods
            repertoire_file["evo_likelihood"] = model.calc_pseudo_likelihood_sequence(list(repertoire_file[sequences_column_HC]),list(starts_hc),list(ends_hc))

        elif model_suffix == "ablang2":
            from ablang2_model import Ablang2
            model = Ablang2()
            #heavy chain
            sequences_column_HC = args.heavy_chain
            starts_hc = [0]*repertoire_file.shape[0]
            ends_hc = repertoire_file[sequences_column_HC].apply(len)
            #light chain
            sequence_column_LC = args.light_chain
            starts = [0]*repertoire_file.shape[0]
            repertoire_file[sequence_column_LC] = repertoire_file[sequence_column_LC].fillna('')
            ends = repertoire_file[sequence_column_LC].apply(len)
            #calculate likelihoods
            hc_evo, lc_evo, paired_evo = model.calc_pseudo_likelihood_sequence(sequences=list(repertoire_file[sequences_column_HC] + "|" + repertoire_file[sequence_column_LC]),
                                                                                                        hc_starts=list(starts_hc),
                                                                                                        hc_ends=list(ends_hc),
                                                                                                        lc_starts=list(starts),
                                                                                                        lc_ends=list(ends))
            repertoire_file["HC_evo_likelihood"] = hc_evo
            repertoire_file["LC_evo_likelihood"] = lc_evo
            repertoire_file["pair_evo_likelihood"] = paired_evo
        #Save output
        repertoire_file.to_csv(os.path.join(save_path,f"{file_name}_evo_likelihood_{model_suffix}.csv"), index=False)  
else:
    data_folder_path = os.path.join("..","data",dataset,"VDJ") ## folder containing all the per-sample folders
    repertoire_file = pd.read_csv(os.path.join("..","data",dataset,"VDJ_"+ dataset +".csv")) # load the VDJ file
    for sample in os.listdir(data_folder_path):
        evo_folder = os.path.join(data_folder_path, sample,"evo_likelihoods")
        if not (os.path.isdir(evo_folder)): ## check if evo_folder already exists
            os.mkdir(evo_folder)

        
        
        for model_suffix in suffixes:
            repertoire_file_sample = repertoire_file[repertoire_file["sample_id"] == sample] # subset the VDJ file per sample
            print(f"Calculating evo likelihoods for {model_suffix} model")

            if model_suffix != "ablang2": # ablang2 calculate the likelihoods for HC and LC together
                print("not ablang2")
                if model_suffix == "ablang": # ablang and sapiens have separate PLMs for heavy and light chains
                    from ablang_model import Ablang
                    HC_model = Ablang(chain="heavy")
                    LC_model = Ablang(chain="light")
                elif model_suffix == "sapiens":
                    from sapiens_model import Sapiens
                    HC_model = Sapiens(chain_type="H")
                    LC_model = Sapiens(chain_type="L")
                elif model_suffix == "protbert":
                    from protbert import ProtBert
                    HC_model = ProtBert()
                    LC_model = ProtBert()
                elif model_suffix == "esm1b":
                    from ESM1b_model import ESM1b
                    HC_model = ESM1b()
                    LC_model = ESM1b()
                elif model_suffix == "esmc":
                    from ESMC_model import ESMc
                    HC_model = ESMc()
                    LC_model = ESMc()
                else:
                    raise ValueError(f"Model {model_suffix} not recognized.")
                
                print("model:")
                print(HC_model)
                print(LC_model)

                if (mode == "full_VDJ"):  
                    save_path = os.path.join(evo_folder,"full_VDJ")
                    if not (os.path.isdir(save_path)):
                        os.mkdir(save_path)

                    repertoire_file_sample["HC_evo_likelihood"] = HC_model.calc_pseudo_likelihood_sequence(list(repertoire_file_sample["VDJ_sequence_aa_trimmed"]),
                                                                                                        list(pd.Series([0]*repertoire_file_sample.shape[0])),
                                                                                                        list(repertoire_file_sample["VDJ_sequence_aa_trimmed"].apply(lambda x: len(x) if not pd.isna(x) else None)))
                    repertoire_file_sample["LC_evo_likelihood"] = LC_model.calc_pseudo_likelihood_sequence(list(repertoire_file_sample["VJ_sequence_aa_trimmed"]),
                                                                                                        list(pd.Series([0]*repertoire_file_sample.shape[0])),
                                                                                                        list(repertoire_file_sample["VJ_sequence_aa_trimmed"].apply(lambda x: len(x) if not pd.isna(x) else None)))
                if (mode == "cdr3_only"):  
                    save_path = os.path.join(evo_folder,"cdr3_only")
                    if not (os.path.isdir(save_path)):
                        os.mkdir(save_path)

                    repertoire_file_sample["HC_evo_likelihood"] = HC_model.calc_pseudo_likelihood_sequence(list(repertoire_file_sample["VDJ_cdr3_aa"]),
                                                                                                        list(pd.Series([0]*repertoire_file_sample.shape[0])),
                                                                                                        list(repertoire_file_sample["VDJ_cdr3_aa"].apply(lambda x: len(x) if not pd.isna(x) else None)))
                    repertoire_file_sample["LC_evo_likelihood"] = LC_model.calc_pseudo_likelihood_sequence(list(repertoire_file_sample["VJ_cdr3_aa"]),
                                                                                                        list(pd.Series([0]*repertoire_file_sample.shape[0])),
                                                                                                        list(repertoire_file_sample["VJ_cdr3_aa"].apply(lambda x: len(x) if not pd.isna(x) else None)))
                if mode == "cdr3_from_VDJ":  ### use only subset of positions defined by "starts" and "ends" - CDR3 region in this case
                    save_path = os.path.join(evo_folder,"cdr3_from_VDJ")
                    if not (os.path.isdir(save_path)):
                        os.mkdir(save_path)

                    hc_x = repertoire_file_sample["VDJ_fwr1_aa"] + repertoire_file_sample["VDJ_cdr1_aa"] + repertoire_file_sample["VDJ_fwr2_aa"] + \
                        repertoire_file_sample["VDJ_cdr2_aa"] + repertoire_file_sample["VDJ_fwr3_aa"]
                    hc_y = repertoire_file_sample["VDJ_fwr1_aa"] + repertoire_file_sample["VDJ_cdr1_aa"] + repertoire_file_sample["VDJ_fwr2_aa"] + \
                        repertoire_file_sample["VDJ_cdr2_aa"] + repertoire_file_sample["VDJ_fwr3_aa"] + repertoire_file_sample["VDJ_cdr3_aa"]  
                    repertoire_file_sample["HC_evo_likelihood"] = HC_model.calc_pseudo_likelihood_sequence(list(repertoire_file_sample["VDJ_sequence_aa_trimmed"]),
                                                                                                        list(hc_x.apply(lambda hc_x: len(hc_x) if not pd.isna(hc_x) else None)),
                                                                                                        list(hc_y.apply(lambda hc_y: len(hc_y) if not pd.isna(hc_y) else None)))
                    
                    lc_x = repertoire_file_sample["VJ_fwr1_aa"] + repertoire_file_sample["VJ_cdr1_aa"] + repertoire_file_sample["VJ_fwr2_aa"] + \
                        repertoire_file_sample["VJ_cdr2_aa"] + repertoire_file_sample["VJ_fwr3_aa"]
                    lc_y = repertoire_file_sample["VJ_fwr1_aa"] + repertoire_file_sample["VJ_cdr1_aa"] + repertoire_file_sample["VJ_fwr2_aa"] + \
                        repertoire_file_sample["VJ_cdr2_aa"] + repertoire_file_sample["VJ_fwr3_aa"] + repertoire_file_sample["VJ_cdr3_aa"]  
                    repertoire_file_sample["LC_evo_likelihood"] = LC_model.calc_pseudo_likelihood_sequence(list(repertoire_file_sample["VJ_sequence_aa_trimmed"]),
                                                                                                        list(lc_x.apply(lambda lc_x: len(lc_x) if not pd.isna(lc_x) else None)),
                                                                                                        list(lc_y.apply(lambda lc_y: len(lc_y) if not pd.isna(lc_y) else None)))
        
            elif model_suffix == "ablang2":
                from ablang2_model import Ablang2
                model = Ablang2()

                print("model:")
                print(model)

                if (mode == "full_VDJ"):  
                    save_path = os.path.join(evo_folder,"full_VDJ")
                    if not (os.path.isdir(save_path)):
                        os.mkdir(save_path)

                    hc_evo, lc_evo, paired_evo = model.calc_pseudo_likelihood_sequence(sequences=list(repertoire_file_sample["VDJ_sequence_aa_trimmed"] + "|" + repertoire_file_sample["VJ_sequence_aa_trimmed"]),
                                                                                                        hc_starts=list(pd.Series([0]*repertoire_file_sample.shape[0])),
                                                                                                        hc_ends=list(repertoire_file_sample["VDJ_sequence_aa_trimmed"].apply(lambda x: len(x) if not pd.isna(x) else None)),
                                                                                                        lc_starts=list(pd.Series([0]*repertoire_file_sample.shape[0])),
                                                                                                        lc_ends=list(repertoire_file_sample["VJ_sequence_aa_trimmed"].apply(lambda x: len(x) if not pd.isna(x) else None)))
                    repertoire_file_sample["HC_evo_likelihood"] = hc_evo
                    repertoire_file_sample["LC_evo_likelihood"] = lc_evo
                    repertoire_file_sample["pair_evo_likelihood"] = paired_evo
                if (mode == "cdr3_only"):  
                    save_path = os.path.join(evo_folder,"cdr3_only")
                    if not (os.path.isdir(save_path)):
                        os.mkdir(save_path)

                    hc_evo, lc_evo, paired_evo = model.calc_pseudo_likelihood_sequence(sequences=list(repertoire_file_sample["VDJ_cdr3_aa"] + "|" + repertoire_file_sample["VJ_cdr3_aa"]),
                                                                                                        hc_starts=list(pd.Series([0]*repertoire_file_sample.shape[0])),
                                                                                                        hc_ends=list(repertoire_file_sample["VDJ_cdr3_aa"].apply(lambda x: len(x) if not pd.isna(x) else None)),
                                                                                                        lc_starts=list(pd.Series([0]*repertoire_file_sample.shape[0])),
                                                                                                        lc_ends=list(repertoire_file_sample["VJ_cdr3_aa"].apply(lambda x: len(x) if not pd.isna(x) else None)))
                    repertoire_file_sample["HC_evo_likelihood"] = hc_evo
                    repertoire_file_sample["LC_evo_likelihood"] = lc_evo
                    repertoire_file_sample["pair_evo_likelihood"] = paired_evo
                if mode == "cdr3_from_VDJ":
                    save_path = os.path.join(evo_folder,"cdr3_from_VDJ")
                    if not (os.path.isdir(save_path)):
                        os.mkdir(save_path)

                    hc_x = repertoire_file_sample["VDJ_fwr1_aa"] + repertoire_file_sample["VDJ_cdr1_aa"] + repertoire_file_sample["VDJ_fwr2_aa"] + \
                        repertoire_file_sample["VDJ_cdr2_aa"] + repertoire_file_sample["VDJ_fwr3_aa"]
                    hc_y = repertoire_file_sample["VDJ_fwr1_aa"] + repertoire_file_sample["VDJ_cdr1_aa"] + repertoire_file_sample["VDJ_fwr2_aa"] + \
                        repertoire_file_sample["VDJ_cdr2_aa"] + repertoire_file_sample["VDJ_fwr3_aa"] + repertoire_file_sample["VDJ_cdr3_aa"]  
                    lc_x = repertoire_file_sample["VJ_fwr1_aa"] + repertoire_file_sample["VJ_cdr1_aa"] + repertoire_file_sample["VJ_fwr2_aa"] + \
                        repertoire_file_sample["VJ_cdr2_aa"] + repertoire_file_sample["VJ_fwr3_aa"]
                    lc_y = repertoire_file_sample["VJ_fwr1_aa"] + repertoire_file_sample["VJ_cdr1_aa"] + repertoire_file_sample["VJ_fwr2_aa"] + \
                        repertoire_file_sample["VJ_cdr2_aa"] + repertoire_file_sample["VJ_fwr3_aa"] + repertoire_file_sample["VJ_cdr3_aa"]  
                    
                    hc_evo, lc_evo, paired_evo = model.calc_pseudo_likelihood_sequence(sequences=list(repertoire_file_sample["VDJ_sequence_aa_trimmed"] + "|" + repertoire_file_sample["VJ_sequence_aa_trimmed"]),
                                                                                                        hc_starts=list(hc_x.apply(lambda hc_x: len(hc_x) if not pd.isna(hc_x) else None)),
                                                                                                        hc_ends=list(hc_y.apply(lambda hc_y: len(hc_y) if not pd.isna(hc_y) else None)),
                                                                                                        lc_starts=list(lc_x.apply(lambda lc_x: len(lc_x) if not pd.isna(lc_x) else None)),
                                                                                                        lc_ends=list(lc_y.apply(lambda lc_y: len(lc_y) if not pd.isna(lc_y) else None)))
                    repertoire_file_sample["HC_evo_likelihood"] = hc_evo
                    repertoire_file_sample["LC_evo_likelihood"] = lc_evo
                    repertoire_file_sample["pair_evo_likelihood"] = paired_evo
                                                                                                                                

            repertoire_file_sample.to_csv(os.path.join(save_path,f"evo_likelihood_{model_suffix}.csv"), index=False)
