import matplotlib
import torch
from evolocity.tools.fb_model import FBModel
import argparse
import pandas as pd
import anndata
import evolocity as evo
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import os
parser = argparse.ArgumentParser()
parser.add_argument("-m", "--model", help="Choose a model: esm, ablang, protbert, sapiens")
parser.add_argument("-i", "--input_source", help="cdr3_only or full_VDJ")

args = parser.parse_args()
model = args.model
input_source = args.input_source

all_embeddings = list()

## Per dataset
human_datasets = ["Bruhn", "Kim", "horns2020a__VDJ_RAW"]
for dataset in human_datasets:
    ## Per sample
    for sample_id in os.listdir(os.path.join("..","data",dataset,"VDJ")):
        #Get embeddings, substring the barcode, only keep heavy chain, and remove NA values from dimensions and barcodes
        embedding_file = os.path.join("..","data",dataset,"VDJ",sample_id,"embeddings",input_source,f"embeddings_{model}.csv.gzip")
        embeddings_df = pd.read_csv(embedding_file, compression="gzip")
        embeddings_df['barcode'] = embeddings_df['barcode'].str[:-2] 
        embeddings_df = embeddings_df[embeddings_df['chain'] == "IGH"]
        embeddings_df = embeddings_df.dropna(subset=["dim_0"])
        embeddings_df = embeddings_df.dropna(subset=["barcode"])

        #Get VDJ information for a specific sample
        vdj_file = os.path.join("..","data",dataset,"vdj_evolike_combine.csv")
        vdj_df = pd.read_csv(vdj_file)
        vdj_df_s = vdj_df[vdj_df['sample_id'] == sample_id]

        #Match embeddings and vdj metadata
        df_merged = pd.merge(embeddings_df, vdj_df_s, on='barcode', how='inner')

        #Set sample names
        if dataset == "horns2020a__VDJ_RAW":
            df_merged["sample_id"] = df_merged["sample_id"].replace("Influenza.vac.11.12.human.S1","Individual1")
            if sample_id != "Influenza.vac.11.12.human.S1":
                continue
            #Only keep sample1          
            #i = df_merged[(df_merged["sample_id"] != 'Individual1')].index
            #df_merged.drop(i)
        if dataset == "Bruhn":
            df_merged["sample_id"] = df_merged["sample_id"].replace("Bruhn","Individual2")
        if dataset == "Kim":
            df_merged["sample_id"] = df_merged["sample_id"].replace("SRR17729703","Individual3").replace("SRR17729692","Individual4").replace("SRR17729726","Individual5")
        
        all_embeddings.append(df_merged)


## Samples combined
embeddings = pd.concat(all_embeddings, ignore_index=True)
df_merged = embeddings

#Get the sequences
if input_source == "full_VDJ":
    sequences = [str(gene) for gene in df_merged['VDJ_sequence_aa_trimmed']]
if input_source == "cdr3_only":
    sequences = [str(gene) for gene in df_merged['VDJ_cdr3_aa']]

#Get other metadata to plot
IgG_subtypes = ["IGHG1","IGHG2","IGHG2B","IGHG2C","IGHG3","IGHG4"]
IgA_subtypes = ["IGHA1","IGHA2", "IGHA"]
df_merged["c_gene"] = df_merged["c_gene"].replace(IgG_subtypes,"IgG").replace(IgA_subtypes,"IgA").replace("IGHM","IgM").replace("IGHD","IgD").replace("IGHE","IgE")
isotype = [str(gene) for gene in df_merged['c_gene']]

shm_count = [float(count) for count in df_merged["SHM_count"]]
sample = [str(s) for s in df_merged['sample_id']]
clonotype = [str(cl) for cl in df_merged['clonotype_id']]
df_merged['v_gene'] = df_merged['v_gene'].apply(lambda x: x.split('-')[0])
v_gene_family = [str(gene) for gene in df_merged['v_gene']]

#Get embedding columns and metadata columns
embedding_cols = [col for col in df_merged.columns if col.startswith('dim')]
metadata_cols = [col for col in df_merged.columns if col not in embedding_cols]

#Create an AnnData object
adata = anndata.AnnData(df_merged[embedding_cols])

#Add sequence and metadata information to .obs slot of AnnData object
adata.obs['seq'] = sequences
adata.obs['isotype'] = isotype
adata.obs["SHM_count"] = shm_count
adata.obs["sample"] = sample
adata.obs["clonotype"] = clonotype
adata.obs["v_gene_family"] = v_gene_family

#Remove duplicate observations
adata = adata[~adata.to_df().duplicated(), :]

#Construct sequence similarity network
evo.pp.neighbors(adata)
sc.tl.umap(adata)
basis = "umap"
if model  == "esm":
    evo.tl.velocity_graph(adata)
else:
    evo.tl.velocity_graph(adata, model_name = model)
if 'model' in adata.uns:
    del adata.uns['model']

#Embed network and velocities in two-dimensions
evo.tl.velocity_embedding(adata, basis = basis)

#Save the processed AnnData object to a file
output_file = os.path.join("..","data",f"adata_AllHumanSamples_{input_source}_{model}.h5ad")
adata.write(output_file)
print(f"Processed data saved to {output_file}")