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
parser.add_argument("-d", "--dataset")

args = parser.parse_args()
model = args.model
input_source = args.input_source
dataset = args.dataset

all_embeddings = list()

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
    all_embeddings.append(df_merged)

    #Get the sequences
    if input_source == "full_VDJ":
        sequences = [str(gene) for gene in df_merged['VDJ_sequence_aa_trimmed']]
    if input_source == "cdr3_only":
        sequences = [str(gene) for gene in df_merged['VDJ_cdr3_aa']]

    #Get other metadata to plot
    IgG_subtypes = ["IGHG", "IGHG1","IGHG2","IGHG2B","IGHG2C","IGHG3","IGHG4"]
    IgA_subtypes = ["IGHA","IGHA1","IGHA2"]
    df_merged["c_gene"] = df_merged["c_gene"].replace(IgG_subtypes,"IgG").replace(IgA_subtypes,"IgA").replace("IGHM","IgM").replace("IGHD","IgD").replace("IGHE","IgE")
    isotype = [str(gene) for gene in df_merged['c_gene']]
    shm_count = [float(count) for count in df_merged["SHM_count"]]
    #sample = [str(s) for s in df_merged['sample_id_x']]
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
    #adata.obs["sample"] = sample
    adata.obs["clonotype"] = clonotype
    adata.obs["v_gene_family"] = v_gene_family

    #Remove duplicate observations
    adata = adata[~adata.to_df().duplicated(), :]

    #Construct sequence similarity network
    evo.pp.neighbors(adata)
    sc.tl.umap(adata)
    basis = "umap"
    print(adata.to_df().duplicated().sum())
    if model  == "esm":
        evo.tl.velocity_graph(adata)
    else:
        evo.tl.velocity_graph(adata, model_name = model)
    if 'model' in adata.uns:
        del adata.uns['model']

    #Embed network and velocities in two-dimensions
    evo.tl.velocity_embedding(adata, basis = basis)

    #Save the processed AnnData object to a file
    output_file = os.path.join("..","data",dataset,"evo-velocity",f"adata_{sample_id}_{input_source}_{model}.h5ad")
    adata.write(output_file)
    print(f"Processed data saved to {output_file}")


## Samples combined
embeddings = pd.concat(all_embeddings, ignore_index=True)
df_merged = embeddings

#Get the sequences
if input_source == "full_VDJ":
    sequences = [str(gene) for gene in df_merged['VDJ_sequence_aa_trimmed']]
if input_source == "cdr3_only":
    sequences = [str(gene) for gene in df_merged['VDJ_cdr3_aa']]

#Get other metadata to plot
#Set isotypes
IgG_subtypes = ["IGHG1","IGHG2","IGHG2B","IGHG2C","IGHG3","IGHG4"]
IgA_subtypes = ["IGHA1","IGHA2", "IGHA"]
df_merged["c_gene"] = df_merged["c_gene"].replace(IgG_subtypes,"IgG").replace(IgA_subtypes,"IgA").replace("IGHM","IgM").replace("IGHD","IgD").replace("IGHE","IgE")
df_merged["sample_id"] = df_merged["sample_id"].replace("S1","Mouse1").replace("S2","Mouse2").replace("S3","Mouse3").replace("S4","Mouse4").replace("S5","Mouse5")
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
output_file = os.path.join("..","data",dataset,"evo-velocity",f"adata_AllSamples_{input_source}_{model}.h5ad")
adata.write(output_file)
print(f"Processed data saved to {output_file}")
