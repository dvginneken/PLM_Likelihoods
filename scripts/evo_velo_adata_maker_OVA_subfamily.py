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
parser.add_argument("-v", "--v_gene", help="V-gene subfamily. Available: IGHV3-1, IGHV9-3, IGHV5-6")

args = parser.parse_args()
model = args.model
subgene = args.v_gene
input_source = "v_gene_only"
dataset = "OVA_V7"

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



## Samples combined
embeddings = pd.concat(all_embeddings, ignore_index=True)
df_merged = embeddings

#Get the v-gene sequences
df_merged['full_sequence'] = df_merged["VDJ_fwr1_aa"] + df_merged["VDJ_cdr1_aa"] + df_merged["VDJ_fwr2_aa"] + df_merged["VDJ_cdr2_aa"] + df_merged["VDJ_fwr3_aa"]

#Set isotypes
IgG_subtypes = ["IGHG1","IGHG2","IGHG2B","IGHG2C","IGHG3","IGHG4"]
IgA_subtypes = ["IGHA1","IGHA2", "IGHA"]
df_merged["c_gene"] = df_merged["c_gene"].replace(IgG_subtypes,"IgG").replace(IgA_subtypes,"IgA").replace("IGHM","IgM").replace("IGHD","IgD").replace("IGHE","IgE")
df_merged["sample_id"] = df_merged["sample_id"].replace("S1","Mouse1").replace("S2","Mouse2").replace("S3","Mouse3").replace("S4","Mouse4").replace("S5","Mouse5")

#Only keep specific v-gene sub family
df_merged = df_merged[df_merged['v_gene'] == subgene]

#Subset columns
embedding_cols = [col for col in df_merged.columns if col.startswith('dim')]
sequence_embeddings = df_merged[embedding_cols]
sequence_metadata = df_merged[["barcode","sample_id","full_sequence","c_gene","v_gene","SHM_count"]]
sequence_df = pd.concat([sequence_metadata, sequence_embeddings], axis=1)

## Germlines
germline_file = os.path.join("..","data","IMGT_germline_embeddings",dataset,f"embeddings_{model}.csv")
germline_df = pd.read_csv(germline_file)
germline_df = germline_df[germline_df['v_gene'] == subgene]
germline_df["sample_id"] = "germline"
germline_df["c_gene"] = "germline"
germline_df["SHM_count"] = 0
germline_df["barcode"] = "germline-" + germline_df["v_gene"]

#Subset columns
embedding_cols = [col for col in germline_df.columns if col.startswith('dim')]
germline_embeddings = germline_df[embedding_cols]
germline_metadata = germline_df[["barcode","sample_id","full_sequence","c_gene","v_gene","SHM_count"]]
germline_df = pd.concat([germline_metadata, germline_embeddings], axis=1)

# Combine germline and sequence dataframe
df_merged = pd.concat([germline_df, sequence_df], axis=0, ignore_index=True)

#Get metadata
sequences = [str(gene) for gene in df_merged['full_sequence']]
isotype = [str(gene) for gene in df_merged['c_gene']]
shm_count = [float(count) for count in df_merged["SHM_count"]]
sample = [str(s) for s in df_merged['sample_id']]
v_gene =  [str(s) for s in df_merged['v_gene']]

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
adata.obs["v_gene"] = v_gene

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
output_file = os.path.join("..","data",dataset,"evo-velocity",f"adata_{subgene}_IncludingGermline_{input_source}_{model}.h5ad")
adata.write(output_file)
print(f"Processed data saved to {output_file}")