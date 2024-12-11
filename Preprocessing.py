import os
import pandas as pd
from scipy.io import mmread
import scanpy as sc #The scanpy.pp library contains all the functions needed to preprocess data 
import numpy as np



############################################################################################################
#### Downloading data and creating the AnnData object
############################################################################################################

# The two samples considered in the study
echant1 = "GSM7109167"
echant2 = "GSM7109168"

# Creating the AnnData object for the first sample
expression_matrix1 = (
    mmread("./Input/" + echant1 + "_PT1_matrix.mtx").tocsc()
).transpose()
cell_barcodes1 = pd.read_csv(
    "./Input/" + echant1 + "_PT1_barcodes.tsv", header=None, index_col=None
)
gene_features1 = pd.read_csv(
    "./Input/" + echant1 + "_PT1_features.tsv", header=None, index_col=None, sep="\t"
)


obs1 = pd.DataFrame(index=(echant1 + "_" + cell_barcodes1[0].values))
var1 = pd.DataFrame(index=(gene_features1[1].values))
adata1 = sc.AnnData(X=expression_matrix1, obs=obs1, var=var1)
adata1.obs["dataset"] = "GSE227828"
adata1.obs["sample_id"] = echant1
adata1.var_names_make_unique()

# Creating the AnnData object for the second sample
expression_matrix2 = (
    mmread("./Input/" + echant2 + "_PT5_matrix.mtx").tocsc()
).transpose()
cell_barcodes2 = pd.read_csv(
    "./Input/" + echant2 + "_PT5_barcodes.tsv", header=None, index_col=None
)
gene_features2 = pd.read_csv(
    "./Input/" + echant2 + "_PT5_features.tsv", header=None, index_col=None, sep="\t"
)

obs2 = pd.DataFrame(index=(echant2 + "_" + cell_barcodes2[0].values))
var2 = pd.DataFrame(index=(gene_features2[1].values))
adata2 = sc.AnnData(X=expression_matrix2, obs=obs2, var=var2)
adata2.obs["dataset"] = "GSE227828"
adata2.obs["sample_id"] = echant2
adata2.var_names_make_unique()

# Concatenation of the two AnnData objects
adata = sc.concat([adata1, adata2], join="outer")

############################################################################################################
#### Adding the genes caracteristics (obs) and filtering genes
############################################################################################################

#Getting the gene ids and feature types from the gene_features matrix
adata.var["gene_ids"] = gene_features1[0].values
adata.var["feature_types"] = gene_features1[2].values

#Using gene names to determine if they are mitochondrial, ribosomal, or associated with hemoglobin
adata.var["mt"] = adata.var.index.str.startswith("MT-")
adata.var["ribo"] = adata.var.index.str.startswith(
    "RPL"
) | adata.var.index.str.startswith("RPS")
adata.var["hb"] = adata.var.index.str.contains("^(HBA|HBB|HBE|HBS|HBG|HBQ)", regex=True)


############################################################################################################
#### Computing the main metrics
############################################################################################################

# Computing the caracteristics for all type of genes, mitochondrial genes, ribosomal genes and hemoglobin genes
sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=[20])
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, inplace=True)
sc.pp.calculate_qc_metrics(adata, qc_vars=["ribo"], percent_top=None, inplace=True)
sc.pp.calculate_qc_metrics(adata, qc_vars=["hb"], percent_top=None, inplace=True)

############################################################################################################
#### Filtering the cells and genes
############################################################################################################

#We use the metrics previously computed to determine the outlier cells
adata.obs["mt_outlier"] = adata.obs["total_counts_mt"] > 0 #If a cell express mitochondrial genes

adata.obs["outlier"] = (
    (adata.obs["total_counts"] < 605 ) #At least 605 reads expressed per cells
    | (adata.obs["n_genes_by_counts"] < 405) # At least 405 genes expressed per cells
    | (adata.obs["n_genes_by_counts"] > 3801) # At most 3801 genes expressed per cells
    | (adata.obs["total_counts_hb"] > 1000)  #No more than 1000 genes related to hemoglobin expressed
    | (adata.obs["pct_counts_hb"] > 12) #No more than 12% of the genes expressed should be related to hemoglobin
    | (adata.obs["pct_counts_in_top_20_genes"] > 27.5) #The top 20% genes should not represent more than 27.5% of total expression
)


#Keep only the cells that are not outliers
adata = adata[adata.obs["outlier"] == False].copy()
adata = adata[adata.obs["mt_outlier"] == False].copy()

#Keep only the genes that are not mitochondrial 
adata = adata[:,adata.var["mt"] == False].copy()


# Keep the cells and genes related to enough genes and cells
sc.pp.filter_genes(adata, min_cells=20)
sc.pp.filter_cells(adata, min_genes=405)


############################################################################################################
####  Computing the cells specific metrics
############################################################################################################


#Computing the size factors
total_counts = adata.obs["log1p_total_counts"]
adata.obs["size_factors"] = total_counts / total_counts.mean()

#Dividing the cells into clusters
#       Agregating the cells in clusters thanks to neibourgh joining algorithm
sc.pp.neighbors(adata, n_neighbors=16,use_rep='X')
adata.X = adata.X.astype('float32')

#   Attributing a cluster to each cells
sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=10)
sc.tl.leiden(adata, resolution=1.0, key_added="leiden_res_1")


############################################################################################################
####  Computing the genes specific metrics
############################################################################################################

#We use a scanpy fonction to get every informations we need
adata.layers['log1p'] = np.log1p(adata.X) #We wish to work with the normalized data
a=sc.pp.highly_variable_genes(
    adata,
    layer = "log1p",
    n_top_genes = 4000, #Number of genes to be considered Highly variables
    flavor="cell_ranger",  # We use the cell_ranger method to be able to compute both dispersions and normalised disperions
    inplace =False 
)

adata.var['means'] =a["means"]
adata.var['dispersions'] =a["dispersions"]
adata.var['dispersions_norm'] = abs(a["dispersions_norm"])
adata.var['highly_variable'] = a['highly_variable']

############################################################################################################
#### Formatting the metrics
############################################################################################################

#Formating the cells metrics
adata.obs["total_counts"] = adata.obs["total_counts"].apply(lambda x: f"{x:.1f}")
adata.obs["log1p_total_counts"] = adata.obs["log1p_total_counts"].round(6)
adata.obs["total_counts_mt"] = adata.obs["total_counts_mt"].apply(lambda x: f"{x:.1f}")
adata.obs["total_counts_ribo"] = adata.obs["total_counts_ribo"].apply(
    lambda x: f"{x:.1f}"
)

#Formating the genes metrics
adata.obs["log1p_total_counts_ribo"] = adata.obs["log1p_total_counts_ribo"].round(7)
adata.obs["pct_counts_ribo"] = adata.obs["pct_counts_ribo"].round(7)

adata.obs["total_counts_hb"] = adata.obs["total_counts_hb"].apply(lambda x: f"{x:.1f}")
adata.obs["log1p_total_counts_hb"] = adata.obs["log1p_total_counts_hb"].round(6)
adata.obs["pct_counts_hb"] = adata.obs["pct_counts_hb"].round(6)

adata.var["mean_counts"] = adata.var["mean_counts"].round(9)
adata.var["log1p_mean_counts"] = adata.var["log1p_mean_counts"].round(9)
adata.var["total_counts"] = adata.var["total_counts"].apply(lambda x: f"{x:.1f}")
adata.var["log1p_total_counts"] = adata.var["log1p_total_counts"].round(7)
adata.var["means"] = adata.var["means"].round(9)
adata.var["dispersions"] = adata.var["dispersions"].round(7)
adata.var["dispersions_norm"] = adata.var["dispersions_norm"].round(7)


############################################################################################################
#### Export outpouts
############################################################################################################

if not os.path.exists("./Output"):
    os.makedirs("./Output")
adata.write("./Output/GSE227828.h5ad")



