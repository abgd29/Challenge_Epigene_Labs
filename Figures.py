import os
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc 

###########################################################################################
### Importing the data to compare and creating the needed directories
###########################################################################################

solution = sc.read_h5ad("Output/GSE227828.h5ad")
exemple = sc.read_h5ad("Output_ex/GSE227828.h5ad")

if not os.path.exists("./Figures"):
    os.makedirs("./Figures")

if not os.path.exists("./Figures/cells"):
    os.makedirs("./Figures/cells")

if not os.path.exists("./Figures/genes"):
    os.makedirs("./Figures/genes")

###########################################################################################
### plotting the graphs for the cells layers 
###########################################################################################


## -------------------- Distribution of the log1p total counts distribution --------------------
fig, ax = plt.subplots(figsize=(8, 6))

sns.histplot(exemple.obs['log1p_total_counts'], kde=True, ax=ax, color='blue', label='Solution', alpha=0.5)
sns.histplot(solution.obs['log1p_total_counts'], kde=True, ax=ax, color='red', label='Exemple', alpha=0.5)

ax.set_title("Comparisons of the log1p of the total counts distribution")
ax.legend()
plt.savefig("./Figures/cells/total_counts")


fig, ax = plt.subplots(figsize=(8, 6))

sns.histplot(exemple.obs['log1p_total_counts'], kde=True, ax=ax, color='blue', label='Solution', alpha=0.5)
sns.histplot(solution.obs['log1p_total_counts'], kde=True, ax=ax, color='red', label='Exemple', alpha=0.5)

ax.set_title("Comparisons of the log1p of the total counts distribution")
ax.legend()
plt.savefig("./Figures/cells/total_counts")



### -------------------- Distribution of the number of genes expressed per cell --------------------
fig, ax = plt.subplots(figsize=(8, 6))

sns.histplot(exemple.obs['n_genes_by_counts'], kde=True, ax=ax, color='blue', label='Solution', alpha=0.5)
sns.histplot(solution.obs['n_genes_by_counts'], kde=True, ax=ax, color='red', label='Exemple', alpha=0.5)

ax.set_title("Comparisons of the number of genes expressed per cell distribution")
ax.legend()
plt.savefig("./Figures/cells/n_genes_by_counts")


### -------------------- Distribution of percentage of counts from the top 20 genes --------------------

fig, ax = plt.subplots(figsize=(8, 6))

sns.histplot(exemple.obs['pct_counts_in_top_20_genes'], kde=True, ax=ax, color='blue', label='Solution', alpha=0.5)
sns.histplot(solution.obs['pct_counts_in_top_20_genes'], kde=True, ax=ax, color='red', label='Exemple', alpha=0.5)

ax.set_title("Comparisons of percentage of counts from the top 20 genes distribution")
ax.legend()
plt.savefig("./Figures/cells/pct_counts_in_top_20_genes")


### -------------------- Distribution of log1p counts from ribosomic genes  --------------------

fig, ax = plt.subplots(figsize=(8, 6))

sns.histplot(exemple.obs['log1p_total_counts_ribo'], kde=True, ax=ax, color='blue', label='Solution', alpha=0.5)
sns.histplot(solution.obs['log1p_total_counts_ribo'], kde=True, ax=ax, color='red', label='Exemple', alpha=0.5)

ax.set_title("Comparisons of log1p counts from ribosomic genes distribution")
ax.legend()
plt.savefig("./Figures/cells/log1p_total_counts_ribo")


### -------------------- Distribution of log1p counts from genes related to hemoglobin  --------------------

fig, ax = plt.subplots(figsize=(8, 6))

sns.histplot(exemple.obs['log1p_total_counts_hb'], kde=True, ax=ax, color='blue', label='Solution', alpha=0.5)
sns.histplot(solution.obs['log1p_total_counts_hb'], kde=True, ax=ax, color='red', label='Exemple', alpha=0.5)

ax.set_title("Comparisons of log1p counts from genes related to hemoglobin distribution")
ax.legend()
plt.savefig("./Figures/cells/log1p_total_counts_hb")


### -------------------- Distribution of the nuber og genes expressed per cells  --------------------

fig, ax = plt.subplots(figsize=(8, 6))

sns.histplot(exemple.obs['n_genes'], kde=True, ax=ax, color='blue', label='Solution', alpha=0.5)
sns.histplot(solution.obs['n_genes'], kde=True, ax=ax, color='red', label='Exemple', alpha=0.5)

ax.set_title("Comparisons of the number of genes expressed per cells distribution")
ax.legend()
plt.savefig("./Figures/cells/n_genes")


### -------------------- Clustering of the cells  --------------------


fig, ax = plt.subplots(1, 2, figsize=(12, 6))  

sc.pl.pca(
    exemple,
    color="leiden_res_1",  
    ax=ax[0],  
    show=False,  
    title="Clusters - Exemple"
)


sc.pl.pca(
    solution,
    color = "leiden_res_1",
    ax=ax[1],  
    show=False,  
    title="Clusters - Solution"
)

plt.tight_layout() 
plt.savefig("./Figures/cells/leiden_res")

###########################################################################################
### plotting the graphs for the genes layers 
###########################################################################################

### -------------------- Log1p mean counts of each genes  --------------------

fig, ax = plt.subplots(figsize=(8, 6))

sns.histplot(exemple.var['log1p_mean_counts'], kde=True, ax=ax, color='blue', label='Solution', alpha=0.5)
sns.histplot(solution.var['log1p_mean_counts'], kde=True, ax=ax, color='red', label='Exemple', alpha=0.5)

ax.set_title("Comparisons of the log1p mean counts of each genes distribution")
ax.legend()
plt.savefig("./Figures/genes/log1p_mean_counts")

### -------------------- Mean of each genes  --------------------

fig, ax = plt.subplots(figsize=(8, 6))

sns.histplot(exemple.var['means'], kde=True, ax=ax, color='blue', label='Solution', alpha=0.5)
sns.histplot(solution.var['means'], kde=True, ax=ax, color='red', label='Exemple', alpha=0.5)

ax.set_title("Comparisons of the means of each genes distribution")
ax.legend()
plt.savefig("./Figures/genes/means")

### -------------------- Dispersion of each genes  --------------------

fig, ax = plt.subplots(figsize=(8, 6))

sns.histplot(exemple.var['dispersions'], kde=True, ax=ax, color='blue', label='Solution', alpha=0.5)
sns.histplot(solution.var['dispersions'], kde=True, ax=ax, color='red', label='Exemple', alpha=0.5)

ax.set_title("Comparisons of the dispersions of each genes distribution")
ax.legend()
plt.savefig("./Figures/genes/dispersions")

### -------------------- Normalized Dispersion of each genes  --------------------


fig, ax = plt.subplots(figsize=(8, 6))

sns.histplot(exemple.var['dispersions_norm'], kde=True, ax=ax, color='blue', label='Solution', alpha=0.5)
sns.histplot(solution.var['dispersions_norm'], kde=True, ax=ax, color='red', label='Exemple', alpha=0.5)

ax.set_title("Comparisons of the normalized dispersions of each genes distribution")
ax.legend()
plt.savefig("./Figures/genes/dispersions_norm")


