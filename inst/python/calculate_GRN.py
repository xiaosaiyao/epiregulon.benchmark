import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import celloracle as co


def calculate_GRN(adata, base_GRN):
	# Instantiate Oracle object
	oracle = co.Oracle()

	# Use the unscaled mRNA count for the input of Oracle object.
	adata.X = adata.layers["raw_count"].copy()

	# Instantiate Oracle object.
	oracle.import_anndata_as_raw_count(adata=adata,
	                                   cluster_column_name="louvain",
	                                   embedding_name="X_umap")

	# Load TF info dataframe 
	oracle.import_TF_data(TF_info_matrix=base_GRN)
	 
	# Perform PCA
	oracle.perform_PCA()

	# Select important PCs
	n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
	n_comps = min(n_comps, 50)
	
	n_cell = oracle.adata.shape[0]
	
	k = int(0.025*n_cell)

	oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
	                      b_maxl=k*4, n_jobs=4)

	links = oracle.get_links(cluster_name_for_GRN_unit="louvain", alpha=10,
	                         verbose_level=10)


	df = pd.DataFrame()
	for cluster in links.links_dict.keys():
	    df_new = links.links_dict[cluster]
	    df_new['cluster'] = cluster
	    df = pd.concat([df, df_new])

	return df

