import scanpy as sc 
import os
import anndata as ad

def get_annotated_data(barcode_tab, data_paths = [""], sample_names = None, 
    n_top_genes=2000, group_variable=""):

    # set directory with data
    adata_list = list(map(lambda x: sc.read_10x_h5(os.path.join(x, "filtered_feature_bc_matrix.h5")), data_paths))
    for i in range(len(sample_names)):
        adata_list[i].obs['sample_id'] = sample_names[i]
        adata_list[i].obs_names = list(map(lambda x: x.replace("-", "."), adata_list[i].obs_names))
        sample_barcodes = barcode_tab.loc[barcode_tab['sample_id'] == sample_names[i]]["barcode"]
        if False in list(map(lambda x: x in adata_list[i].obs_names, sample_barcodes)):
            raise Exception("All barcodes should be present in the gene expression data")
        adata_list[i] = adata_list[i][sample_barcodes,]
        adata_list[i].obs[group_variable] = barcode_tab.loc[barcode_tab['sample_id'] == sample_names[i]][group_variable].values
        adata_list[i].obs_names = list(map(lambda x: x[0]+"___"+x[1], zip(adata_list[i].obs_names, adata_list[i].obs['sample_id'])))
        adata_list[i].var_names_make_unique()
    adatas = dict(zip(sample_names, adata_list))
    adata = ad.concat(adatas, label = "dataset")
    adata.var_names_make_unique()
    # Only consider genes with more than 1 count
    sc.pp.filter_genes(adata, min_counts=1)
    
    # Normalize gene expression matrix with total UMI count per cell
    sc.pp.normalize_per_cell(adata, key_n_counts='n_counts_all')

    # Select top highly-variable genes
    filter_result = sc.pp.filter_genes_dispersion(adata.X,
                                                  flavor='cell_ranger',
                                                  n_top_genes=n_top_genes,
                                                  log=False)

    # Subset the genes
    adata = adata[:, filter_result.gene_subset]

    # Renormalize after filtering
    sc.pp.normalize_per_cell(adata)

    # Data normalization
    # save raw counts
    adata.raw = adata 
    adata.layers["raw_count"] = adata.raw.X.copy() 

    # Log transformation and scaling
    sc.pp.log1p(adata)
    sc.pp.scale(adata)

    # dimensionality reduction
    # PCA
    sc.tl.pca(adata, svd_solver='arpack')

    # Diffusion map
    sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)

    sc.tl.diffmap(adata, n_comps = 15)
    # Calculate neihbors again based on diffusionmap
    sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_diffmap')

    # cell clustering
    sc.tl.louvain(adata, resolution=0.8)

    #UMAP
    sc.tl.umap(adata)

    return adata
 
