import scanpy as sc 
import os
import anndata as ad

def get_annotated_data(barcode_tab, data_paths = [""], sample_names = None, 
    n_top_genes=2000, mode = "cellOracle"):

    # set directory with data
    adata_list = list(map(lambda x: sc.read_10x_h5(os.path.join(x, "filtered_feature_bc_matrix.h5")), data_paths))
    for i in range(len(sample_names)):
        adata_list[i].obs['sample_id'] = sample_names[i]
        adata_list[i].obs_names = list(map(lambda x: x.replace("-", "."), adata_list[i].obs_names))
        sample_barcodes = barcode_tab.loc[barcode_tab['sample_id'] == sample_names[i]]["barcode"]
        if False in list(map(lambda x: x in adata_list[i].obs_names, sample_barcodes)):
            raise Exception("All barcodes should be present in the gene expression data")
        adata_list[i] = adata_list[i][sample_barcodes,]
        adata_list[i].obs["HTO"] = barcode_tab.loc[barcode_tab['sample_id'] == sample_names[i]]["hash_assignment"].values
        adata_list[i].obs_names = list(map(lambda x: x[0]+"___"+x[1], zip(adata_list[i].obs_names, adata_list[i].obs['sample_id'])))
        adata_list[i].var_names_make_unique()
    adatas = dict(zip(sample_names, adata_list))
    adata = ad.concat(adatas, label = "dataset")
    adata.var_names_make_unique()
    if mode == "cellOracle":
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
    
    else:
        # filtering skipped
        # adata = adata[adata.obs.n_genes_by_counts < n_counts_filter, :]
        # adata = adata[adata.obs.pct_counts_mt < mito_filter, :]
        
        # skip normalization since it's only used for visualization purposes
        # adata.raw = adata
        # sc.pp.normalize_total(adata, target_sum=1e4)
        # sc.pp.log1p(adata)
        # sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        # adata = adata[:, adata.var.highly_variable]
        # sc.pp.scale(adata, max_value=10)
        # sc.tl.pca(adata, svd_solver='arpack')
        # sc.pl.pca_variance_ratio(adata, log=True)
        # sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
        # sc.tl.umap(adata)
        # sc.pl.umap(adata, color = "HTO")
        
        return adata
    
