import scanpy as sc 
import os
import anndata as ad

def get_annotated_data(barcodes_list, HTO_list, data_paths = [""], sample_names = None, n_top_genes=2000):

    # set directory with data
    adata_list = map(lambda x: sc.read_10x_h5(os.path.join(dataPath_1, "filtered_feature_bc_matrix.h5")), data_paths)
    adata_list = list(map(lambda x: x.var_names_make_unique(), adata_list))
    if samples_name:
        for(i in range(length(adata_list))):
            adata_list[i].obs_names = list(map(lambda x: sample_names[i]+"#"+x, adata_1.obs_names))
    if len(data_paths)>1:
        adatas = dict(zip(sample_names, adata_list))
        adata = ad.concat(adatas, label = "dataset")
    if False in list(map(lambda x: x in adata.obs_names, barcodes_list)):
        raise Exception("All barcodes should be present in the gene expression data")
    
    adata = adata[barcodes_list,]
    adata.obs["HTO"] = HTO_list
    
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
