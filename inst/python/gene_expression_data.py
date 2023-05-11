import scanpy as sc 
import os
import anndata as ad

def get_annotated_data(barcodes_list, HTO_list, data_path = ""):

    # indicate directory with the data
    dataPath = "/gne/data/lab-shares/xie-lab/Sequencing_Data/2022/mapping/20220124_ReprogramSeq_Multiome/JT65_67/outs/"
    adata = sc.read_10x_h5(os.path.join(dataPath, "filtered_feature_bc_matrix.h5"))
    adata.var_names_make_unique()

    barcodes_updated = list(set(barcodes_list).intersection(set(adata.obs_names)))
    adata = adata[barcodes_updated,]
    adata.obs["HTO"] = HTO_list
    
    # remove genes not being expressed in any cell; otherwise calling sc.pp.filter_genes_dispersion might raise error
    non_zero_genes = adata.X.sum(axis=0)>0
    adata = adata[:,non_zero_genes]

    # Select top 2000 highly-variable genes
    filter_result = sc.pp.filter_genes_dispersion(adata.X,
                                                  flavor='cell_ranger',
                                                  n_top_genes=2000,
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
