import scanpy as sc 
import os
import anndata as ad

def get_annotated_data_2(barcodes_list, HTO_list, data_path = ""):

    # set directory with data
    dataPath_1 = "/gstore/data/genomics/congee_rest_runs/62e2c81cd6d7e0bd49c579dd/SAM24418230/croo_output/"
    dataPath_2 = "/gstore/data/genomics/congee_rest_runs/62e2c81cd6d7e0bd49c579dd/SAM24418231/croo_output/"

    adata_1 = sc.read_10x_h5(os.path.join(dataPath_1, "filtered_feature_bc_matrix.h5"))
    adata_1.var_names_make_unique()
    adata_2 = sc.read_10x_h5(os.path.join(dataPath_2, "filtered_feature_bc_matrix.h5"))
    adata_2.var_names_make_unique()

    adata_1.obs_names = list(map(lambda x: "SAM24418230#"+x, adata_1.obs_names))
    adata_2.obs_names = list(map(lambda x: "SAM24418231#"+x, adata_2.obs_names))

    adatas = {"SAM24418230" : adata_1, "SAM24418231": adata_2}
    adata = ad.concat(adatas, label = "dataset")

    barcodes_updated = list(set(barcodes_list).intersection(set(adata.obs_names)))
    adata = adata[barcodes_updated,]
    adata.obs["HTO"] = HTO_list

    # Use 5000 instead of 2000 genes since otherwise there an issue caused by the lack of the convergence raised by the sc.tl.diffmap
    # function or a waring about many disconnected components in the transformation matrix

    # Select top 5000 highly-variable genes
    filter_result = sc.pp.filter_genes_dispersion(adata.X,
                                                  flavor='cell_ranger',
                                                  n_top_genes=5000,
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
