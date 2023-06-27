# preparation of the data to be used in scenic plus

# data retreived using python
# adata = scanpy.read_10x_h5("/gne/data/lab-shares/xie-lab/Sequencing_Data/2022/mapping/20220124_ReprogramSeq_Multiome/JT65_67/outs/filtered_feature_bc_matrix.h5")
# adata2.var_names_make_unique()
# df = pd.DataFrame(adata2.var)
# df.to_csv("adata_vars_reprogram_seq.csv")
data("adata_vars_reprogram_seq")
write.csv(adata_vars_reprogram_seq, "adata_vars_reprogram_seq.csv")
