# Benchmark of epiregulon package

## Usage
To generate report showing the ROC curves for the different weight methods use the function
`render_report_reprogram_seq` with several arguments:
- `regulon` The regulon object with `tf`, `target` and `idxATAC` columns. Weights, if computed,
will be disregarded. 
- `treatments` Named list with names being transcription factors included in the `hasg_assignment` column in
the `GeneExpressionMatrix`. List elements are the vectors with groups which should be used as a control.

```
render_report_reprogram_seq(data_file_paths = "/gne/data/lab-shares/xie-lab/Sequencing_Data/2022/mapping/20220124_ReprogramSeq_Multiome/JT65_67/outs/",
                            work_dir = "/gstore/scratch/u/wlodarct",
                            temp_folder = "/gstore/scratch/u/wlodarct/",
                            conda_exe = "/gstore/home/wlodarct/anaconda3/bin/conda",
                            virtual_env = "/gstore/home/wlodarct/anaconda3/envs/scenic_plus",
                            n_cpu = 6,
                            output_file = file.path(getwd(), "Reprogram-seq_report_2.html"))
```


## Benchmarking against other tools
The package uses `basilisk` to set up virtual environments to run the `python` code for the `cellOracle` and `scenicplus` worfklow.

Run the benchmark for the reprogram seq data using `render_report_reprogram_seq` function. Provide the
following parameters:
- `treatments` the list with names corresponding to transcription factors being overexpressed and list elements to the strings present in the hash tags of the background group. If elements are `NULL` objects
then all hash tag groups apart from the target one are used as a background.
- `data_file_paths`: paths to the Cell Ranger output (feature_matrix.h5, fragment files); each path should correspond to a single sample
- `temp_dir`: path to the folder into which temporary from `scenicplus` analysis will be saved
- `work_dir`: path to the folder into which intermediate files from `scenicplus` analysis will be saved
- `n_cpu`: number of cores for parallel jobs
- `GRaNIE_tfbs`: path to the folder with `.bam` files containing information on motif data for the `GRaNIE` workflow
- `save_path`: path to which result of the benchmark will be saved, including plots

## To do:
- ArchR project for NGS4557 dataset is needed
- remove path to ArchR project in pruned_regulon_NGS4557.Rmd
- remove vignettes from .Rbuildignore
- use atacseq_counts instead of atacseq_data
- add paths to data files with scores/rankings/tf motifs as report parameters  
- making gene names unique after concatenation of adata objects
