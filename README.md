# Benchmark of `epiregulon` package

## Usage
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

To run the benchmark for the AR dataset use the function `render_report_AR`. Provide the
following parameters:
- `path_to_ArchR_proj` path to the ArchR project created based on the AR dataset
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
