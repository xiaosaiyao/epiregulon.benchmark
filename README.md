# Benchmark of `epiregulon` package

The package allows for the replication of the benchmark analysis used for assessment of the performance of Epiregulon against other tools.
The user can run the benchmark report which is created from the scratch so that all the calculations need to be performed. Two datasets are used for the benchmark: reprogram-seq and AR. Each of them is analyzed by each benchmarked tool and the code is organized into separate R markdown child files sourced by the main file. These files are stored in the `/inst/reports` folder. 

# Preparation
If you want to run `scenicplus` workflow you have to download motif databases which are use by `pycistarget` as a part of the data pre-processing. By default they are downloaded into the package installation directory by running `download_databases` function from https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/. 
The same is true for the HOCOMOCO motif database used by `GRaNIE`. It can be downloaded and extracted using `download_HOCOMOCO` function. The file url is https://www.embl.de/download/zaugg/diffTF/TFBS/TFBS_hg38_FIMO_HOCOMOCOv11.tar.gz.

# Usage

Run the benchmark for the **reprogram seq** data using `render_report_reprogram_seq` function. Provide the
following parameters:
- `treatments` the list with names corresponding to transcription factors being overexpressed and list elements to the strings present in the hash tags of the background group. If elements are `NULL` objects
then all hash tag groups apart from the target one are used as a background.
- `data_file_paths`: paths to the Cell Ranger output (feature_matrix.h5, fragment files); each path should correspond to a single sample;
defaults are set to Rosalind file location.
- `temp_dir`: path to the folder in which intermediates from the `scenicplus` analysis will be saved
- `work_dir`: path to the folder into which intermediate files from `scenicplus` analysis will be saved
- `n_cpu`: number of cores for parallel jobs
- `save_path`: path to which result of the benchmark will be saved, including plots

To run the benchmark for the **AR dataset** use the function `render_report_AR`. Provide the
following parameters:
- `path_to_ArchR_proj` path to the ArchR project created based on the AR dataset
- `data_file_paths`: paths to the Cell Ranger output (feature_matrix.h5, fragment files); each path should correspond to a single sample; defaults are set to Rosalind file location.
- `temp_dir`: path to the folder in which intermediates from the `scenicplus` analysis will be saved
- `work_dir`: path to the folder into which intermediate files from `scenicplus` analysis will be saved
- `n_cpu`: number of cores for parallel jobs
- `save_path`: path to which result of the benchmark will be saved, including plots

# To do:
- ArchR project for NGS4557 dataset is needed
- remove path to ArchR project in pruned_regulon_NGS4557.Rmd
- remove vignettes from .Rbuildignore
- use atacseq_counts instead of atacseq_data
- making gene names unique after concatenation of adata objects
