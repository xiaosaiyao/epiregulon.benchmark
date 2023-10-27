# Benchmark of `epiregulon` package

The package allows for the replication of the benchmark analysis used for assessment of the performance of Epiregulon compared to other tools.
The user can run the benchmark report which is created from the scratch so that all the calculations need to be performed. Two datasets are used for the benchmark: reprogram-seq and AR. Each of them is analyzed by each tool and the code for that is organized into separate R markdown child files sourced by the main R markdown file. These files are stored in the `/inst/reports` folder. 

# Installation 
The package uses `basilisk` to set up virtual environments to run the `python` code for the `cellOracle` and `scenicplus` worfklow. Since `scenicplus` is not available through conda or pip, it should be copied into the package source folder. This step has to be done by the user due to space limitations in Roche GitLab.
```
git clone git@ssh.code.roche.com:grn/epiregulon/epiregulon_benchmark.git
cd epiregulon.benchmark/inst
git clone https://github.com/aertslab/scenicplus
```

After this step the package can be installed.

# Additional data
For the `GRaNIE` workflow the path to files with the transcriptpion factor motifs should be provided. This should be downloaded from https://www.embl.de/download/zaugg/diffTF/TFBS/TFBS_hg38_FIMO_HOCOMOCOv11.tar.gz and extracted. The storage folder will be used as argument to the function generating benchmark report.

## Usage

Run the benchmark for the **reprogram seq** data using `render_report_reprogram_seq` function. Provide the
following parameters:
- `treatments` the list with names corresponding to transcription factors being overexpressed and list elements to the strings present in the hash tags of the background group. If elements are `NULL` objects
then all hash tag groups apart from the target one are used as a background.
- `data_file_paths`: paths to the Cell Ranger output (feature_matrix.h5, fragment files); each path should correspond to a single sample;
defaults are set to Rosalind file location.
- `temp_dir`: path to the folder in which intermediates from the `scenicplus` analysis will be saved
- `work_dir`: path to the folder into which intermediate files from `scenicplus` analysis will be saved
- `n_cpu`: number of cores for parallel jobs
- `GRaNIE_tfbs`: path to the folder with `.bam` files containing information on motif data for the `GRaNIE` workflow
- `save_path`: path to which result of the benchmark will be saved, including plots

To run the benchmark for the **AR dataset** use the function `render_report_AR`. Provide the
following parameters:
- `path_to_ArchR_proj` path to the ArchR project created based on the AR dataset
- `data_file_paths`: paths to the Cell Ranger output (feature_matrix.h5, fragment files); each path should correspond to a single sample; defaults are set to Rosalind file location.
- `temp_dir`: path to the folder in which intermediates from the `scenicplus` analysis will be saved
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
