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
The package uses `basilisk` to set up virtual environments to run the `python` code for the `cellOracle` worfklow.
However, to run the `scenicplus` analysis user has to set up the environment manually. Here are the shell 
commands for that:
```
conda create -n scenic_plus python=3.8
conda activate scenic_plus
git clone https://github.com/aertslab/scenicplus
cd scenicplus
pip install -e .
conda install -c anaconda requests
conda install -c numba numba
```

## To do:
- ArchR project for NGS4557 dataset is needed
- change the name of NGS4557 
- remove path to ArchR project in pruned_regulon_NGS4557.Rmd
- remove hard coded path to data matrices in scanpy code for cellOracle
- remove calls to load_all() preceded by setting environmental variable
- remove vignettes from .Rbuildignore
- plot LDA models to choose opitimal number of topics for NGS4557 dataset
- install the latest version of the scMultiSim package
- use atacseq_counts instead of atacseq_data
- visualization in the scenic plus workflow 
- add paths to data files with scores/rankings/tf motifs as report parameters  
- making gene names unique after concadetation of adata objects
- add ncpu parameter to find_topics
