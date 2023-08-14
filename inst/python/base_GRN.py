import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


import os, sys, shutil, importlib, glob
from tqdm.notebook import tqdm
from celloracle import motif_analysis as ma
import celloracle as co
from celloracle.utility import save_as_pickled_object

def prepare_base_GRN(cicero_connections, peaks, association_cutoff=0.8):
    #print(peaks)
    peaks = peaks.x.values
    tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome="hg38")
    integrated = ma.integrate_tss_peak_with_cicero(tss_peak=tss_annotated, cicero_connections=cicero_connections)
    peak = integrated[integrated.coaccess >= association_cutoff]
    peak = peak[["peak_id", "gene_short_name"]].reset_index(drop=True)
    ref_genome = "hg38"
    genome_installation = ma.is_genome_installed(ref_genome=ref_genome)
    print(ref_genome, "installation: ", genome_installation)
    if not genome_installation:
        import genomepy
        genomepy.install_genome(name=ref_genome, provider="UCSC")

    def decompose_chrstr(peak_str):
        """
        Args:
            peak_str (str): peak_str. e.g. 'chr1_3094484_3095479'

        Returns:
            tuple: chromosome name, start position, end position
        """

        *chr_, start, end = peak_str.split("_")
        chr_ = "_".join(chr_)
        return chr_, start, end

    from genomepy import Genome

    def check_peak_format(peaks_df, ref_genome):
        """
        Check peak format.
         (1) Check chromosome name.
         (2) Check peak size (length) and remove sort DNA sequences (<5bp)

        """

        df = peaks_df.copy()

        n_peaks_before = df.shape[0]

        # Decompose peaks and make df
        decomposed = [decompose_chrstr(peak_str) for peak_str in df["peak_id"]]
        df_decomposed = pd.DataFrame(np.array(decomposed), index=peaks_df.index)
        df_decomposed.columns = ["chr", "start", "end"]
        df_decomposed["start"] = df_decomposed["start"].astype(int)
        df_decomposed["end"] = df_decomposed["end"].astype(int)

        # Load genome data
        genome_data = Genome(ref_genome)
        all_chr_list = list(genome_data.keys())


        # DNA length check
        lengths = np.abs(df_decomposed["end"] - df_decomposed["start"])


        # Filter peaks with invalid chromosome name
        n_threshold = 5
        df = df[(lengths >= n_threshold) & df_decomposed.chr.isin(all_chr_list)]

        # DNA length check
        lengths = np.abs(df_decomposed["end"] - df_decomposed["start"])

        # Data counting
        n_invalid_length = len(lengths[lengths < n_threshold])
        n_peaks_invalid_chr = n_peaks_before - df_decomposed.chr.isin(all_chr_list).sum()
        n_peaks_after = df.shape[0]


        return df


    peaks = check_peak_format(peak, ref_genome)

    # Instantiate TFinfo object
    tfi = ma.TFinfo(peak_data_frame=peaks,
                    ref_genome=ref_genome)

    # Scan motifs. !!CAUTION!! This step may take several hours if you have many peaks!
    tfi.scan(fpr=0.02,
             motifs=None,  # If you enter None, default motifs will be loaded.
             verbose=True)

    # Save tfinfo object
    tfi.to_hdf5(file_path="/gstore/project/lineage/tomek/method-benchmark/analysis/output_data/reprogram-seq/cellOracle/test1.celloracle.tfinfo")

    # Check motif scan results
    tfi.scanned_df.head()

    # Reset filtering
    tfi.reset_filtering()

    # Do filtering
    tfi.filter_motifs_by_score(threshold=10)

    # Format post-filtering results.
    tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)

    return tfi.to_dataframe()

