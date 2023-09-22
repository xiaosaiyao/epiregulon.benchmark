import os
import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
from pycisTopic.cistopic_class import *
from pycisTopic.lda_models import *
from pycisTopic.clust_vis import *
from pycisTopic.topic_binarization import *
from pycisTopic.diff_features import *
from scenicplus.wrappers.run_pycistarget import run_pycistarget

def find_topics(adata, sample_names, paths_to_fragments, work_dir, tmp_dir, paths_to_peak_matrix, n_cpu, group_variable):
    scRNA_bc = adata.obs_names
    cell_data = adata.obs
    cell_data[group_variable] = cell_data[group_variable].astype(str) # set data type of the celltype column to str, otherwise the export_pseudobulk function will complain.
    #cell_data.index = list(map(lambda x: x[0]+"___"+x[1], zip(cell_data.index, cell_data['sample_id'])))
    fragments_dict = dict(zip(sample_names, paths_to_fragments))
    matrices_dict = dict(zip(sample_names, paths_to_peak_matrix))
    cistopic_obj_list = [create_cistopic_object_from_matrix_file(fragment_matrix_file = matrices_dict[key],
                                                      path_to_fragments=fragments_dict[key],
                                                      project = key) for key in fragments_dict.keys()]
    cistopic_obj = merge(cistopic_obj_list)
    cistopic_obj.add_cell_data(cell_data)
    models=run_cgs_models(cistopic_obj,
                n_topics=[2,4,8,15,20,25,32,38,48],
                n_cpu=8,
                n_iter=500,
                random_state=555,
                alpha=50,
                alpha_by_topic=True,
                eta=0.1,
                eta_by_topic=False,
                save_path=None,
                _temp_dir = '/gstore/scratch/u/wlodarct/temp/scenic/')
                    
    model = evaluate_models(models,
                       select_model=25,
                       return_model=True,
                       metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                       plot_metrics=False)

    cistopic_obj.add_LDA_model(model)
    run_umap(cistopic_obj, target  = 'cell', scale=True)
    region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')
    region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop = 3000)
    imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
    normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
    variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)
    markers_dict = find_diff_features(cistopic_obj, imputed_acc_obj, variable='HTO', var_features=variable_regions)
    region_sets = {}
    region_sets['topics_otsu'] = {}
    region_sets['topics_top_3'] = {}
    region_sets['DARs'] = {}
    for topic in region_bin_topics_otsu.keys():
        regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
        region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
    for topic in region_bin_topics_top3k.keys():
        regions = region_bin_topics_top3k[topic].index[region_bin_topics_top3k[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
        region_sets['topics_top_3'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
    for DAR in markers_dict.keys():
        # account for clusters without differetnially accessible regions
        if markers_dict[DAR].size == 0: continue
        regions = markers_dict[DAR].index[markers_dict[DAR].index.str.startswith('chr')] #only keep regions on known chromosomes
        region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))
    # use databased downloaded from https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/
    rankings_db = '/gstore/scratch/u/wlodarct/GREDXFER-4710/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather'
    scores_db = '/gstore/scratch/u/wlodarct/GREDXFER-4710/hg38_screen_v10_clust.regions_vs_motifs.scores.feather'
    motif_annotation = '/gstore/scratch/u/wlodarct/GREDXFER-4710/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl' #motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
    if not os.path.exists(os.path.join(work_dir, 'motifs')):
        os.makedirs(os.path.join(work_dir, 'motifs'))
    run_pycistarget(
        region_sets = region_sets,
        species = 'homo_sapiens',
        save_path = os.path.join(work_dir, 'motifs'),
        ctx_db_path = rankings_db,
        dem_db_path = scores_db,
        path_to_motif_annotations = motif_annotation,
        run_without_promoters = True,
        n_cpu = n_cpu,
        _temp_dir = os.path.join(tmp_dir, 'ray_spill'),
        annotation_version = 'v10nr_clust'
    )
    return cistopic_obj
