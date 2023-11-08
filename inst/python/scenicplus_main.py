import dill
import scanpy as sc
import os
import warnings
import pandas
import pyranges
import sys
import requests
import numpy as np
import pybiomart as pbm
from scenicplus.scenicplus_class import create_SCENICPLUS_object
from scenicplus.wrappers.run_scenicplus import run_scenicplus
from scenicplus.preprocessing.filtering import apply_std_filtering_to_eRegulons
from scenicplus.eregulon_enrichment import score_eRegulons

ensembl_version_dict = {'109': 'http://www.ensembl.org',
                                '108': 'http://oct2022.archive.ensembl.org/',
                                '107': 'http://jul2022.archive.ensembl.org/',
                                '106': 'http://apr2022.archive.ensembl.org/',
                                '105': 'http://dec2021.archive.ensembl.org/',
                                '104': 'http://may2021.archive.ensembl.org/',
                                '103': 'http://feb2021.archive.ensembl.org/',
                                '102': 'http://nov2020.archive.ensembl.org/',
                                '101': 'http://aug2020.archive.ensembl.org/',
                                '100': 'http://apr2020.archive.ensembl.org/',
                                '99': 'http://jan2020.archive.ensembl.org/',
                                '98': 'http://sep2019.archive.ensembl.org/',
                                '97': 'http://jul2019.archive.ensembl.org/',
                                '96': 'http://apr2019.archive.ensembl.org/',
                                '95': 'http://jan2019.archive.ensembl.org/',
                                '94': 'http://oct2018.archive.ensembl.org/',
                                '93': 'http://jul2018.archive.ensembl.org/',
                                '92': 'http://apr2018.archive.ensembl.org/',
                                '91': 'http://dec2017.archive.ensembl.org/',
                                '90': 'http://aug2017.archive.ensembl.org/',
                                '89': 'http://may2017.archive.ensembl.org/',
                                '88': 'http://mar2017.archive.ensembl.org/',
                                '87': 'http://dec2016.archive.ensembl.org/',
                                '86': 'http://oct2016.archive.ensembl.org/',
                                '80': 'http://may2015.archive.ensembl.org/',
                                '77': 'http://oct2014.archive.ensembl.org/',
                                '75': 'http://feb2014.archive.ensembl.org/',
                                '54': 'http://may2009.archive.ensembl.org/'}


def test_ensembl_host(scplus_obj, host, species):
    dataset = pbm.Dataset(name=species+'_gene_ensembl',  host=host)
    annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name', 'transcript_biotype'])
    annot.columns = ['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    annot['Chromosome'] = annot['Chromosome'].astype('str')
    filter = annot['Chromosome'].str.contains('CHR|GL|JH|MT')
    annot = annot[~filter]
    annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    gene_names_release = set(annot['Gene'].tolist())
    ov=len([x for x in scplus_obj.gene_names if x in gene_names_release])
    return ov

def run_scenicplus_analysis(work_dir, adata, cistopic_obj, n_cpu, group_variable,
save_results, res_file, save_path):
    warnings.filterwarnings("ignore")
    # Set stderr to null to avoid strange messages from ray
    _stderr = sys.stderr
    null = open(os.devnull,'wb')
    tmp_dir = '/scratch/'
    menr = dill.load(open(os.path.join(work_dir, 'motifs/menr.pkl'), 'rb'))
    scplus_obj = create_SCENICPLUS_object(
        GEX_anndata = adata,
        cisTopic_obj = cistopic_obj,
        menr = menr
    )
    scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())
    
    n_overlap = {}
    for version in ensembl_version_dict.keys():
        try:
            n_overlap[version] =  test_ensembl_host(scplus_obj, ensembl_version_dict[version], 'hsapiens')
        except:
            pass
    v = sorted(n_overlap.items(), key=lambda item: item[1], reverse=True)[0][0]
    biomart_host = ensembl_version_dict[v]
    if not os.path.exists(os.path.join(work_dir, 'data')):
        os.makedirs(os.path.join(work_dir, 'data'))
     
    if not os.path.exists(os.path.join(work_dir, 'data', 'utoronto_human_tfs_v_1.01.txt')): 
        response = requests.get('http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt')
        output_file = open(os.path.join(work_dir, 'data', 'utoronto_human_tfs_v_1.01.txt'), 'wb')
        output_file.write(response.content)
        output_file.close()

    try:
        run_scenicplus(
            scplus_obj = scplus_obj,
            variable = ['GEX_' + group_variable],
            species = 'hsapiens',
            assembly = 'hg38',
            tf_file = os.path.join(work_dir, 'data', 'utoronto_human_tfs_v_1.01.txt'),
            save_path = os.path.join(work_dir, 'scenicplus'),
            biomart_host = biomart_host,
            upstream = [1000, 250000],
            downstream = [1000, 250000],
            calculate_TF_eGRN_correlation = True,
            calculate_DEGs_DARs = True,    
            export_to_loom_file = True,   
            export_to_UCSC_file = False,   #non-default, otherwise there is an error
            n_cpu = n_cpu,
            _temp_dir = os.path.join(tmp_dir, 'ray_spill'))
    except Exception as e:
        #in case of failure, still save the object
        dill.dump(scplus_obj, open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'wb'), protocol=-1)
        raise(e)
  
  # Kepp only eRegulons with an extended annotation if there is no direct annotation available 
  # (given that the confidence of direct motif annotations is in genral higher). 
  # Discard eRegulons for which the region-to-gene correlation is negative (these are often noisy). 
  # Rename the eRegulons so that eRegulons with the suffix TF_+_+ become TF_+ and those with TF_-_+ become TF_-.
    scplus_obj = dill.load(open(os.path.join(work_dir, 'scenicplus', 'scplus_obj.pkl'), 'rb'))
    apply_std_filtering_to_eRegulons(scplus_obj)

    region_ranking = dill.load(open(os.path.join(work_dir, 'scenicplus/region_ranking.pkl'), 'rb')) #load ranking calculated using the wrapper function
    gene_ranking = dill.load(open(os.path.join(work_dir, 'scenicplus/gene_ranking.pkl'), 'rb')) #load ranking calculated using the wrapper function
    score_eRegulons(scplus_obj,
                    ranking = region_ranking,
                    eRegulon_signatures_key = 'eRegulon_signatures_filtered',
                    key_added = 'eRegulon_AUC_filtered',
                    enrichment_type= 'region',
                    auc_threshold = 0.05,
                    normalize = False,
                    n_cpu = n_cpu)
    score_eRegulons(scplus_obj,
                    gene_ranking,
                    eRegulon_signatures_key = 'eRegulon_signatures_filtered',
                    key_added = 'eRegulon_AUC_filtered',
                    enrichment_type = 'gene',
                    auc_threshold = 0.05,
                    normalize= False,
                    n_cpu = n_cpu)
    scplus_obj.uns['eRegulon_AUC_filtered']['Region_based'].columns = scplus_obj.uns['eRegulon_AUC_filtered']['Region_based'].columns.str.replace("+","act")
    scplus_obj.uns['eRegulon_AUC_filtered']['Region_based'].columns = scplus_obj.uns['eRegulon_AUC_filtered']['Region_based'].columns.str.replace("-","supp")
    scplus_obj.uns['eRegulon_AUC_filtered']['Gene_based'].columns = scplus_obj.uns['eRegulon_AUC_filtered']['Gene_based'].columns.str.replace("+","act")
    scplus_obj.uns['eRegulon_AUC_filtered']['Gene_based'].columns = scplus_obj.uns['eRegulon_AUC_filtered']['Gene_based'].columns.str.replace("-","supp")
    if save_results:
        scplus_obj.uns['eRegulon_AUC_filtered']['Gene_based'].to_csv(os.path.join(save_path, res_file+"_genes.csv"))
        scplus_obj.uns['eRegulon_AUC_filtered']['Region_based'].to_csv(os.path.join(save_path, res_file+"_region.csv"))
    return list(scplus_obj.uns['eRegulon_AUC_filtered']['Gene_based'], 
    scplus_obj.uns['eRegulon_AUC_filtered']['Region_based'])


