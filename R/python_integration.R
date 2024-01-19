#' @import basilisk
#' @export
get_base_GRN <- function(all_peaks, conns, association_cutoff = 0.8){
    proc = basiliskStart(venv1)
    on.exit(basiliskStop(proc))
    base_GRN <- basiliskRun(proc, function(peaks, connections, cutoff){
        reticulate::source_python(system.file("python/base_GRN.py", package = "epiregulon.benchmark"))
        output_GRN <- prepare_base_GRN(peaks, connections, cutoff)
        output_GRN
    }, peaks = all_peaks, connections = conns, cutoff = association_cutoff)
    base_GRN
}

#' @export
get_annotated_data <- function(barcode_tab, paths_to_data_files, sample_names=NULL,
                               n_top_genes=2000, variable){
    proc = basiliskStart(venv1)
    on.exit(basiliskStop(proc))
    adata <- basiliskRun(proc, function(barcode_tab, paths, samples, n_top_genes, variable){
        reticulate::source_python(system.file("python/gene_expression_data.py", package = "epiregulon.benchmark"))
        adata <- get_annotated_data(barcode_tab, as.list(paths), as.list(samples), as.integer(n_top_genes),
                                    variable)
        adata
    }, barcode_tab = barcode_tab,
    paths = paths_to_data_files,
    samples = sample_names,
    n_top_genes = n_top_genes,
    variable = variable)
    adata
}

#' @export
calculate_GRN <- function(adata, base_GRN){
    proc = basiliskStart(venv1)
    on.exit(basiliskStop(proc))
    cellorcacle_GRN <- basiliskRun(proc, function(geneExpr, GRN){
        reticulate::source_python(system.file("python/calculate_GRN.py", package = "epiregulon.benchmark"))
        complete_GRN <- calculate_GRN(geneExpr, GRN)
        complete_GRN
    }, geneExpr = adata, GRN = base_GRN)
    cellorcacle_GRN
}

find_topics <- function(proc, barcode_tab, sample_names,
                        paths_to_fragments, work_dir,
                        tmp_dir, paths_to_peak_matrix, n_cpu,
                        group_variable,
                        save_results, file_name, save_path,
                        motif_db_dir, dataset, n_top_genes){
    obj_list <- basiliskRun(proc, function(barcode_tab, sample_names,
                               paths_to_fragments, work_dir,
                               tmp_dir, paths_to_peak_matrix, n_cpu,
                               group_variable,
                               save_results, file_name, save_path,
                               motif_db_dir, dataset, n_top_genes){
        reticulate::source_python(system.file("python/topics.py", package = "epiregulon.benchmark"))
        find_topics_python(barcode_tab, sample_names,
                    paths_to_fragments, work_dir,
                    tmp_dir, paths_to_peak_matrix, n_cpu,
                    group_variable,
                    save_results, file_name, save_path,
                    motif_db_dir, dataset, n_top_genes)
    }, barcode_tab=barcode_tab, sample_names = sample_names, paths_to_fragments = paths_to_fragments,
    work_dir = work_dir, tmp_dir = tmp_dir,
    paths_to_peak_matrix = paths_to_peak_matrix, n_cpu = n_cpu, group_variable = group_variable,
    save_results = save_results, file_name = file_name, save_path = save_path, motif_db_dir=motif_db_dir,
    dataset = dataset, n_top_genes=n_top_genes)
    obj_list
}

#' @export
run_scenicplus <- function(barcode_tab, sample_names, paths_to_fragments, work_dir,
                           tmp_dir, paths_to_peak_matrix, n_cpu, group_variable,
                           save_results, file_name, save_path, motif_db_dir,
                           dataset, n_top_genes, scenicplus_res_file){
    proc = basiliskStart(venv2)
    on.exit(basiliskStop(proc))
    test <- try(
        {obj_list <- find_topics(proc=proc, barcode_tab=barcode_tab, sample_names=sample_names,
                                 paths_to_fragments=paths_to_fragments, work_dir=work_dir,
                                 tmp_dir=tmp_dir,
                                 paths_to_peak_matrix=paths_to_peak_matrix, n_cpu=n_cpu,
                                 group_variable=group_variable,
                                 save_results=save_results, file_name=file_name, save_path=save_path,
                                 motif_db_dir=motif_db_dir,
                                 dataset=dataset, n_top_genes=n_top_genes)})
    if(is(test, "try-error")){
        msg <- attr(test, "condition")$message
        if(grepl("ImportError:.*version.*LIB.*not found", msg)) {
            setBasiliskForceFallback(TRUE)
            obj_list<-find_topics(barcode_tab, sample_names, paths_to_fragments, work_dir, tmp_dir,
                        paths_to_peak_matrix, n_cpu, group_variable,
                        save_results, file_name, save_path, motif_db_dir,
                        dataset, n_top_genes)
        }
    }
    scenicplus_res <- basiliskRun(proc, function(work_dir,
                                                 cistopic_obj, n_cpu,
                                                 group_variable, save_results,
                                                 res_file, save_path){
        reticulate::source_python(system.file("python/scenicplus_main.py", package = "epiregulon.benchmark"))
        run_scenicplus_analysis(work_dir, adata, cistopic_obj, n_cpu, group_variable,
                                save_results, res_file, save_path)
    }, work_dir = work_dir, adata = obj_list[[2]], cistopic_obj = obj_list[[1]],
    n_cpu = n_cpu, group_variable = group_variable, save_results = save_results,
    res_file = res_file, save_path = save_path)
    scenicplus_res
}
