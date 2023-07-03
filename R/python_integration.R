#' @import basilisk
#' @export
get_base_GRN <- function(all_peaks, conns, association_cutoff = 0.8){
    proc = basiliskStart(venv)
    on.exit(basiliskStop(proc))
    base_GRN <- basiliskRun(proc, function(peaks, connections, cutoff){
        source_python(system.file("python/base_GRN.py", package = "epiregulon.benchmark"))
        output_GRN <- prepare_base_GRN(peaks, connections, cutoff)
        output_GRN
    }, peaks = all_peaks, connections = conns, cutoff = association_cutoff)
    base_GRN
}

#' @export
get_annotated_data <- function(barcodes_list, HTO_list, paths_to_data_files, sample_names=NULL, n_top_genes=2000){
    proc = basiliskStart(venv)
    on.exit(basiliskStop(proc))
    adata <- basiliskRun(proc, function(barcodes, HTO){
        source_python(system.file("python/gene_expression_data.py", package = "epiregulon.benchmark"))
        adata <- get_annotated_data(barcodes, HTO, paths, samples, n_top_genes)
        adata
    }, barcodes = barcodes_list,
    HTO = HTO_list,
    paths = paths_to_data_files,
    samples = sample_names,
    n_top_genes = n_top_genes)
    adata
}

#' @export
calculate_GRN <- function(adata, base_GRN){
    proc = basiliskStart(venv)
    on.exit(basiliskStop(proc))
    cellorcacle_GRN <- basiliskRun(proc, function(geneExpr, GRN){
        source_python(system.file("python/calculate_GRN.py", package = "epiregulon.benchmark"))
        complete_GRN <- calculate_GRN(geneExpr, GRN)
        complete_GRN
    }, geneExpr = adata, GRN = base_GRN)
    cellorcacle_GRN
}
