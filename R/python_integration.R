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
get_annotated_data <- function(barcode_tab, paths_to_data_files, sample_names=NULL,
                               n_top_genes=2000, python_library = "Scenic_plus",
                               variable){
    proc = basiliskStart(venv)
    on.exit(basiliskStop(proc))
    adata <- basiliskRun(proc, function(barcode_tab, paths, samples, n_top_genes, python_lib){
        source_python(system.file("python/gene_expression_data.py", package = "epiregulon.benchmark"))
        adata <- get_annotated_data(barcode_tab, as.list(paths), as.list(samples), as.integer(n_top_genes), python_lib,
                                    variable)
        adata
    }, barcode_tab = barcode_tab,
    paths = paths_to_data_files,
    samples = sample_names,
    n_top_genes = n_top_genes,
    python_lib = python_library,
    variale = variable)
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
