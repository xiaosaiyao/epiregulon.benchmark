#' @export
get_base_GRN <- function(all_peaks, conns){
    proc = basiliskStart(venv)
    on.exit(basiliskStop(proc))
    base_GRN <- basiliskRun(proc, function(peaks, connections){
        source_python(system.file("python/base_GRN.py", package = "epiregulon.benchmark"))
        output_GRN <- prepare_base_GRN(peaks, connections)
        output_GRN
    }, peaks = all_peaks, connections = conns)
    base_GRN
}

#' @export
get_annotated_data <- function(barcodes_list, HTO_list){
    proc = basiliskStart(venv)
    on.exit(basiliskStop(proc))
    adata <- basiliskRun(proc, function(barcodes, HTO){
        source_python(system.file("python/gene_expression_data.py", package = "epiregulon.benchmark"))
        adata <- get_annotated_data(barcodes, HTO)
        adata
    }, barcodes = barcodes_list, HTO = HTO_list)
    adata
}

#' @export
calculate_GRN <- function(adata, HTO_list){
    proc = basiliskStart(venv)
    on.exit(basiliskStop(proc))
    cellorcacle_GRN <- basiliskRun(proc, function(geneExpr, HTO){
        source_python(system.file("python/calculate_GRN.py", package = "epiregulon.benchmark"))
        GRN <- calculate_GRN(geneExpr, HTO)
        GRN
    }, geneExpr = adata, HTO = HTO_list)
    lapply(cellorcacle_GRN, py_to_r)
}
