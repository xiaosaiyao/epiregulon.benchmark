#' Generate the report with the use of reprogram-seq data
#'
#' @param treatments A named list of character vectors. Each name corresponds to
#' one treatment and each vector to respective control groups. Treatments and control groups
#' are matched to the values in the `hash_assignment` column from the `GeneExpressionMatrix`.
#' If the name points to the `NULL` object then control set of cell is created from all the groups
#' apart from that indicated by the element name. List names should contain transcription factors' names.
#' The default is set to `list(GATA6 = NULL, NKX2.1 = NULL)`.
#' @param motif_score A logical indicating whether motif scores should be added. If so, regulon rows with zero motif score
#' will have also 0 value for the weight.
#' @param motif_db_dir Path to the direcotry with motif databases used by `pycistarget`. Can be downloaded from
#' https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/.
#' @export
render_report_reprogram_seq <- function(output_file="Reprogram-seq_benchmark.html",
                                        data_file_paths="/gne/data/lab-shares/xie-lab/Sequencing_Data/2022/mapping/20220124_ReprogramSeq_Multiome/JT65_67/outs/",
                                        motif_db_dir = find.package("epiregulon.benchmark"),...){
    rmarkdown::render(system.file("reports/Reprogram-seq.Rmd", package = "epiregulon.benchmark"), output_format="all", output_file=output_file,
                      params=c(list(data_file_paths = data_file_paths, motif_db_dir=motif_db_dir), list(...)))
}

#' @param path_to_ArchR_proj Absolute path the to folder with saved ArchR project in which AR dataset was analysed
#' @export
render_report_AR <- function(output_file="AR.html", path_to_ArchR_proj,
                             data_file_paths = c("/gstore/data/genomics/congee_rest_runs/62e2c81cd6d7e0bd49c579dd/SAM24418230/croo_output/",
                                                 "/gstore/data/genomics/congee_rest_runs/62e2c81cd6d7e0bd49c579dd/SAM24418231/croo_output/"),
                             motif_db_dir = find.package("epiregulon.benchmark"),...){
    rmarkdown::render(system.file("reports/AR.Rmd", package="epiregulon.benchmark"), output_format="all", output_file=output_file,
                      params=c(list(data_file_paths = data_file_paths, path_to_ArchR_proj=path_to_ArchR_proj,
                                    motif_db_dir=motif_db_dir), list(...)))
}

#' @rdname render_report_reprogram_seq
#' @export
render_report_simulation <- function(output_file = "AR.html", ...){
    rmarkdown::render(system.file("reports/Simulation_benchmark.Rmd", package="epiregulon.benchmark"), output_format="all", output_file=output_file,
                      params=list(...))
}
