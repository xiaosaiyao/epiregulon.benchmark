#' Generate the report with the use of reprogram-seq data
#'
#' @param include_other_packages A logical indicating whether the raport should present
#' benchmarking the results using other tools i.e. FigR, Pando, cellOracle and scenic+.
#' @param regulon A data frame being the output from epiregulon workflow. The weight
#' column, if present, will be disregarded.
#' @param regulon_script A full path to the R script file with the code generating the regulon. It
#' will not be evaluated.
#' @param treatments A named list of character vectors. Each name corresponds to
#' one treatment and each vector to respective control groups. Treatments and control groups
#' are matched to the values in the `hash_assignment` column from the `GeneExpressionMatrix`.
#' If the name points to the `NULL` object then control set of cell is created from all the groups
#' apart from that indicated by the element name. List names should contain transcription factors' names.
#' The default is set to `list(GATA6 = NULL, NKX2.1 = NULL)`.
#' @param cluster_column A character specifying the column in the `GeneExpressionMatrix` in which cluster information
#' is stored. By default is is `Clusters`.
#' @param output_file_summary A path to the .html file to which report will be saved.
#' @param output_file_summary A path to the .csv file to which data frame with the summary of the analysis output will be saved.
#' @export
render_report_reprogram_seq <- function(...){
    rmarkdown::render(system.file("reports/Reprogram-seq.Rmd", package = "epiregulon.benchmark"), output_format = "all", output_file = output_file,
                      params = list(...))
}

#' @export
render_report_NGS4557 <- function(data_file_paths = "",
                                  temp_folder = "",
                                  sample_names = "",
                                  conda_exe = "",
                                  virtual_env = "",
                                  work_dir = ""){
    rmarkdown::render("NGS4557.Rmd", output_format = "all", output_file = "NGS4557.html",
                      params = list(data_file_paths = data_file_paths,
                                    temp_folder = temp_folder,
                                    sample_names = sample_names,
                                    conda_exe = conda_exe,
                                    virtual_env = virtual_env,
                                    work_dir = work_dir))
}
