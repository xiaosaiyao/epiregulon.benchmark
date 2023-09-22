#' Generate the report with the use of reprogram-seq data
#'
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

#' @param motif_score A logical indicating whether motif scores should be added. If so, regulon rows with zero motif score
#' will have also 0 value for the weight.
#' @export
render_report_reprogram_seq <- function(output_file = file.path(getwd(),"Reprogram-seq_benchmark.html"), ...){
    rmarkdown::render(system.file("reports/Reprogram-seq.Rmd", package = "epiregulon.benchmark"), output_format = "all", output_file = output_file,
                      params = list(...))
}


#' @param path_to_ArchR_proj Absolute path the to folder with saved ArchR project in which AR dataset was analysed
#' @export
render_report_AR <- function(output_file = file.path(getwd(), "AR.html"), path_to_ArchR_proj, ...){
    rmarkdown::render(system.file("reports/Reprogram-seq.Rmd", package = "epiregulon.benchmark"), output_format = "all", output_file = output_file,
                      params = c(list(path_to_ArchR_proj = path_to_ArchR_proj), list(...)))
}
