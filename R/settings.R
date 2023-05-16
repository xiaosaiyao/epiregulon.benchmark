store_intermediates <- function(path_to_intermediates){
    res = NULL
    if(!dir.exists(path_to_intermediates))
        while(! res %in% c("1","2")){
            res = readline(prompt=sprintf("Directory %s does not exist. Create a new directory?\n1 : no\n2: yes", path_to_intermediates))
        }
    if(res  == "2") dir.create(path_to_intermediates, recursive = TRUE)
    else return()
    # protect from overwriting existing files
    if(length(list.files(path_to_intermediates, all.files = TRUE))!=0 & is.null(getOption("epiregulon.benchmark_new_intermediates")))
        stop("Folder is non empty. Consider setting 'options(epiregulon.benchmark_new_intermediates = TRUE)' to allow for the overwriting existing files.")
    options(epiregulon.benchmark_path_to_intermediates = path_to_intermediates)
}

.check_intermediate <- function(file_path){
    overwrite = getOption("epiregulon.benchmark_new_intermediates")
    if(overwrite) return(FALSE)
    intermediates_folder = getOption("epiregulon.benchmark_path_to_intermediates")
    if(is.null(intermediates_folder)) return(FALSE)
    else return(file.exists(file.path(intermediates_folder, file_path)))
}
