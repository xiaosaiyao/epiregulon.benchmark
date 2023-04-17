#' @export
calculate_AUC <- function(x,y){
 y <- y[order(x)]
 x <- x[order(x)]
 x_intervals <- diff(x)
 pair_mean_y <- (y[1:(length(y)-1)] + y[2:length(y)])/2
 sum(x_intervals*pair_mean_y)
}

#' @export
calculate_accuracy_metrics <- function(values, positive_elements_ind, negative_elements_ind, n_steps = 1e3){
    values <- as.vector(values)
    values <- values[unique(c(positive_elements_ind, negative_elements_ind))]
    positive_elements_ind <- which(positive_elements_ind %in% unique(c(positive_elements_ind, negative_elements_ind)))
    max_val <- max(values)
    min_val <- min(values)
    steps <- seq(min_val, max_val, length.out = 1e3)
    is_positive <- rep(FALSE, length(values))
    is_positive[positive_elements_ind] <- TRUE
    res_list <- list()
    all_combinations <- data.frame(threshold_reached=as.logical(c(1,1,0,0)), category = as.logical(c(1,0,1,0)))
    for(i in seq_along(steps)){
        threshold_reached <- values >= steps[i]
        confusion_matrix <- data.frame(threshold_reached =threshold_reached, category = is_positive)
        confusion_matrix <- rbind(confusion_matrix, all_combinations)
        tab <- table(confusion_matrix)
        # account for adding one observation for combination
        tab <- tab - 1
        res_list[[i]] <- c(TP = tab[2,2], FP = tab[2, 1], TN = tab[1, 1],
                           FN = tab[1, 2], cutoff = steps[i])
    }
    res <- as.data.frame(do.call(rbind, res_list))
    TPR <- rev(res$TP/(res$TP + res$FN))
    FPR <- rev(res$FP/(res$FP + res$TN))
    list(TPR = TPR, FPR = FPR, cutoff = res$cutoff, confusion_matrix_data = res)
}

#' @export
prepare_plot_data <- function(regulon, weight.args, group_combinations, geneExprMatrix.sce,
                              motif_score = TRUE, ArchR_proj_path = NULL,
                              weight_clusters,
                              activity_clusters,
                              use_cell_line_weights = TRUE, group_column = "Cellline",
                              path_to_archR_project = NULL){
    weight.args$regulon <- regulon
    weight.args$expMatrix <- geneExprMatrix.sce
    weight.args$clusters <- geneExprMatrix.sce[[weight_clusters]]
    regulon.w <- do.call(addWeights, weight.args)
    if (motif_score){
        if (is.null(path_to_archR_project)){
            regulon.w <- addMotifScore(regulon.w,
                                       peaks = rowRanges(weight.args$peakMatrix),
                                       species="human",
                                       genome="hg38")

        }
        else{
            regulon.w <- addMotifScore(regulon.w,
                                       species="human",
                                       genome="hg38",
                                       archr_path = path_to_archR_project)

        }


        # adjust action to whether weight is a matrix or vector
        if(is.null(dim(regulon.w$weight)))
            regulon.w$weight[regulon.w$motif == 0] <- 0
        else
            regulon.w$weight[regulon.w$motif==0,] <- 0
    }
    activity_cluster_labels = NULL
    if(!is.na(activity_clusters))
        activity_cluster_labels <- geneExprMatrix.sce[[activity_clusters]]
    activity.matrix <- calculateActivity(regulon = regulon.w,
                                         expMatrix = geneExprMatrix.sce,
                                         clusters  = activity_cluster_labels)
    activity.matrix <- activity.matrix[,colnames(geneExprMatrix.sce)]
    plot_data <- data.frame()
    for(i in 1:nrow(group_combinations)){
        selected_row <- group_combinations[i,]
        current_combination <- as.list(selected_row)
        cell_line_ind <- grep(current_combination$cell_line, geneExprMatrix.sce$hash_assignment)
        activity_values <- activity.matrix[current_combination$transcription_factor,]
        experimental_treatment_ind <- grep(current_combination$experimental_treatment, geneExprMatrix.sce$hash_assignment)
        control_treatment_ind <- grep(current_combination$control_treatment, geneExprMatrix.sce$hash_assignment)
        if (current_combination$effect_direction == "up"){
            positive_group <- intersect(cell_line_ind, experimental_treatment_ind)
            negative_group <- intersect(cell_line_ind, control_treatment_ind)
        }
        else{
            positive_group <- intersect(cell_line_ind, control_treatment_ind)
            negative_group <- intersect(cell_line_ind, experimental_treatment_ind)
        }
        accuracy_stats <- calculate_accuracy_metrics(activity_values, positive_group, negative_group)
        partial_data <- data.frame(TPR = accuracy_stats$TPR,
                                   FPR = accuracy_stats$FPR,
                                   cell_line = current_combination$cell_line,
                                   treatment = current_combination$experimental_treatment,
                                   tf = current_combination$transcription_factor)
        partial_data[["method"]] <- weight.args$method

        plot_data <- rbind(plot_data, partial_data)
    }
    plot_data
}

#' @export

get_activity_matrix <- function(geneExprMatrix.sce,
                                 method = "FigR",
                                 GRN = NULL,
                                 Seurat_obj = NULL,
                                 n_bin =24,
                                 tfs = NULL){
    stopifnot(method %in% c("FigR", "Epiregulon", "cellOracle", "Pando", "Scenic"))
    library(epiregulon)
    if(method == "FigR"){
        # adjust gene names to FigR
        GRN <- GRN[,c("Motif", "DORC", "Score")]
        colnames(GRN) <- c("tf", "target", "weight")
        return(epiregulon::calculateActivity(expMatrix = geneExprMatrix.sce,
                                                         regulon = GRN))
    }
    else if(method == "Epiregulon"){
        return(epiregulon::calculateActivity(expMatrix = geneExprMatrix.sce,
                                                         regulon = GRN))
    }
    else if(method == "Pando"){
        library(Signac)
        library(Seurat)
        library(Pando)
        test_srt <- find_modules(Seurat_obj, rsq_thresh = 0.05)
        TFmodules <- NetworkModules(test_srt)
        DefaultAssay(Seurat_obj) <- "RNA"
        tfs <- tfs[tfs %in% names(TFmodules@features$genes_pos)]
        activity.matrix <- matrix(nrow = length(tfs), ncol = dim(Seurat_obj)[2])
        for(i in seq_along(tfs)){
            tf <- tfs[i]
            # extract target genes (only upregulated)
            targets <- TFmodules@features$genes_pos[[tf]]
            # adjust gene name to Pando requirements
            tf <- tfs[i] <- names(positive_clusters)[i] <- gsub("-", "", tf)
            Seurat_obj <- AddModuleScore(Seurat_obj, features = list(targets), name = paste0(tf, "_activity"),
                                         nbin = n_bin)
            activity.matrix[i,] <- Seurat_obj@meta.data[[paste0(tf, "_activity1")]]
        }
        rownames(activity.matrix) <- tfs
        return(activity.matrix)
    }
    else if(method == "cellOracle"){
        regulon <- GRN[,c("target","source", "coef_mean", "cluster")]
        n_tf <- length(unique(regulon$source))
        colnames(regulon) <- c("target", "tf", "weight", "cluster")
        tfs <- tfs[tfs %in% regulon$tf]
        regulon <- split(regulon, regulon$cluster)
        colnames(clusters) <- c("barcode","cluster_id")
        unique_clusters <- as.character(unique(clusters$cluster_id))
        clusters <- split(clusters, clusters$cluster_id)
        activity.matrix <- matrix(ncol = ncol(assay(geneExprMatrix.sce)), nrow  = n_tf)
        cell_ind <- 1
        matrix_columns <- c()
        for (cluster_id in unique_clusters){
            cluster_cells <- clusters[[cluster_id]]$barcode
            activity_part <- epiregulon::calculateActivity(expMatrix = geneExprMatrix.sce[,cluster_cells],
                                                           regulon = regulon[[cluster_id]])
            activity.matrix[,cell_ind:(cell_ind+length(cluster_cells)-1)] <- as.matrix(activity_part)
            if(is.null(rownames(activity.matrix))) rownames(activity.matrix) <- rownames(activity_part)
            cell_ind <- cell_ind + length(cluster_cells)
            matrix_columns <- c(matrix_columns, cluster_cells)
        }
        # adjust row order to sce object
        colnames(activity.matrix) <- matrix_columns
        return(activity.matrix)
    }
}
