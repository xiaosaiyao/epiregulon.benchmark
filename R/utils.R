#' @export
calculate_AUC <- function(x,y){
    y <- y[order(x)]
    x <- x[order(x)]
    non_unique_x_ind <- which(duplicated(x))
    non_unique_x_ind <- sort(unique(c(non_unique_x_ind, non_unique_x_ind-1)))
    y[non_unique_x_ind] <- sort(y[non_unique_x_ind])
    x_intervals <- diff(x)
    pair_mean_y <- (y[1:(length(y)-1)] + y[2:length(y)])/2
    sum(x_intervals*pair_mean_y)
}

#' @export
calculate_accuracy_metrics <- function(values, positive_elements_ind, negative_elements_ind, n_steps = 1e3){
    values <- as.vector(values)
    values <- values[unique(c(positive_elements_ind, negative_elements_ind))]
    positive_elements_ind <- which(unique(c(positive_elements_ind, negative_elements_ind)) %in% positive_elements_ind)
    max_val <- max(values)
    min_val <- min(values)
    max_val <- max_val+(max_val-min_val)/n_steps #add to account for ">=" used for threshold
    steps <- seq(min_val, max_val, length.out = n_steps)
    is_positive <- rep(FALSE, length(values))
    is_positive[positive_elements_ind] <- TRUE
    res_list <- list()
    all_combinations <- data.frame(threshold_reached=as.logical(c(1,1,0,0)), category = as.logical(c(1,0,1,0)))
    for(i in seq_along(steps)){
        threshold_reached <- values >= steps[i]
        confusion_matrix <- data.frame(threshold_reached =threshold_reached, category = is_positive)
        confusion_matrix <- rbind(all_combinations, confusion_matrix)
        tab <- table(confusion_matrix)
        # account for adding one observation for combination
        tab <- tab - 1
        res_list[[i]] <- c(TP = tab[2,2], FP = tab[2, 1], TN = tab[1, 1],
                           FN = tab[1, 2], cutoff = steps[i])
    }
    res <- as.data.frame(do.call(rbind, res_list))
    TPR <- res$TP/(res$TP + res$FN)
    FPR <- res$FP/(res$FP + res$TN)
    list(TPR = TPR, FPR = FPR, cutoff = res$cutoff, confusion_matrix_data = res)
}

#' @export
prepare_plot_data <- function(regulon, weight.args, group_combinations, geneExprMatrix.sce,
                              motif_score = TRUE, archr_path = NULL,
                              weight_clusters,
                              activity_clusters=NULL,
                              use_cell_line_weights = TRUE, group_column = "Cellline",
                              treatment_column = "hash_assignment",
                              ...){
    weight.args$regulon <- regulon
    weight.args$expMatrix <- geneExprMatrix.sce
    if (!is.null(weight_clusters))
        weight.args$clusters <- as.vector(geneExprMatrix.sce[[weight_clusters]])
    if(weight.args$method == "logFC"){
        assays(weight.args$expMatrix, withDimnames = TRUE)[[weight.args$exp_assay]] <- log2(assays(weight.args$expMatrix, withDimnames = TRUE)[[weight.args$exp_assay]]+1)
    }
    regulon.w <- do.call(addWeights, weight.args)
    if (motif_score){
        if (is.null(archr_path)){
            regulon.w <- addMotifScore(regulon.w,
                                       peaks = rowRanges(weight.args$peakMatrix),
                                       species="human",
                                       genome="hg38")

        }
        else{
            regulon.w <- addMotifScore(regulon.w,
                                       species="human",
                                       genome="hg38",
                                       archr_path = archr_path)
        }
        # adjust action to whether weight is a matrix or vector
        if(is.null(dim(regulon.w$weight)))
            regulon.w$weight[regulon.w$motif == 0 | is.na(regulon.w$motif)] <- 0
        else
            regulon.w$weight[regulon.w$motif==0 | is.na(regulon.w$motif),] <- 0
    }
    activity_cluster_labels = NULL
    if(!is.na(activity_clusters))
        activity_cluster_labels <- as.vector(geneExprMatrix.sce[[activity_clusters]])

    if(!is.null(weight.args$exp_assay)) exp_assay <- weight.args$exp_assay
    else exp_assay <- formals(calculateActivity)$exp_assay
    activity.matrix <- calculateActivity(regulon = regulon.w,
                                         expMatrix = geneExprMatrix.sce,
                                         clusters  = activity_cluster_labels,
                                         exp_assay = exp_assay,...)

    activity.matrix <- activity.matrix[,colnames(geneExprMatrix.sce),drop=FALSE]
    plot_data <- data.frame()
    for(i in 1:nrow(group_combinations)){
        selected_row <- group_combinations[i,]
        current_combination <- as.list(selected_row)
        if (!as.character(current_combination$transcription_factor) %in% rownames(activity.matrix)) next
        if(!is.na(current_combination$cell_line))
            cell_line_ind <- grep(current_combination$cell_line, geneExprMatrix.sce$cell_line)
        else cell_line_ind <- seq_len(dim(geneExprMatrix.sce)[2])
        activity_values <- activity.matrix[as.character(current_combination$transcription_factor),]
        experimental_treatment_ind <- grep(current_combination$experimental_treatment, geneExprMatrix.sce[[treatment_column]])
        control_treatment_ind <- grep(current_combination$control_treatment, geneExprMatrix.sce[[treatment_column]])
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

#' @import epiregulon.archr
#' @export
get_activity_matrix <- function(method = NULL,
                                GRN = NULL,
                                tfs = NULL,
                                geneExprMatrix.sce = NULL,
                                exp_assay = "logcounts"){
    if(!method %in% c("FigR", "Epiregulon", "cellOracle", "Pando","GRaNIE")) return(NULL)
    library(epiregulon)
    if(method == "FigR"){
        # adjust gene names to FigR
        GRN <- GRN[,c("Motif", "DORC", "Score")]
        colnames(GRN) <- c("tf", "target", "weight")
        if(length(intersect(tfs, GRN$tf))==0) {
            warning("Tfs not found in the FigR output.")
            return(NULL)
        }
        return(calculateActivity(expMatrix = geneExprMatrix.sce,
                                 regulon = GRN,
                                 exp_assay = exp_assay))
    }
    else if(method == "GRaNIE"){
        GRN <- GRN[,c("TF.name", "gene.name")]
        GRN$weight <- 1 # set weights to 1
        colnames(GRN) <- c("tf", "target", "weight")
        if(length(intersect(tfs, GRN$tf))==0) {
            warning("Tfs not found in the GRaNIE output.")
            return(NULL)
        }
        return(calculateActivity(expMatrix = geneExprMatrix.sce,
                                 regulon = GRN,
                                 exp_assay = exp_assay))
    }
    else if(method == "Epiregulon"){
        return(calculateActivity(expMatrix = geneExprMatrix.sce,
                                 regulon = GRN,
                                 exp_assay = exp_assay))
    }
    else if(method == "Pando"){
        library(Signac)
        library(Seurat)
        library(Pando)
        Seurat_obj <- GRN
        test_srt <- find_modules(Seurat_obj, rsq_thresh = 0.05)
        TFmodules <- NetworkModules(test_srt)
        DefaultAssay(Seurat_obj) <- "RNA"
        tfs <- tfs[tfs %in% names(TFmodules@features$genes_pos)]
        if(length(tfs)==0) {
            warning("Tfs not found in the Pando output.")
            return(NULL)
        }
        activity.matrix <- matrix(nrow = length(tfs), ncol = dim(Seurat_obj)[2])
        for(i in seq_along(tfs)){
            tf <- tfs[i]
            # extract target genes (only upregulated)
            targets <- TFmodules@features$genes_pos[[tf]]
            # adjust gene name to Pando requirements
            tf <- tfs[i] <- gsub("-", "", tf)
            Seurat_obj_2 <- NULL
            n_bin=24
            while(is.null(Seurat_obj_2)){
                Seurat_obj_2 <- tryCatch({AddModuleScore(Seurat_obj, features = list(targets), name = paste0(tf, "_activity"),
                                                         nbin = n_bin)},
                                         error = function(cond){message(cond); message("Trying with ", n_bin-1, " bins")},
                                         finally = {n_bin = n_bin - 1})

            }

            activity.matrix[i,] <- Seurat_obj_2@meta.data[[paste0(tf, "_activity1")]]
        }
        rownames(activity.matrix) <- tfs
        colnames(activity.matrix) <- colnames(Seurat_obj_2)
        return(activity.matrix)
    }
    else if(method == "cellOracle"){
        clusters <- colData(geneExprMatrix.sce)$cluster_cellOracle
        regulon <- GRN[,c("target","source", "coef_mean", "cluster")]
        colnames(regulon) <- c("target", "tf", "weight", "cluster")
        tfs <- tfs[tfs %in% regulon$tf]
        if(length(tfs)==0) {
            warning("Tfs not found in the cellOracle output.")
            return(NULL)
        }
        regulon <- regulon[regulon$tf %in% tfs,]
        regulon <- split(regulon, regulon$cluster)
        unique_clusters <- as.character(unique(clusters))
        cluster_cells <- split(colnames(geneExprMatrix.sce), clusters)
        activity.matrix <- matrix(ncol = ncol(assay(geneExprMatrix.sce)), nrow  = length(tfs), dimnames = list(tfs))
        cell_ind <- 1
        matrix_columns <- c()
        for (cluster_id in unique_clusters){
            selected_cells <- cluster_cells[[cluster_id]]
            activity_part <- calculateActivity(expMatrix = geneExprMatrix.sce[,selected_cells],
                                               regulon = regulon[[cluster_id]],
                                               exp_assay = exp_assay)
            if(is.null(activity_part))
                activity_part <- matrix(0, ncol = length(selected_cells),
                                        nrow  = length(tfs),
                                        dimnames = list(tfs))
            # add values for tfs which have 0 activity in the cluster
            if (nrow(activity_part) < nrow(activity.matrix)){
                missing_tfs <- setdiff(tfs, rownames(activity_part))
                activity_part <- rbind(activity_part, matrix(0, ncol = ncol(activity_part), nrow = length(missing_tfs), dimnames = list(missing_tfs)))
            }
            activity_part <- activity_part[tfs,,drop = FALSE]
            activity.matrix[,cell_ind:(cell_ind+length(selected_cells)-1)] <- as.matrix(activity_part)
            cell_ind <- cell_ind + length(selected_cells)
            matrix_columns <- c(matrix_columns, selected_cells)
        }
        # adjust row order to sce object
        colnames(activity.matrix) <- matrix_columns
        activity.matrix <- activity.matrix[Matrix::rowSums(activity.matrix)!=0,,drop= FALSE]
        return(activity.matrix)
    }
}

#' @export
getResultsFromActivity <- function(activity.matrix, add_plot=FALSE, tf,
                                   labels,
                                   positive_elements_label,
                                   negative_elements_label, ...){
    if(length(positive_elements_label) > 1)
        positive_elements_label <- paste(positive_elements_label, collapse = "|", sep ="")
    if(length(negative_elements_label) > 1)
        negative_elements_label <- paste(negative_elements_label, collapse = "|", sep ="")
    positive_elements_ind <- grep(positive_elements_label, labels)
    negative_elements_ind <- grep(negative_elements_label, labels)
    activity_values <- activity.matrix[as.character(tf), ]
    res <- calculate_accuracy_metrics(activity_values, positive_elements_ind,
                                      negative_elements_ind)
    if(add_plot)
        lines(res$TPR~res$FPR, ...)
    else
        plot(res$TPR~res$FPR, type ="l", xlab = "False postive rate",
             ylab = "True positive rate", ...)
    return(calculate_AUC(res$FPR, res$TPR))
}

#' @export
map_treatment_to_tf <- function(treatment, regulon){
    all_tfs <- unique(regulon$tf)
    filter_ind <- unlist(lapply(all_tfs, function(x) grepl(x, gsub("\\.", "-", treatment))))
    if(length(all_tfs[filter_ind])!=1) stop(sprintf("Treatment group %s should be hit by one transcription factor but is hit by %d",
                                                    treatment,length(all_tfs[filter_ind])))
    all_tfs[filter_ind]
}

#' @export
get_complementary_names <- function(pattern, group_names){
    complementary_names <- grep(pattern, group_names, invert = TRUE, value = TRUE)
    paste0(complementary_names, collapse = "|")
}

#' @export
plotDataFromActivity <- function(matrices_list, tf,
                                   labels,
                                   positive_elements_label,
                                   negative_elements_label,
                                 GeneExpressionMatrix,
                                 experimental_treatment,
                                 title = "",
                                 ...){
    colors <- c(Epiregulon = "red", FigR = "black", Pando = "blue", GRaNIE = "green1",
                cellOracle = "purple", scenic_plus_neg_gene = "grey", scenic_plus_pos_gene = "orange",
                scenic_plus_pos_region = "green4", scenic_plus_neg_region = "tan4",
                "Gene expression" = "darkslategray2")
    pos_ind <-  grep(paste(positive_elements_label, collapse = "|", sep =""), labels)
    neg_ind <-  grep(paste(negative_elements_label, collapse = "|", sep =""), labels)
    plot_data <- data.frame()
    matrices_list <- c(matrices_list, list("Gene expression" = assay(GeneExpressionMatrix, "normalizedCounts")[as.character(tf),,drop=FALSE]))
    colnames(matrices_list[[length(matrices_list)]]) <- colnames(GeneExpressionMatrix)
    for(i in seq_along(matrices_list))
    {
        if(nrow(matrices_list[[i]])==0) next
        if(!tf %in% rownames(matrices_list[[i]])) next
        activity_values <- matrices_list[[i]][as.character(tf),colnames(GeneExpressionMatrix)]
        res <- calculate_accuracy_metrics(activity_values, pos_ind, neg_ind)
        partial_data <- data.frame(FPR = res$FPR, TPR = res$TPR)
        partial_data$tf <- tf
        partial_data$package <- names(matrices_list)[i]
        if(!is.null(experimental_treatment)) res$treatment <- experimental_treatment
        plot_data <- rbind(plot_data, partial_data)
    }
    plot_data$package <- factor(plot_data$package, levels = c("Epiregulon", setdiff(sort(unique(plot_data$package)), "Epiregulon")))
    packages <- levels(plot_data$package)
    colors <- colors[packages]
    AUC_data <- plot_data %>% group_by(package) %>%
        group_map(function(x,y) cbind(y, data.frame(AUC = calculate_AUC(x$FPR, x$TPR))))
    AUC_data <- do.call(rbind, AUC_data)
    AUC_data$x <- 0.85
    AUC_data$y <- seq(0.5, 0.55-(nrow(AUC_data)*0.05), by =-0.05)
    AUC_data$AUC <- format(round(AUC_data$AUC,3), nsmall=3)
    library(ggplot2)
    ggp <- ggplot(data = plot_data, aes(x= FPR, y=TPR, color = package))+
        geom_line()+
        labs(title = substitute(title))+
        scale_color_manual(values = colors)+
        theme(panel.border = element_rect(color = "black", fill = NA),
              panel.background = element_rect(fill = "white", colour="white"),
              plot.title = element_text(hjust = 0.5))+
        geom_text(aes(x=0.85, y=0.55), label="AUC:", col="black")+
        geom_text(data=AUC_data, aes(x=x, y=y, label=AUC), col=colors, show.legend = FALSE)
    return(list(ggp,plot_data))
}
