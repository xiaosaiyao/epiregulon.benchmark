# take simulation output and build reuglon
#' @import reshape2
build_regulon <- function(sim_res){
    re_target <- reshape2::melt(sim_res$region_to_gene)
    re_target <- re_target[re_target[,3]>0,1:2]
    colnames(re_target) <- c("idxATAC", "target")
    tf_tg <- reshape2::melt(sim_res$.grn$geff)
    tf_tg <- tf_tg[tf_tg[,3]!=0,1:2]
    colnames(tf_tg) <- c("target", "tf")
    regulon <- merge(re_target, tf_tg)
    regulon[,c("tf", "target")] <- lapply(regulon[,c("tf", "target")] , as.character)
    regulon[,c("tf", "idxATAC", "target")]
}

#' @import Matrix
normalize_counts <- function(count_matrix){
    library_sizes <- Matrix::colSums(count_matrix)
    t(apply(count_matrix,1,function(x) x/library_sizes)*1e4)
}

# return objects needed for epiregulon flow
#' @import SingleCellExperiment
#' @export
processSimResults <- function(sim_res, seed=23143, technical_noise = TRUE){
    regulon <- build_regulon(sim_res)
    dimnames(sim_res$counts_obs$counts) <- dimnames(sim_res$counts)
    if(technical_noise){
        peakMatrix <- SingleCellExperiment(list(peak = sim_res$atac_counts, peak_obs = sim_res$atacseq_obs),
                                           colData = DataFrame(label=sim_res$cell_meta$pop),
                                           rowData=DataFrame(idxATAC=seq_len(nrow(sim_res$atacseq_data))))
        logcounts <- log2(sim_res$counts+1)
        logcounts_obs <- log2(sim_res$counts_obs$counts+1)
        norm_counts <- normalize_counts(sim_res$counts)
        lognorm_counts <- log2(norm_counts+1)
        norm_counts_obs <- normalize_counts(sim_res$counts_obs$counts)
        lognorm_counts_obs <- log2(norm_counts_obs+1)
        geneExpMatrix <- SingleCellExperiment(list(counts = sim_res$counts, counts_obs = sim_res$counts_obs$counts,
                                                   lognorm_counts = lognorm_counts, lognorm_counts_obs = lognorm_counts_obs,
                                                   norm_counts = norm_counts, norm_counts_obs = norm_counts_obs,
                                                   logcounts = logcounts, logcounts_obs = logcounts_obs),
                                              colData = DataFrame(label=sim_res$cell_meta$pop),
                                              rowData=DataFrame(gene=seq_len(nrow(sim_res$counts))))
        rownames(geneExpMatrix) <- seq_len(nrow(geneExpMatrix))
        rownames(peakMatrix) <- seq_len(nrow(peakMatrix))
        set.seed(seed)
        geneExpMatrix <- scran::fixedPCA(geneExpMatrix, name="IterativeLSI_ATAC", subset.row=NULL) # use this name to be consistent with the default epiregulon settings
        rD <- reducedDim(geneExpMatrix)
        rownames(rD) <- colnames(peakMatrix)
        reducedDim(peakMatrix) <- rD
        reducedDimNames(peakMatrix) <- "UMAP_ATAC"
    }
    else{
        peakMatrix <- SingleCellExperiment(list(peak = sim_res$atac_counts),
                                           colData = DataFrame(label=sim_res$cell_meta$pop),
                                           rowData=DataFrame(idxATAC=seq_len(nrow(sim_res$atacseq_data))))
        logcounts <- log2(sim_res$counts+1)
        norm_counts <- normalize_counts(sim_res$counts)
        lognorm_counts <- log2(norm_counts+1)
        geneExpMatrix <- SingleCellExperiment(list(counts = sim_res$counts,
                                                   lognorm_counts = lognorm_counts,
                                                   norm_counts = norm_counts,
                                                   logcounts = logcounts),
                                              colData = DataFrame(label=sim_res$cell_meta$pop),
                                              rowData=DataFrame(gene=seq_len(nrow(sim_res$counts))))
        rownames(geneExpMatrix) <- seq_len(nrow(geneExpMatrix))
        rownames(peakMatrix) <- seq_len(nrow(peakMatrix))
        set.seed(seed)
        geneExpMatrix <- scran::fixedPCA(geneExpMatrix, name="IterativeLSI_ATAC", subset.row=NULL) # use this name to be consistent with the default epiregulon settings
        rD <- reducedDim(geneExpMatrix)
        rownames(rD) <- colnames(peakMatrix)
        reducedDim(peakMatrix) <- rD
        reducedDimNames(peakMatrix) <- "UMAP_ATAC"
    }

    list(regulon = regulon, peakMatrix = peakMatrix, geneExpMatrix = geneExpMatrix)
}

calculate_interchanges <- function(vec){
    vec <- rank(vec, ties.method = "random")
    interchanges <- 0
    for(i in (seq_len(length(vec)-1))){
        current.pos = which(vec==i)
        interchanges <- interchanges + current.pos - 1
        vec <- vec[vec!= i]
    }
    interchanges
}

#' @export
removeTF <- function(basic_sim_res, BPPARAM = BiocParallel::MulticoreParam(),
                     batch_size = 50, sim_options, tfs = NULL, seed=10110){
    if(is.null(tfs)) tfs <- unique(sim_options$GRN[,2])
    remaining_tfs <- tfs
    tf_removal <- list()
    while(length(remaining_tfs) != 0){
        chosen_tfs <- unlist(sample(list(remaining_tfs), min(batch_size, length(remaining_tfs))))
        tryCatch({tf_removal <- c(tf_removal, BiocParallel::bplapply(X = chosen_tfs,
                                                                     FUN = tf_effect_bp,
                                                                     sim_options = sim_options,
                                                                     seed = seed,
                                                                     BPPARAM = BPPARAM))
        tf_removal <- tf_removal[!sapply(tf_removal, function(x) is.null(x$counts))]
        message(Sys.time())
        message(sprintf("Well done. Number of analysed tfs %d", length(tf_removal)))
        message(names(tf_removal))
        },
        error = function(cond) {
            message(Sys.time())
            message(sprintf("Error. Number of analysed tfs %d", length(tf_removal)))
            message(cond)},
        finally = {
            message(sprintf("Remaining tfs: %s", remaining_tfs))
            message(sprintf("chosen tfs: %s", chosen_tfs))
            remaining_tfs <- setdiff(remaining_tfs, chosen_tfs)
            message(sprintf("Remaining tfs after setdiff: %s", remaining_tfs))})
    }
    message(Sys.time())
    names(tf_removal) <- sapply(tf_removal, function(x) x$tf)
    lapply(tf_removal, function(x) x$counts)
}

#' @export
calculateTrueActivityRanks <- function(basic_sim_res, counts_removal, seed = 1010,
                                       sim_options, transcription_effect = FALSE){
    set.seed(seed)
    true_ranks_matrix <- matrix(nrow = length(counts_removal), ncol = sim_options$num.cells,
                                dimnames = list(names(counts_removal), colnames(basic_sim_res$counts)))
    for(i in seq_along(counts_removal)){
        tf <- names(counts_removal)[i]
        targets <- unique(sim_options$GRN[sim_options$GRN[,2] == tf, 1])
        counts_diff <- basic_sim_res$counts[targets,,drop = FALSE] - counts_removal[[i]][targets,,drop = FALSE]
        if(transcription_effect) true_ranks_matrix[tf,] <- colSums(counts_diff)
        else true_ranks_matrix[tf,] <- rank(colSums(counts_diff), ties.method = "random")
    }
    true_ranks_matrix
}

#' @export
assessActivityAccuracy <- function(activities_obs, transcription_difference_matrix) {
    if(is.null(activities_obs)) return(NULL)
    correlations <- unlist(lapply(rownames(activities_obs), function(tf, m1, m2) cor(m1[tf, ], m2[tf, ]),
                           activities_obs, transcription_difference_matrix))
    names(correlations) <- rownames(activities_obs)
    correlations
}

#' @import scMultiSim
tf_effect_bp <- function(tf, sim_options, seed){
    # set weights for a given transcription factor to the values close to 0
    sim_options$GRN[sim_options$GRN$regulator.gene==tf,"regulator.effect"] = 1e-10
    counts <- tryCatch({
        set.seed(seed)
        sim_res <- scMultiSim::sim_true_counts(sim_options)
        sim_res <- as.list(sim_res)
        sim_res$counts
    },
    error = function(cond){
        message(cond)
        NULL
    })
    list(tf = tf, counts = counts)
}

#' @param x integer vector corresponding to the true ranks arranged according to the expected ranks
#' @export
plot_recovery_curve <- function(x, ...){
    x <- mapply(x, 1:length(x), FUN=max)
    x <- table(c(1:length(x), x))
    x <- x - 1
    y <- cumsum(x)
    x <- 0:(length(y)-1)
    x <- c(rep(x, each =2), length(y))
    y <- c(0, rep(y, each = 2))
    x <- x/max(x)
    y <- y/max(y)
    plot(x, y, xlim = c(0,1), ylim = c(0,1), xlab = "Rank cutoff", ylab = "Proportion of retrieved elements",
         type = "l", ...)
}

#' @export
addFalseConnections <- function(regulon, fraction_false = 0.5, seed = 10010){
    set.seed(seed)
    if(fraction_false<=0 | fraction_false>=1) stop("franction_false argument shoubd be a number from the open set (0,1)")
    false_rows_n <- round(nrow(regulon)/(1-fraction_false)) - nrow(regulon)
    false_connections <- data.frame()
    while((false_rows_n - nrow(false_connections)) > 0){
        to_be_added <- false_rows_n - nrow(false_connections)
        tf <- sample(regulon$tf, to_be_added, replace = TRUE)
        re <- sample(regulon$idxATAC, to_be_added, replace = TRUE)
        target <- sample(regulon$target, to_be_added, replace = TRUE)
        new_false_connections <- data.frame(tf = tf, idxATAC = re, target = target)
        duplicated_rows <- rev(duplicated(rbind(false_connections, regulon, new_false_connections)))[1:nrow(new_false_connections)]
        new_false_connections <- new_false_connections[!rev(duplicated_rows),]
        false_connections <- rbind(false_connections, new_false_connections)
    }
    regulon$connection_type <- TRUE
    false_connections$connection_type <- FALSE
    if (any(duplicated(rbind(regulon, false_connections)))) stop("Duplicated rows in the regulon")
    rbind(regulon, false_connections)
}

#' @export
addWeightsVariousMethods <- function(regulon, input_objects, exp_assay = "norm_counts_obs",
                                     peak_assay = "peak_obs",
                                     weightMethods =  c("wilcoxon", "corr", "MI"), ...){
     clusters_list = list("corr" = input_objects$geneExpMatrix$label,
                          "MI" = input_objects$geneExpMatrix$label)
    regulon_list <- list()
    geneExpMatrix <- input_objects$geneExpMatrix
    peakMatrix <- input_objects$peakMatrix
     for(method in weightMethods){
         clusters <- clusters_list[[method]]
         regulon.w <- addWeights(regulon,
                                 peakMatrix = peakMatrix,
                                 expMatrix = geneExpMatrix,
                                 method = method,
                                 clusters = clusters,
                                 peak_assay = peak_assay,
                                 exp_assay = exp_assay,
                                 ...)
         if(is.null(regulon.w)) next
         regulon_list <- setNames(c(regulon_list, list(regulon.w)), c(names(regulon_list), method))
     }
    regulon_list
}

#' @export
getActivity <- function(regulon_list, input_objects, exp_assay = "norm_counts_obs", ...){
    res <- list()
    for(method in names(regulon_list)){
        score.combine <- calculateActivity(regulon = regulon_list[[method]],
                                           mode = "weight",
                                           method = "weightedMean",
                                           expMatrix = input_objects$geneExpMatrix,
                                           exp_assay = exp_assay,
                                           ...)
        res <- setNames(c(res, list(score.combine)), c(names(res), method))
    }
    res
}

#' @export
accuracyComparisonPruning <- function(regulon, input_objects,
                                      transcription_difference_matrix,
                                      regulon_cutoffs = c(2, 0.05, 0.001, 0.0001),
                                      weights_exp_cutoff = NULL,
                                      weights_peak_cutoff = NULL,
                                      peak_assay = "peak",
                                      exp_assay = "norm_counts",
                                      only_clusters = FALSE,
                                      aggregateCells = FALSE,
                                      tf_re.merge = TRUE,
                                      TF_expression = FALSE,
                                      return_intermediates = FALSE,
                                      ...){
    regulon <-  pruneRegulon(regulon = regulon,
                             expMatrix = input_objects$geneExpMatrix,
                             exp_assay = exp_assay,
                             peakMatrix = input_objects$peakMatrix,
                             peak_assay = peak_assay,
                             test = "chi.sq",
                             clusters = input_objects$geneExpMatrix$label,
                             regulon_cutoff = 2,
                             ...)
    regulon_list <- addWeightsVariousMethods(regulon, input_objects,
                                             exp_cutoff = weights_exp_cutoff,
                                             peak_cutoff = weights_peak_cutoff,
                                             exp_assay = exp_assay,
                                             peak_assay = peak_assay,
                                             aggregateCells = aggregateCells,
                                             tf_re.merge = tf_re.merge)
    df <- data.frame()
        for(p_val_cutoff in regulon_cutoffs){
            regulon_list_pruned <- lapply(regulon_list, function(x) filterRegulonPVal(x, p_val_cutoff, only_clusters = only_clusters))
            activity_list <- getActivity(regulon_list_pruned, input_objects, exp_assay = exp_assay)
            for(i in seq_along(activity_list)){
                method_res <- assessActivityAccuracy(activities_obs = activity_list[[i]],
                                                     transcription_difference_matrix = transcription_difference_matrix)
                # account for no tf being recovered
                if(is.null(method_res)) next
                df_part <- data.frame(tf = names(method_res),
                                      method = names(activity_list)[i],
                                      correlations = method_res)
                df_part$n_tf <- nrow(activity_list[[i]])
                df_part$pval_cutoff <- p_val_cutoff
                df <- rbind(df, df_part)
            }
        }
        if(TF_expression){
            tf_expression_matrix <- assay(input_objects$geneExpMatrix, exp_assay)[rownames(transcription_difference_matrix),]
            method_res <- assessActivityAccuracy(activities_obs = tf_expression_matrix,
                                                 transcription_difference_matrix = transcription_difference_matrix)
            # account for no tf being recovered
            if(is.null(method_res)) next
            df_part <- data.frame(tf = names(method_res),
                                  method = "TF expression",
                                  correlations = method_res)
            df_part$n_tf <- nrow(transcription_difference_matrix)
            df_part$pval_cutoff <- NA
            df <- rbind(df, df_part)

        }
    if(return_intermediates) return(list(df = df, regulon_list = regulon_list))
    df
}

filterRegulonPVal <- function(regulon, cutoff, only_clusters = FALSE){
    if(!"pval" %in% colnames(regulon)) return(regulon)
    if(is.na(cutoff)) return(regulon)
    if(only_clusters){
        min_pval <- apply(regulon$pval, 1, function (x){
            if (sum(is.na(x[2:length(x)])) == length(x)-1)
                1
            else
                min(x[2:length(x)], na.rm = TRUE)
        })
        regulon <- regulon[min_pval < cutoff, ]
    }
    else{
        min_pval <- apply(regulon$pval, 1, function (x){
            if (sum(is.na(x)) == length(x))
                1
            else
                min(x, na.rm = TRUE)
        })
        regulon <- regulon[min_pval < cutoff, ]
    }
    regulon
}


#' @export
runEpiregulonWorkflows <- function(regulon, input_objects,
                                      transcription_difference_matrix,
                                      regulon_cutoffs = NULL,
                                      weights_exp_cutoff = NULL,
                                      weights_peak_cutoff = NULL,
                                      peak_assay = "peak",
                                      exp_assay = "norm_counts",
                                      only_clusters = FALSE,
                                      aggregateCells = FALSE,
                                      tf_re.merge = TRUE,
                                      TF_expression = FALSE,
                                      return_intermediates = FALSE,
                                      ...){
    if(!is.null(regulon_cutoffs)){
        regulon <-  pruneRegulon(regulon = regulon,
                                 expMatrix = input_objects$geneExpMatrix,
                                 exp_assay = exp_assay,
                                 peakMatrix = input_objects$peakMatrix,
                                 peak_assay = peak_assay,
                                 test = "chi.sq",
                                 clusters = input_objects$geneExpMatrix$label,
                                 regulon_cutoff = 2,
                                 ...)
    }
    else regulon_cutoffs <- NA
    regulon_list <- addWeightsVariousMethods(regulon, input_objects,
                                             exp_cutoff = weights_exp_cutoff,
                                             peak_cutoff = weights_peak_cutoff,
                                             exp_assay = exp_assay,
                                             peak_assay = peak_assay,
                                             aggregateCells = aggregateCells,
                                             tf_re.merge = tf_re.merge)
    df <- data.frame()
    for(p_val_cutoff in regulon_cutoffs){
        regulon_list_pruned <- lapply(regulon_list, function(x) filterRegulonPVal(x, p_val_cutoff, only_clusters = only_clusters))
        activity_list <- getActivity(regulon_list_pruned, input_objects, exp_assay = exp_assay)
        for(i in seq_along(activity_list)){
            method_res <- assessActivityAccuracy(activities_obs = activity_list[[i]],
                                                 transcription_difference_matrix = transcription_difference_matrix)
            # account for no tf being recovered
            if(is.null(method_res)) next
            df_part <- data.frame(tf = names(method_res),
                                  method = names(activity_list)[i],
                                  correlations = method_res)
            df_part$n_tf <- nrow(activity_list[[i]])
            df_part$pval_cutoff <- p_val_cutoff
            df <- rbind(df, df_part)
        }
    }
    if(TF_expression){
        tf_expression_matrix <- assay(input_objects$geneExpMatrix, exp_assay)[rownames(transcription_difference_matrix),]
        method_res <- assessActivityAccuracy(activities_obs = tf_expression_matrix,
                                             transcription_difference_matrix = transcription_difference_matrix)
        # account for no tf being recovered
        if(is.null(method_res)) next
        df_part <- data.frame(tf = names(method_res),
                              method = "TF expression",
                              correlations = method_res)
        df_part$n_tf <- nrow(transcription_difference_matrix)
        df_part$pval_cutoff <- NA
        df <- rbind(df, df_part)

    }
    if(return_intermediates) return(list(df = df, regulon_list = regulon_list))
    df
}


#' @export
check_pruning_effect <- function(df){
    df_new <- data.frame()
    reference_all <- df[is.na(df$pval_cutoff),]
    df <- df[!is.na(df$pval_cutoff), ]
    df <- split(df, df$pval_cutoff)
    for(i in 1:length(df)){
        reference <- reference_all[reference_all$tf %in% df[[i]]$tf,]
        reference <- reference[match(df[[i]]$tf, reference$tf),]
        df[[i]]$correlation_delta <- df[[i]]$correlations - reference$correlations
        df[[i]]$n_tf <- length(unique(df[[i]]$tf))
        df_new <- rbind(df_new, df[[i]])
    }
    df_new$pval_cutoff <- factor(df_new$pval_cutoff, levels = c(0.05, 0.001, 0.0001))
    ggplot(df_new, aes(x = pval_cutoff, y = correlation_delta))+
        geom_boxplot(lwd= 0.1, outlier.size = 0.1)+
        #stat_summary(fun = median, geom = "label")+
        coord_flip()+
        facet_wrap(~method, ncol = 1)+
        xlab("p-value cutoff")+
        ylab("Correlation change")+
        ggtitle("Effect of the pruning on the activity assessment accuracy")+
        theme(plot.title = element_text(hjust = 0.5),
              panel.background = element_rect(fill = "white", colour="white"))+
        geom_hline(yintercept = 0, linetype = "dashed", lwd = 0.2)
}

#' @export
plot_accuracy_curves <- function(regulon.mixed_fixed){
    RP_table <- data.frame()
    for(cutoff in seq(0,3, by = 0.01)){
        TP_fixed <- sum((abs(regulon.mixed_fixed$stats[,"all"]) >= cutoff) & (regulon.mixed_fixed$connection_type == TRUE))
        FP_fixed <- sum((abs(regulon.mixed_fixed$stats[,"all"]) >= cutoff) & (regulon.mixed_fixed$connection_type == FALSE))
        FN_fixed <- sum((abs(regulon.mixed_fixed$stats[,"all"]) < cutoff) & (regulon.mixed_fixed$connection_type == TRUE))
        RP_table <- rbind(RP_table, data.frame(cutoff = cutoff,
                                               TP_fixed = TP_fixed, FP_fixed = FP_fixed, FN_fixed = FN_fixed))
    }
    y1 <- (RP_table$TP_fixed/(RP_table$TP_fixed + RP_table$FN_fixed))
    y2 <- (RP_table$TP_fixed/(RP_table$TP_fixed + RP_table$FP_fixed))
    ylim <- c(min(y1,y2)-0.01, max(y1,y2)+0.01)
    pdf(save_paths[1], width  = 6, height = 6)
    plot(seq(0,3, by = 0.01), y1, xlab = "z-score cutoff (absolute value)",
         ylab = "Recall/precision", type = "l", main = "Accuracy of the pruning, observed counts", ylim=ylim)
    lines(seq(0,3, by = 0.01), y2, col = "red")
    legend(x = 0.2, y=1, c("Recall", "Precision"), col=c("red", "black"), lty = c(1,1))
}

#' @export
addTfTargetCorr <- function(regulon, geneExpMatrix, assay_name = "norm_counts_obs"){
    geneMatrix <- SummarizedExperiment::assay(geneExpMatrix, assay_name)
    regulon$tf_tg_corr <- unlist(lapply(seq_len(nrow(regulon)), function(i, geneMatrix) cor(geneMatrix[regulon$tf[i],], geneMatrix[regulon$target[i], ]), geneMatrix))
    regulon
}

#' @export
addTgNumb<- function(regulon, sim_res){
    n_tg <- table(sim_res$.grn$params[,"regulator"])
    regulon$tg_numb <- n_tg[as.character(regulon$tf)]
    regulon
}

#' @export
addReTargetCorr <- function(regulon, geneExpMatrix, peakMatrix, assay_name_gene = "norm_counts", assay_name_peak="peak"){
    geneMatrix <- SummarizedExperiment::assay(geneExpMatrix, assay_name_gene)
    peakMatrix <- SummarizedExperiment::assay(peakMatrix, assay_name_peak)
    regulon$re_tg_corr <- unlist(lapply(seq_len(nrow(regulon)), function(i, geneMatrix, peakMatrix) cor(peakMatrix[regulon$idxATAC[i],], geneMatrix[regulon$target[i], ]), geneMatrix, peakMatrix))
    regulon
}

