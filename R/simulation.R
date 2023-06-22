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
    regulon <- do.call(cbind,regulon) |> as.data.frame()
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
processSimResults <- function(sim_res){
    regulon <- build_regulon(sim_res)
    dimnames(sim_res$counts_obs) <- dimnames(sim_res$counts)
    peakMatrix <- SingleCellExperiment(list(peak = sim_res$atacseq_data, peak_obs = sim_res$atacseq_obs),
                                       colData = DataFrame(label=sim_res$cell_meta$pop),
                                       rowData=DataFrame(idxATAC=seq_len(nrow(sim_res$atacseq_data))))
    
    norm_counts <- normalize_counts(sim_res$counts)
    logcounts <- log2(norm_counts+1)
    norm_counts_obs <- normalize_counts(sim_res$counts_obs)
    logcounts_obs <- log2(norm_counts_obs+1)
    geneExpMatrix <- SingleCellExperiment(list(counts = sim_res$counts, counts_obs = sim_res$counts_obs,
                                               logcounts = logcounts, logcounts_obs = logcounts_obs,
                                               norm_counts = norm_counts, norm_counts_obs = norm_counts_obs),
                                          colData = DataFrame(label=sim_res$cell_meta$pop),
                                          rowData=DataFrame(gene=seq_len(nrow(sim_res$counts))))
    rownames(geneExpMatrix) <- seq_len(nrow(geneExpMatrix))
    rownames(peakMatrix) <- seq_len(nrow(peakMatrix))
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
                     batch_size = 50, sim_options, tfs = NULL){
    if(is.null(tfs)) tfs <- unique(sim_options$GRN[,2])
    remaining_tfs <- tfs
    tf_removal <- list()
    while(length(remaining_tfs) != 0){
        chosen_tfs <- sample(remaining_tfs, min(batch_size, length(remaining_tfs)))
        tryCatch({tf_removal <- c(tf_removal, BiocParallel::bplapply(X = chosen_tfs,
                                                                     FUN = tf_effect_bp,
                                                                     sim_options = sim_options,
                                                                     BPPARAM = BPPARAM))
        tf_removal <- tf_removal[!sapply(tf_removal, function(x) is.null(x$counts))]
        remaining_tfs <- setdiff(remaining_tfs, chosen_tfs)
        message(Sys.time())
        message(sprintf("Well done. Number of analysed tfs %d", length(tf_removal)))
        },
        error = function(cond) {
            message(Sys.time())
            message(sprintf("Error. Number of analysed tfs %d", length(tf_removal)))
            message(cond)},
        finally = {remaining_tfs <- setdiff(remaining_tfs, chosen_tfs)})
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

#' @import BiocParallel
#' @export
assessActivityAccuracy <- function(sim_options, activities_obs,
                                   seed = 1010, true_ranks_matrix) {
    correct_ranks <- list()
    interchanges <- c()
    tfs <- as.character(unique(sim_options$GRN[,2]))
    set.seed(seed)
    for(i in seq_along(tfs)){
        ranks_true <- true_ranks_matrix[tfs[i],]
        if (tfs[i] %in% rownames(activities_obs))
            ranks_obs <- rank(activities_obs[as.character(tfs[i]),], ties.method = "random")
        # all activities equal to 0 so the rank order is random
        else
            ranks_obs <- sample(1:ncol(activities_obs))
        # use observed ranks in the true order
        interchanges <- c(interchanges, calculate_interchanges(ranks_obs[order(ranks_true)]))
        # use true ranks in the observed order
        correct_ranks[[i]] <- ranks_true[order(ranks_obs)]
    }
    names(interchanges) <- names(correct_ranks) <- tfs
    list(interchange_number = interchanges, correct_ranks = correct_ranks)
}

#' @import scMultiSim
tf_effect_bp <- function(tf, sim_options){
    # set weights for a given transcription factor to the values close to 0
    sim_options$GRN[sim_options$GRN$regulator.gene==tf,"regulator.effect"] = 1e-10
    counts <- tryCatch({
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
    while(nrow(false_connections) < false_rows_n){
        tf <- sample(regulon$tf,1)
        re <- sample(regulon$idxATAC,1)
        target <- sample(regulon$target,1)
        new_false_row <- data.frame(tf = tf, idxATAC = re, target = target)
        if(tail(duplicated(rbind(regulon, new_false_row)), 1)) next
        false_connections <- rbind(false_connections, new_false_row)
    }
    regulon$connection_type <- TRUE
    false_connections$connection_type <- FALSE
    rbind(regulon, false_connections)
}

#' @export
getActivity <- function(regulon, geneExpMatrix, peakMatrix,
                        weightMethods = c("wilcoxon"),
                        clusters_list = list(), ...){
    res_list <- list()
    for(method in weightMethods){
        clusters <- clusters_list[[method]]
        regulon.w <- addWeights(regulon,
                                peakMatrix = peakMatrix,
                                expMatrix = geneExpMatrix,
                                method = method,
                                clusters = clusters,
                                ...)

        exp_assay <- ifelse(is.null(list(...)$exp_assay), "counts", list(...)$exp_assay)
        score.combine <- calculateActivity(regulon = regulon.w,
                                           mode = "weight",
                                           method = "weightedMean",
                                           expMatrix = geneExpMatrix,
                                           exp_assay = exp_assay)
        res_list[[length(res_list)+1]] <- score.combine
        names(res_list)[length(res_list)] <- method
    }
    res_list
}
