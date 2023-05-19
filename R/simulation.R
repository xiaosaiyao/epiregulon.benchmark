# take simulation output and build reuglon
build_regulon <- function(sim_res){
    re_target <- reshape2::melt(sim_res$region_to_gene)
    re_target <- re_target[re_target[,3]>0,1:2]
    colnames(re_target) <- c("idxATAC", "target")
    tf_tg <- reshape2::melt(sim_res$.grn$geff)
    tf_tg <- tf_tg[tf_tg[,3]!=0,1:2]
    colnames(tf_tg) <- c("target", "tf")
    regulon <- merge(re_target, tf_tg)
    regulon <- lapply(regulon, as.character)
    regulon <- do.call(cbind,regulon) |> as.data.frame()
    regulon[,c("tf", "idxATAC", "target")]
}

# return objects needed for epiregulon flow
#' @export
processSimResults <- function(sim_res){
    regulon <- build_regulon(sim_res)
    dimnames(sim_res$counts_obs) <- dimnames(sim_res$counts)
    peakMatrix_obs <- SingleCellExperiment::SingleCellExperiment(list(peak = sim_res$atacseq_obs),
                                                                 colData = DataFrame(label=sim_res$cell_meta$pop),
                                                                 rowData=DataFrame(idxATAC=seq_len(nrow(sim_res$atacseq_obs))))
    geneExpMatrix_obs <- SingleCellExperiment::SingleCellExperiment(list(counts = sim_res$counts_obs),
                                                                    colData = DataFrame(label=sim_res$cell_meta$pop),
                                                                    rowData=DataFrame(gene=seq_len(nrow(sim_res$counts_obs))))
    peakMatrix <- SingleCellExperiment::SingleCellExperiment(list(peak = sim_res$atacseq_data),
                                                             colData = DataFrame(label=sim_res$cell_meta$pop),
                                                             rowData=DataFrame(idxATAC=seq_len(nrow(sim_res$atacseq_data))))
    geneExpMatrix <- SingleCellExperiment::SingleCellExperiment(list(counts = sim_res$counts),
                                                                colData = DataFrame(label=sim_res$cell_meta$pop),
                                                                rowData=DataFrame(gene=seq_len(nrow(sim_res$counts))))

    rownames(geneExpMatrix) <- seq_len(nrow(geneExpMatrix))
    rownames(peakMatrix) <- seq_len(nrow(peakMatrix))
    list(regulon = regulon, peakMatrix = peakMatrix, geneExpMatrix = geneExpMatrix,
         peakMatrix_obs = peakMatrix_obs, geneExpMatrix_obs = geneExpMatrix_obs)
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
assessActivityAccuracy <- function(GRN, sim_options, basic_sim_res, activities_obs) {
    correct_ranks <- list()
    interchanges <- c()
    tfs <- unique(GRN[,2])
    for(i in seq_along(tfs)){
        GRN_removal <- GRN
        # set weights for a given transcription factor to the values close to 0
        GRN_removal[GRN_removal$regulator.gene==tfs[i],"regulator.effect"] = 0.00001
        sim_options$GRN <- GRN_removal
        sim_res <- scMultiSim::sim_true_counts(sim_options)
        targets <- unique(GRN[GRN[,2] == tfs[i],1])
        counts_diff <- basic_sim_res$counts[targets,,drop = FALSE] - sim_res$counts[targets,,drop = FALSE]
        ranks_true <- rank(colSums(counts_diff), ties.method = "random")
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
