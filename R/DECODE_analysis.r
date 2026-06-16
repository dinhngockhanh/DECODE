#' Infer the probability for each mutation to belong to a SFS component
#' as inferred by DECODE.
#'
#' @param DECODE_result Object returned by function \code{\link{DECODE}}.
#' @param mutation_table Mutational dataframe.
#' Each row corresponds to a mutation, which can be associated with different copy number states
#' (unlike the mutation table for \code{\link{DECODE}}, which should include only mutations with the same copy number background).
#' \code{mutation_table} must contain column \code{normalized_VAF},
#' consisting of VAFs normalized to match the copy number state analyzed by \code{\link{DECODE}}.
#' @param mode String to specify which DECODE mode to use:
#' \code{"inference_A"}, \code{"inference_B"}, or \code{"validation"} (\code{"inference_A"} by default).
#' @param with_tail Logical variable for whether to apply DECODE fit with or without tail:
#' \code{TRUE}, \code{FALSE}, or \code{NA} (\code{NA} by default).
#' If \code{with_tail = NA}, the best overall fit chosen by \code{\link{DECODE}} is used.
#' @param N_clusters Integer to select the cluster count in the DECODE fit to apply (\code{NULL} by default).
#' If \code{N_clusters = NULL}, the best overall fit chosen by \code{\link{DECODE}} is used.
#' @return Dataframe extending \code{mutation_table} with probability columns:
#' \code{prob_tail} = probability that each mutation belongs to the tail (if a DECODE fit with tail is selected).
#' \code{prob_cluster_*} = probability that each mutation belongs to cluster \code{*}.
#' @export
DECODE_mutation_assignment <- function(DECODE_result,
                                       mutation_table,
                                       mode = "inference_A",
                                       with_tail = NA,
                                       N_clusters = NULL) {
    suppressPackageStartupMessages(library(crayon))
    if (!"normalized_VAF" %in% colnames(mutation_table)) {
        stop("mutation_table must contain 'normalized_VAF' column")
    }
    if (!mode %in% c("inference_A", "inference_B", "validation")) {
        stop("Invalid mode. Must be one of: 'inference_A', 'inference_B', 'validation'")
    }
    #---Determine which fit configuration to use
    if (is.na(with_tail)) {
        with_tail <- DECODE_result$best_with_tail
    }
    if (is.null(N_clusters)) {
        if (with_tail) {
            N_clusters <- DECODE_result$fits_with_tail$best_N_clusters
            fit_results <- DECODE_result$fits_with_tail$all_fits[[paste0(N_clusters, "_clusters")]]
        } else {
            N_clusters <- DECODE_result$fits_without_tail$best_N_clusters
            fit_results <- DECODE_result$fits_without_tail$all_fits[[paste0(N_clusters, "_clusters")]]
        }
    } else {
        if (with_tail) {
            fit_results <- DECODE_result$fits_with_tail$all_fits[[paste0(N_clusters, "_clusters")]]
        } else {
            fit_results <- DECODE_result$fits_without_tail$all_fits[[paste0(N_clusters, "_clusters")]]
        }
    }
    report <- paste0("\n", bold(red("Assign mutations to DECODE components with configuration ")), bold(yellow(paste0(ifelse(with_tail, "with tail", "without tail"), " + ", N_clusters, " clusters"))), bold(red("...")), "\n")
    cat(report)
    #---Extract DECODE components (each as a length-B vector)
    sfs_components <- list()
    if (with_tail) {
        sfs_components[["tail"]] <- colMeans(fit_results[[paste0("SFS_", mode, "_tail")]], na.rm = TRUE)
    }
    for (i in 1:N_clusters) {
        sfs_components[[paste0("cluster_", i)]] <- colMeans(fit_results[[paste0("SFS_", mode, "_cluster_", i)]], na.rm = TRUE)
    }
    component_IDs <- names(sfs_components)
    prob_cols <- paste0("prob_", component_IDs)
    #---Build (B x C) component proportions matrix and row-normalize in one shot
    P <- do.call(cbind, sfs_components)
    if (!is.matrix(P)) {
        P <- matrix(P,
            ncol = length(sfs_components),
            dimnames = list(NULL, component_IDs)
        )
    }
    row_sums <- rowSums(P, na.rm = TRUE)
    P <- P / row_sums
    P[!is.finite(P)] <- 0
    #---Vectorized bin assignment: preserves original (lo, hi] semantics via left.open = TRUE
    bin_edges <- c(-Inf, DECODE_result$SFS_frequencies[1:(length(DECODE_result$SFS_frequencies) - 1)], Inf)
    B <- length(bin_edges) - 1L
    vaf <- mutation_table$normalized_VAF
    bin_idx <- findInterval(vaf, bin_edges, left.open = TRUE)
    valid <- !is.na(bin_idx) & bin_idx >= 1L & bin_idx <= B
    prob_mat <- matrix(NA_real_, nrow = length(vaf), ncol = length(component_IDs))
    if (any(valid)) {
        prob_mat[valid, ] <- P[bin_idx[valid], , drop = FALSE]
    }
    #---Append probability columns in one pass (one column-copy per component, not per bin)
    for (j in seq_along(prob_cols)) {
        mutation_table[[prob_cols[j]]] <- prob_mat[, j]
    }
    #---Return mutation table extended with DECODE component probabilities
    return(mutation_table)
}
