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

MOBSTER_summary_statistics <- function(mobster_df, MOBSTER_result) {
    mobster_df[nrow(mobster_df) + 1, "Sample"] <- MOBSTER_result$best$description
    mobster_df[nrow(mobster_df), "Nmut"] <- MOBSTER_result$best$N
    mobster_df[nrow(mobster_df), "Tail"] <- MOBSTER_result$best$fit.tail
    mobster_df[nrow(mobster_df), "Tail_Nmut"] <- MOBSTER_result$best$N.k[[1]]
    mobster_df[nrow(mobster_df), "Tail_power"] <- MOBSTER_result$best$shape + 1
    mobster_df[nrow(mobster_df), "Tail_Pareto_shape"] <- MOBSTER_result$best$shape
    mobster_df[nrow(mobster_df), "Tail_Pareto_scale"] <- MOBSTER_result$best$scale
    mobster_df[nrow(mobster_df), "Cluster_count"] <- MOBSTER_result$best$Kbeta
    for (k in 1:MOBSTER_result$best$Kbeta) {
        mobster_df[nrow(mobster_df), paste0("Cluster_Nmut_", k)] <- MOBSTER_result$best$N.k[[k + 1]]
        mobster_df[nrow(mobster_df), paste0("Cluster_VAF_", k)] <- MOBSTER_result$best$a[[k]] / (MOBSTER_result$best$a[[k]] + MOBSTER_result$best$b[[k]])
        mobster_df[nrow(mobster_df), paste0("Cluster_Beta_a_", k)] <- MOBSTER_result$best$a[[k]]
        mobster_df[nrow(mobster_df), paste0("Cluster_Beta_b_", k)] <- MOBSTER_result$best$b[[k]]
    }
    return(mobster_df)
}

DECODE_summary_statistics <- function(decode_df, DECODE_result) {
    parameters_inference_A <- DECODE_result$final_fit$best_fit$parameters_inference_A
    parameters_inference_B <- DECODE_result$final_fit$best_fit$parameters_inference_B
    if (DECODE_result$final_fit$best_fit$tail_status) {
        Tail_Nmut_inference_A <- parameters_inference_A[1] * sum(DECODE_result$SFS_data_inference_A)
        Tail_Nmut_inference_B <- parameters_inference_B[1] * sum(DECODE_result$SFS_data_inference_B)
        Tail_power <- parameters_inference_A[2]
        Cluster_Nmut_inference_A <- parameters_inference_A[seq(3, length(parameters_inference_A), 2)] * sum(DECODE_result$SFS_data_inference_A)
        Cluster_Nmut_inference_B <- parameters_inference_B[seq(3, length(parameters_inference_B), 2)] * sum(DECODE_result$SFS_data_inference_B)
        Cluster_VAF <- parameters_inference_A[seq(4, length(parameters_inference_A), 2)]
    } else {
        Tail_Nmut_inference_A <- 0
        Tail_Nmut_inference_B <- 0
        Tail_power <- NA
        Cluster_Nmut_inference_A <- parameters_inference_A[seq(1, length(parameters_inference_A), 2)] * sum(DECODE_result$SFS_data_inference_A)
        Cluster_Nmut_inference_B <- parameters_inference_B[seq(1, length(parameters_inference_B), 2)] * sum(DECODE_result$SFS_data_inference_B)
        Cluster_VAF <- parameters_inference_A[seq(2, length(parameters_inference_A), 2)]
    }
    decode_df[nrow(decode_df) + 1, "Sample"] <- sample
    decode_df[nrow(decode_df), "Nmut"] <- nrow(DECODE_result$mutational_table)
    decode_df[nrow(decode_df), "Nmut_inference_A"] <- sum(DECODE_result$SFS_data_inference_A)
    decode_df[nrow(decode_df), "Nmut_inference_B"] <- sum(DECODE_result$SFS_data_inference_B)
    decode_df[nrow(decode_df), "Tail"] <- DECODE_result$final_fit$best_fit$tail_status
    decode_df[nrow(decode_df), "Tail_Nmut_inference_A"] <- Tail_Nmut_inference_A
    decode_df[nrow(decode_df), "Tail_Nmut_inference_B"] <- Tail_Nmut_inference_B
    ####################################################################
    if (DECODE_result$final_fit$best_fit$tail_status) {
        decode_df[nrow(decode_df), "Tail_Nmut_predict"] <-
            Tail_Nmut_inference_A /
                rowSums(DECODE_result$final_fit$best_fit$component_distributions_inference_A$SFS_expected)[1] *
                rowSums(DECODE_result$final_fit$best_fit$component_distributions_inference_A$SFS_exact)[1]
    } else {
        decode_df[nrow(decode_df), "Tail_Nmut_predict"] <- NA
    }
    ####################################################################
    decode_df[nrow(decode_df), "Tail_power"] <- Tail_power
    decode_df[nrow(decode_df), "Cluster_count"] <- length(Cluster_VAF)
    for (k in 1:(length(Cluster_VAF))) {
        decode_df[nrow(decode_df), paste0("Cluster_Nmut_inference_A_", k)] <- Cluster_Nmut_inference_A[k]
        decode_df[nrow(decode_df), paste0("Cluster_Nmut_inference_B_", k)] <- Cluster_Nmut_inference_B[k]
        ################################################################
        tmp <- ifelse(DECODE_result$final_fit$best_fit$tail_status, 1, 0)
        decode_df[nrow(decode_df), paste0("Cluster_Nmut_predict_", k)] <-
            Cluster_Nmut_inference_A[k] /
                rowSums(DECODE_result$final_fit$best_fit$component_distributions_inference_A$SFS_expected)[k + tmp] *
                rowSums(DECODE_result$final_fit$best_fit$component_distributions_inference_A$SFS_exact)[k + tmp]
        ################################################################
        decode_df[nrow(decode_df), paste0("Cluster_VAF_", k)] <- Cluster_VAF[k]
    }
    return(decode_df)
}

plot_analysis <- function(results,
                          algorithm,
                          cohort,
                          algorithm_colors,
                          cluster_shapes,
                          cluster_labels,
                          cluster_colors,
                          cohort_colors,
                          folder_workplace) {
    library(dplyr)
    library(ggplot2)
    #---------------Plot distribution of cancer-specific mutation counts
    if ("Cancer_type" %in% colnames(results)) {
        plot_df <- results %>%
            group_by(Cancer_type) %>%
            mutate(median_Nmut = median(Nmut, na.rm = TRUE)) %>%
            ungroup() %>%
            arrange(desc(median_Nmut))
        plot_df$Cancer_type <- factor(plot_df$Cancer_type, levels = unique(plot_df$Cancer_type))
        p <- ggplot(plot_df, aes(x = Cancer_type, y = Nmut, fill = Cancer_type)) +
            geom_boxplot(color = "black", size = 2) +
            scale_y_log10() +
            scale_fill_manual(values = ICGC_cohort_colors) +
            xlab(NULL) +
            ylab("Mutation count") +
            labs(NULL) +
            theme_bw() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1, size = 50, margin = margin(t = -30)),
                axis.ticks.length = unit(0, "cm"),
                legend.position = "none",
                text = element_text(size = 50),
                panel.background = element_rect(fill = "white", colour = "white"),
                panel.grid.major = element_line(colour = "white"),
                panel.grid.minor = element_line(colour = "white"),
                panel.border = element_blank(),
                plot.margin = margin(t = 10, r = 10, b = 10, l = 80)
            )
        filename <- paste0(folder_workplace, cohort, "_", algorithm, "_mutation_count.png")
        png(filename, res = 150, width = 30, height = 15, units = "in", pointsize = 12)
        print(p)
        dev.off()
        cancer_type_levels <- levels(plot_df$Cancer_type)
    }
    #------------------Plot distribution of cancer-specific sample sizes
    if ("Cancer_type" %in% colnames(results)) {
        plot_df <- results %>%
            group_by(Cancer_type) %>%
            summarise(sample_count = n()) %>%
            arrange(desc(sample_count))
        plot_df$Cancer_type <- factor(plot_df$Cancer_type, levels = cancer_type_levels)
        p <- ggplot(plot_df, aes(x = Cancer_type, y = sample_count, fill = Cancer_type)) +
            geom_bar(stat = "identity", color = "black", size = 2) +
            scale_fill_manual(values = ICGC_cohort_colors) +
            xlab(NULL) +
            ylab("Sample count") +
            labs(NULL) +
            theme_bw() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1, size = 50, margin = margin(t = -30)),
                axis.ticks.length = unit(0, "cm"),
                legend.position = "none",
                text = element_text(size = 50),
                panel.background = element_rect(fill = "white", colour = "white"),
                panel.grid.major = element_line(colour = "white"),
                panel.grid.minor = element_line(colour = "white"),
                panel.border = element_blank(),
                plot.margin = margin(t = 10, r = 10, b = 10, l = 80)
            )
        filename <- paste0(folder_workplace, cohort, "_", algorithm, "_sample_size.png")
        png(filename, res = 150, width = 30, height = 15, units = "in", pointsize = 12)
        print(p)
        dev.off()
    }
    #----------------Plot distribution of cancer-specific tail detection
    if ("Cancer_type" %in% colnames(results)) {
        plot_df <- results %>%
            group_by(Cancer_type) %>%
            summarise(
                total_count = n(),
                tail_count = sum(Tail == TRUE)
            ) %>%
            mutate(proportion = tail_count / total_count)
        plot_df$Cancer_type <- factor(plot_df$Cancer_type, levels = cancer_type_levels)
        p <- ggplot(plot_df, aes(x = Cancer_type, y = proportion, fill = Cancer_type)) +
            geom_bar(stat = "identity", color = "black", size = 2) +
            geom_text(aes(label = paste0("n=", total_count)), vjust = -0.5, size = 10) +
            scale_y_continuous(labels = scales::percent) +
            scale_fill_manual(values = ICGC_cohort_colors) +
            xlab(NULL) +
            ylab("% detected neutral component") +
            labs(NULL) +
            theme_bw() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1, size = 50, margin = margin(t = -30)),
                axis.ticks.length = unit(0, "cm"),
                legend.position = "none",
                text = element_text(size = 50),
                panel.background = element_rect(fill = "white", colour = "white"),
                panel.grid.major = element_line(colour = "white"),
                panel.grid.minor = element_line(colour = "white"),
                panel.border = element_blank(),
                plot.margin = margin(t = 10, r = 10, b = 10, l = 80)
            )
        filename <- paste0(folder_workplace, cohort, "_", algorithm, "_tail_detection.png")
        png(filename, res = 150, width = 30, height = 15, units = "in", pointsize = 12)
        print(p)
        dev.off()
    }
    #---------------------Plot cluster count distribution by cancer type
    if ("Cancer_type" %in% colnames(results)) {
        cluster_distribution_df <- results %>%
            group_by(Cancer_type, Cluster_count) %>%
            summarise(count = n(), .groups = "drop") %>%
            ungroup()

        total_counts_df <- cluster_distribution_df %>%
            group_by(Cancer_type) %>%
            summarise(total_count = sum(count), .groups = "drop")
        cluster_distribution_df <- merge(cluster_distribution_df, total_counts_df, by = "Cancer_type")
        cluster_distribution_df <- cluster_distribution_df %>%
            mutate(percentage = (count / total_count) * 100)
        cluster_distribution_df$Cancer_type <- factor(cluster_distribution_df$Cancer_type, levels = cancer_type_levels)
        cluster_distribution_df$Cluster_count <- factor(cluster_distribution_df$Cluster_count, levels = rev(sort(unique(cluster_distribution_df$Cluster_count))))

        p <- ggplot(cluster_distribution_df, aes(x = Cancer_type, y = percentage, fill = as.factor(Cluster_count))) +
            geom_bar(stat = "identity", position = "stack", width = 1, color = "black", size = 2) +
            xlab(NULL) +
            ylab(NULL) +
            labs(fill = "Cluster count") +
            scale_y_continuous(labels = scales::percent_format(scale = 1)) +
            scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
            scale_fill_manual(values = cluster_colors) +
            guides(fill = guide_legend(keywidth = 2.5, keyheight = 2, reverse = TRUE)) +
            theme_bw() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1, size = 50, margin = margin(t = -30)),
                axis.ticks.length = unit(0, "cm"),
                legend.position = "top",
                legend.justification = "left",
                text = element_text(size = 50),
                panel.background = element_rect(fill = "white", colour = "white"),
                panel.grid.major = element_line(colour = "white"),
                panel.grid.minor = element_line(colour = "white"),
                panel.border = element_blank(),
                plot.margin = margin(t = 0, r = -30, b = 10, l = 95)
            )
        filename <- paste0(folder_workplace, cohort, "_", algorithm, " _cluster_count_distribution_by_cancer_type.png")
        png(filename, res = 150, width = 30, height = 15, units = "in", pointsize = 12)
        print(p)
        dev.off()
    }
    #-------------Plot joint distribution of truncal VAF & sample purity
    if ("Purity" %in% colnames(results)) {
        plot_df <- results %>%
            rowwise() %>%
            mutate(max_vaf = max(c_across(starts_with("Cluster_VAF_")), na.rm = TRUE)) %>%
            # mutate(max_vaf = max(c_across(starts_with("Cluster_VAF_"))[c_across(starts_with("Cluster_VAF_")) < 0.5], na.rm = TRUE)) %>%
            ungroup() %>%
            mutate(max_vaf_scaled = max_vaf * 2) %>%
            mutate(within_bounds = ifelse(max_vaf_scaled >= Purity - 0.1 & max_vaf_scaled <= Purity + 0.1, "Correct", "Wrong"))
        within_bounds_percentage <- 100 * sum(plot_df$within_bounds == "Correct") / nrow(plot_df)
        out_of_bounds_percentage <- 100 * sum(plot_df$within_bounds == "Wrong") / nrow(plot_df)
        p <- ggplot() +
            geom_abline(intercept = 0.1, slope = 1, color = "grey", linewidth = 2) +
            geom_abline(intercept = -0.1, slope = 1, color = "grey", linewidth = 2) +
            geom_point(data = plot_df, aes(x = Purity, y = max_vaf_scaled, color = within_bounds), size = 10, alpha = 0.5, stroke = 2) +
            scale_color_manual(
                values = c("Correct" = "#D55E00", "Wrong" = "#56B4E9"),
                labels = c(
                    paste0("Within bounds (", round(within_bounds_percentage, 1), "%)"),
                    paste0("Out of bounds (", round(out_of_bounds_percentage, 1), "%)")
                )
            ) +
            scale_x_continuous(name = "Sample purity", breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
            scale_y_continuous(name = "Predicted purity", breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
            theme_bw() +
            theme(
                legend.position = "top",
                legend.justification = c(0, 1),
                legend.title = element_blank(),
                plot.title = element_blank(),
                aspect.ratio = 1,
                text = element_text(size = 50),
                panel.background = element_rect(fill = "white", colour = "white"),
                panel.grid.major = element_line(colour = "white"),
                panel.grid.minor = element_line(colour = "white"),
                panel.border = element_blank()
            )
        filename <- paste0(folder_workplace, cohort, "_", algorithm, "_predicted_vs_true_purity.png")
        png(filename, res = 150, width = 15, height = 15, units = "in", pointsize = 12)
        print(p)
        dev.off()
    }
}
